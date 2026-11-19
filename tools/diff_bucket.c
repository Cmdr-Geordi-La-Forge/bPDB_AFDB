#define _POSIX_C_SOURCE 200809L

#include "diff_bucket.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <stdbool.h>

static void mkdir_if_needed(const char* dir) {
    if (!dir || dir[0] == '\0') return;
    if (mkdir(dir, 0755) != 0 && errno != EEXIST) {
        fprintf(stderr, "Warning: could not create dump dir '%s': %s\n", dir, strerror(errno));
    }
}

void diff_tracker_init(DiffTracker* dt, const char* dump_dir, size_t max_samples) {
    memset(dt, 0, sizeof(*dt));
    pthread_mutex_init(&dt->mu, NULL);
    dt->dump_dir = dump_dir;
    dt->max_samples_per_bucket = max_samples;
    if (dump_dir) mkdir_if_needed(dump_dir);
}

void diff_tracker_destroy(DiffTracker* dt) {
    pthread_mutex_destroy(&dt->mu);
    free(dt->items);
    memset(dt, 0, sizeof(*dt));
}

static size_t find_line_start(const char* buf, size_t pos) {
    while (pos > 0 && buf[pos - 1] != '\n') pos--;
    return pos;
}

static size_t find_line_end(const char* buf, size_t len, size_t pos) {
    while (pos < len && buf[pos] != '\n') pos++;
    return pos;
}

static size_t count_line_number(const char* buf, size_t pos) {
    size_t line = 1;
    for (size_t i = 0; i < pos; i++) {
        if (buf[i] == '\n') line++;
    }
    return line;
}

static void classify_cif_diff(const char* orig, size_t orig_len, size_t diff_pos,
                              char out_cat[64]) {
    if (diff_pos >= orig_len) {
        snprintf(out_cat, 64, "eof_length_mismatch");
        return;
    }

    size_t lstart = find_line_start(orig, diff_pos);
    size_t lend   = find_line_end(orig, orig_len, lstart);
    size_t llen   = lend - lstart;

    if (llen >= 5 && memcmp(orig + lstart, "data_", 5) == 0) {
        snprintf(out_cat, 64, "data_header");
        return;
    }
    if (llen >= 1 && orig[lstart] == '#') {
        snprintf(out_cat, 64, "block_separator");
        return;
    }
    if (llen >= 1 && orig[lstart] == '_') {
        size_t k = 0;
        while (k < llen && k < 63 && orig[lstart + k] != ' ' && orig[lstart + k] != '\t' && orig[lstart + k] != '\r') {
            out_cat[k] = orig[lstart + k];
            k++;
        }
        out_cat[k] = '\0';
        return;
    }

    // Data row or multiline value: walk backward to find enclosing _category
    size_t cur = lstart;
    size_t tag_start = 0, tag_len = 0;
    bool found_tag = false;
    bool in_loop   = false;

    while (cur > 0) {
        size_t prev_end   = cur - 1; // points at '\n'
        size_t prev_start = find_line_start(orig, prev_end);
        size_t plen       = prev_end - prev_start;

        if (!found_tag && plen > 0 && orig[prev_start] == '_') {
            tag_start = prev_start;
            tag_len   = plen;
            found_tag = true;
        } else if (found_tag) {
            if (plen >= 5 && memcmp(orig + prev_start, "loop_", 5) == 0) {
                in_loop = true;
                break;
            }
            if (plen == 0 || orig[prev_start] != '_') break;
        } else if (plen > 0 && orig[prev_start] == '#') {
            break;
        }
        cur = prev_start;
    }

    if (!found_tag) {
        snprintf(out_cat, 64, "unknown_section");
        return;
    }

    if (in_loop) {
        // Extract "_category.*"
        size_t k = 0;
        while (k < tag_len && k < 60 && orig[tag_start + k] != '.' && orig[tag_start + k] != ' ') {
            out_cat[k] = orig[tag_start + k];
            k++;
        }
        memcpy(out_cat + k, ".*", 3);
    } else {
        size_t k = 0;
        while (k < tag_len && k < 63 && orig[tag_start + k] != ' ' && orig[tag_start + k] != '\t' && orig[tag_start + k] != '\r') {
            out_cat[k] = orig[tag_start + k];
            k++;
        }
        out_cat[k] = '\0';
    }
}

static void sanitize_filename(const char* src, char* dst, size_t cap) {
    size_t d = 0;
    for (size_t i = 0; src[i] != '\0' && d + 1 < cap; i++) {
        char c = src[i];
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
            (c >= '0' && c <= '9') || c == '_' || c == '-' || c == '.') {
            dst[d++] = c;
        } else {
            dst[d++] = '_';
        }
    }
    dst[d] = '\0';
}

static const char* entry_basename(const char* name) {
    const char* slash = strrchr(name, '/');
    return slash ? (slash + 1) : name;
}

static DiffBucket* get_or_create_bucket_locked(DiffTracker* dt, const char* category, const char* example) {
    for (size_t i = 0; i < dt->count; i++) {
        if (strcmp(dt->items[i].category, category) == 0) {
            return &dt->items[i];
        }
    }
    if (dt->count >= dt->cap) {
        dt->cap = (dt->cap == 0) ? 16 : dt->cap * 2;
        dt->items = realloc(dt->items, dt->cap * sizeof(DiffBucket));
    }
    DiffBucket* b = &dt->items[dt->count++];
    memset(b, 0, sizeof(*b));
    snprintf(b->category, sizeof(b->category), "%.63s", category);
    snprintf(b->first_example, sizeof(b->first_example), "%.127s", example);
    return b;
}

static void write_file_bytes(const char* path, const char* data, size_t len) {
    FILE* f = fopen(path, "wb");
    if (!f) return;
    if (len > 0) fwrite(data, 1, len, f);
    fclose(f);
}

void diff_tracker_record_error(DiffTracker* dt, const char* entry_name,
                               const char* category,
                               const char* orig, size_t orig_len) {
    const char* base = entry_basename(entry_name);
    bool should_dump = false;

    pthread_mutex_lock(&dt->mu);
    DiffBucket* b = get_or_create_bucket_locked(dt, category, base);
    b->count++;
    if (dt->dump_dir && b->saved_samples < dt->max_samples_per_bucket) {
        b->saved_samples++;
        should_dump = true;
    }
    pthread_mutex_unlock(&dt->mu);

    if (should_dump && orig && orig_len > 0) {
        char safe_cat[64], safe_base[96], path[512];
        sanitize_filename(category, safe_cat, sizeof(safe_cat));
        sanitize_filename(base, safe_base, sizeof(safe_base));
        snprintf(path, sizeof(path), "%.200s/%.63s__%.95s.orig.cif", dt->dump_dir, safe_cat, safe_base);
        write_file_bytes(path, orig, orig_len);
    }
}

void diff_tracker_record_mismatch(DiffTracker* dt, const char* entry_name,
                                  const char* orig, size_t orig_len,
                                  const char* dec, size_t dec_len) {
    size_t min_len = (orig_len < dec_len) ? orig_len : dec_len;
    size_t pos = 0;
    while (pos < min_len && orig[pos] == dec[pos]) pos++;

    char category[64];
    classify_cif_diff(orig, orig_len, pos, category);

    const char* base = entry_basename(entry_name);
    bool should_dump = false;

    pthread_mutex_lock(&dt->mu);
    DiffBucket* b = get_or_create_bucket_locked(dt, category, base);
    b->count++;
    if (dt->dump_dir && b->saved_samples < dt->max_samples_per_bucket) {
        b->saved_samples++;
        should_dump = true;
    }
    pthread_mutex_unlock(&dt->mu);

    if (!should_dump) return;

    char safe_cat[64], safe_base[96], path[512];
    sanitize_filename(category, safe_cat, sizeof(safe_cat));
    sanitize_filename(base, safe_base, sizeof(safe_base));

    snprintf(path, sizeof(path), "%.200s/%.63s__%.95s.orig.cif", dt->dump_dir, safe_cat, safe_base);
    write_file_bytes(path, orig, orig_len);

    snprintf(path, sizeof(path), "%.200s/%.63s__%.95s.dec.cif", dt->dump_dir, safe_cat, safe_base);
    write_file_bytes(path, dec, dec_len);

    snprintf(path, sizeof(path), "%.200s/%.63s__%.95s.diff.txt", dt->dump_dir, safe_cat, safe_base);
    FILE* f = fopen(path, "w");
    if (f) {
        size_t o_start = find_line_start(orig, (pos < orig_len) ? pos : orig_len);
        size_t o_end   = find_line_end(orig, orig_len, o_start);
        size_t d_start = find_line_start(dec, (pos < dec_len) ? pos : dec_len);
        size_t d_end   = find_line_end(dec, dec_len, d_start);

        fprintf(f, "entry:    %s\n", entry_name);
        fprintf(f, "category: %s\n", category);
        fprintf(f, "byte_off: %zu (line %zu)\n", pos, count_line_number(orig, (pos < orig_len) ? pos : orig_len));
        fprintf(f, "orig_len: %zu, dec_len: %zu\n\n", orig_len, dec_len);
        fprintf(f, "ORIG: %.*s\n", (int)(o_end - o_start), orig + o_start);
        fprintf(f, "DEC : %.*s\n", (int)(d_end - d_start), dec + d_start);
        fclose(f);
    }
}

static int cmp_buckets_desc(const void* a, const void* b) {
    const DiffBucket* ba = (const DiffBucket*)a;
    const DiffBucket* bb = (const DiffBucket*)b;
    if (bb->count > ba->count) return 1;
    if (bb->count < ba->count) return -1;
    return strcmp(ba->category, bb->category);
}

void diff_tracker_print_report(DiffTracker* dt) {
    pthread_mutex_lock(&dt->mu);
    if (dt->count == 0) {
        pthread_mutex_unlock(&dt->mu);
        return;
    }
    qsort(dt->items, dt->count, sizeof(DiffBucket), cmp_buckets_desc);

    printf("\nMismatch Breakdown (%zu distinct buckets):\n", dt->count);
    printf("  %-38s  %10s  %s\n", "CATEGORY", "COUNT", "FIRST EXAMPLE");
    printf("  %-38s  %10s  %s\n", "--------------------------------------", "----------", "------------------------------");
    for (size_t i = 0; i < dt->count; i++) {
        printf("  %-38s  %10llu  %s\n",
               dt->items[i].category,
               (unsigned long long)dt->items[i].count,
               dt->items[i].first_example);
    }
    pthread_mutex_unlock(&dt->mu);
}
