#define _POSIX_C_SOURCE 200809L
#define _FILE_OFFSET_BITS 64

#include "tar_stream.h"
#include <stdlib.h>
#include <string.h>

#define TAR_BLOCK_SIZE   512
#define MAX_ENTRY_BYTES  (256ULL * 1024ULL * 1024ULL)

static bool is_zero_block(const uint8_t blk[TAR_BLOCK_SIZE]) {
    for (size_t i = 0; i < TAR_BLOCK_SIZE; i++) {
        if (blk[i] != 0) return false;
    }
    return true;
}

static uint64_t parse_tar_size(const uint8_t field[12]) {
    if (field[0] & 0x80) {
        // GNU base-256 encoding
        uint64_t val = field[0] & 0x7F;
        for (int i = 1; i < 12; i++) {
            val = (val << 8) | field[i];
        }
        return val;
    }
    uint64_t val = 0;
    for (int i = 0; i < 12; i++) {
        if (field[i] == '\0' || field[i] == ' ') {
            if (val > 0) break;
            continue;
        }
        if (field[i] >= '0' && field[i] <= '7') {
            val = (val << 3) | (uint64_t)(field[i] - '0');
        }
    }
    return val;
}

static bool skip_bytes(TarStream* ts, uint64_t n) {
    if (n == 0) return true;
    if (ts->is_seekable && fseeko(ts->fp, (off_t)n, SEEK_CUR) == 0) {
        return true;
    }
    uint8_t sink[16384];
    while (n > 0) {
        size_t chunk = (n > sizeof(sink)) ? sizeof(sink) : (size_t)n;
        if (fread(sink, 1, chunk, ts->fp) != chunk) return false;
        n -= chunk;
    }
    return true;
}

static bool read_exact(FILE* fp, void* dst, size_t n) {
    return (n == 0) || (fread(dst, 1, n, fp) == n);
}

static bool has_suffix(const char* s, const char* suffix) {
    size_t slen = strlen(s);
    size_t suflen = strlen(suffix);
    return (slen >= suflen) && (memcmp(s + slen - suflen, suffix, suflen) == 0);
}

static bool is_target_cif_name(const char* name) {
    const char* base = strrchr(name, '/');
    base = base ? (base + 1) : name;
    if (base[0] == '\0' || (base[0] == '.' && base[1] == '_')) return false;
    return has_suffix(base, ".cif.gz") || has_suffix(base, ".cif");
}

static void parse_pax_path(const char* buf, size_t len, char* out_name, size_t out_cap, bool* found) {
    size_t pos = 0;
    while (pos < len) {
        size_t rec_len = 0;
        size_t p = pos;
        while (p < len && buf[p] >= '0' && buf[p] <= '9') {
            rec_len = rec_len * 10 + (size_t)(buf[p] - '0');
            p++;
        }
        if (rec_len == 0 || pos + rec_len > len || p >= len || buf[p] != ' ') break;
        p++; // skip space
        size_t kv_len = (pos + rec_len) - p;
        if (kv_len > 6 && memcmp(buf + p, "path=", 5) == 0) {
            size_t vlen = kv_len - 5;
            if (vlen > 0 && buf[p + 5 + vlen - 1] == '\n') vlen--;
            if (vlen >= out_cap) vlen = out_cap - 1;
            memcpy(out_name, buf + p + 5, vlen);
            out_name[vlen] = '\0';
            *found = true;
        }
        pos += rec_len;
    }
}

bool tar_stream_open(TarStream* ts, const char* path) {
    memset(ts, 0, sizeof(*ts));
    if (strcmp(path, "-") == 0) {
        ts->fp = stdin;
        ts->close_on_destroy = false;
        ts->is_seekable = false;
    } else {
        ts->fp = fopen(path, "rb");
        if (!ts->fp) return false;
        ts->close_on_destroy = true;
        ts->is_seekable = (fseeko(ts->fp, 0, SEEK_CUR) == 0);
    }
    return true;
}

void tar_stream_close(TarStream* ts) {
    if (ts->fp && ts->close_on_destroy) {
        fclose(ts->fp);
    }
    memset(ts, 0, sizeof(*ts));
}

bool tar_stream_next_cif(TarStream* ts, TarEntry* out_entry, bool* out_error) {
    *out_error = false;
    uint8_t hdr[TAR_BLOCK_SIZE];

    while (true) {
        size_t got = fread(hdr, 1, TAR_BLOCK_SIZE, ts->fp);
        if (got == 0) return false; // clean EOF
        if (got != TAR_BLOCK_SIZE) { *out_error = true; return false; }

        if (is_zero_block(hdr)) {
            // Consume optional second zero trailer block
            (void)fread(hdr, 1, TAR_BLOCK_SIZE, ts->fp);
            return false;
        }

        uint64_t size = parse_tar_size(hdr + 124);
        uint64_t pad  = (TAR_BLOCK_SIZE - (size % TAR_BLOCK_SIZE)) % TAR_BLOCK_SIZE;
        char typeflag = (char)hdr[156];

        // GNU LongLink filename
        if (typeflag == 'L' && size < sizeof(ts->pending_name)) {
            if (!read_exact(ts->fp, ts->pending_name, (size_t)size) || !skip_bytes(ts, pad)) {
                *out_error = true;
                return false;
            }
            ts->pending_name[size] = '\0';
            ts->has_pending_name = true;
            continue;
        }

        // PAX extended header for next file
        if (typeflag == 'x' && size < 65536) {
            char* pax = malloc((size_t)size + 1);
            if (!pax || !read_exact(ts->fp, pax, (size_t)size) || !skip_bytes(ts, pad)) {
                free(pax);
                *out_error = true;
                return false;
            }
            pax[size] = '\0';
            parse_pax_path(pax, (size_t)size, ts->pending_name, sizeof(ts->pending_name), &ts->has_pending_name);
            free(pax);
            continue;
        }

        // Build full entry path
        char name[4096];
        if (ts->has_pending_name) {
            snprintf(name, sizeof(name), "%s", ts->pending_name);
            ts->has_pending_name = false;
        } else {
            char raw_name[101] = {0};
            char prefix[156]   = {0};
            memcpy(raw_name, hdr, 100);
            memcpy(prefix, hdr + 345, 155);
            if (prefix[0] != '\0') {
                snprintf(name, sizeof(name), "%.155s/%.100s", prefix, raw_name);
            } else {
                snprintf(name, sizeof(name), "%.100s", raw_name);
            }
        }

        bool is_regular = (typeflag == '0' || typeflag == '\0');
        if (!is_regular || size == 0 || size > MAX_ENTRY_BYTES || !is_target_cif_name(name)) {
            ts->skipped_entries++;
            if (!skip_bytes(ts, size + pad)) {
                *out_error = true;
                return false;
            }
            continue;
        }

        uint8_t* data = malloc((size_t)size);
        if (!data) {
            *out_error = true;
            return false;
        }
        if (!read_exact(ts->fp, data, (size_t)size) || !skip_bytes(ts, pad)) {
            free(data);
            *out_error = true;
            return false;
        }

        snprintf(out_entry->name, sizeof(out_entry->name), "%s", name);
        out_entry->data = data;
        out_entry->size = (size_t)size;
        return true;
    }
}
