#define _POSIX_C_SOURCE 200809L
#define _FILE_OFFSET_BITS 64

#include "bpdb.h"
#include "bpdb_internal.h"
#include "pack_format.h"

#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#define MERGE_FANIN 256U

typedef struct {
    char**   paths;
    uint32_t count;
    uint32_t cap;
} PathList;

static bool has_suffix(const char* s, const char* suf) {
    size_t sl = strlen(s), sufl = strlen(suf);
    return (sl >= sufl) && (memcmp(s + sl - sufl, suf, sufl) == 0);
}

static bool pathlist_push(PathList* pl, const char* rel) {
    if (pl->count >= pl->cap) {
        uint32_t ncap = (pl->cap == 0) ? 256U : (pl->cap * 2U);
        char** np = realloc(pl->paths, (size_t)ncap * sizeof(char*));
        if (!np) return false;
        pl->paths = np;
        pl->cap   = ncap;
    }
    char* dup = strdup(rel);
    if (!dup) return false;
    pl->paths[pl->count++] = dup;
    return true;
}

static void pathlist_free(PathList* pl) {
    for (uint32_t i = 0; i < pl->count; i++) free(pl->paths[i]);
    free(pl->paths);
    memset(pl, 0, sizeof(*pl));
}

static bool walk_packs(const char* root_dir, const char* rel_sub, PathList* out) {
    char full_dir[1024];
    if (rel_sub && rel_sub[0] != '\0') snprintf(full_dir, sizeof(full_dir), "%s/%s", root_dir, rel_sub);
    else                               snprintf(full_dir, sizeof(full_dir), "%s", root_dir);

    DIR* d = opendir(full_dir);
    if (!d) return false;

    struct dirent* ent;
    while ((ent = readdir(d)) != NULL) {
        if (ent->d_name[0] == '.') continue;
        char child_rel[1024];
        char child_full[2048];
        if (rel_sub && rel_sub[0] != '\0') snprintf(child_rel, sizeof(child_rel), "%s/%s", rel_sub, ent->d_name);
        else                               snprintf(child_rel, sizeof(child_rel), "%s", ent->d_name);
        snprintf(child_full, sizeof(child_full), "%s/%s", root_dir, child_rel);

        struct stat st;
        if (stat(child_full, &st) != 0) continue;
        if (S_ISDIR(st.st_mode)) walk_packs(root_dir, child_rel, out);
        else if (S_ISREG(st.st_mode) && has_suffix(ent->d_name, ".bpdbpack")) pathlist_push(out, child_rel);
    }
    closedir(d);
    return true;
}

typedef struct {
    FILE*              fp;
    uint32_t           shard_id;
    uint32_t           priority; // 0 = v6 overlay, 1 = v4 base
    bool               is_raw_run;
    uint64_t           remaining;
    BpdbGlobalIdxEntry cur;
    bool               valid;
} MergeStream;

typedef struct {
    MergeStream* streams;
    size_t       stream_count;
    size_t*      heap;
    size_t       heap_size;
} MergeHeap;

static int cmp_merge_items(const MergeStream* a, const MergeStream* b) {
    int c = bpdb_acc_frag_cmp(a->cur.accession, a->cur.frag_num, b->cur.accession, b->cur.frag_num);
    if (c != 0) return c;
    if (a->priority != b->priority) return (a->priority < b->priority) ? -1 : 1;
    return (a->cur.shard_id < b->cur.shard_id) ? -1 : (a->cur.shard_id > b->cur.shard_id);
}

static bool merge_stream_advance(MergeStream* ms) {
    if (ms->remaining == 0) { ms->valid = false; return false; }
    if (ms->is_raw_run) {
        if (fread(&ms->cur, sizeof(BpdbGlobalIdxEntry), 1, ms->fp) != 1) { ms->valid = false; return false; }
        ms->priority = 0;
    } else {
        BpdbPackTocEntry te;
        if (fread(&te, sizeof(te), 1, ms->fp) != 1) { ms->valid = false; return false; }
        memcpy(ms->cur.accession, te.accession, 12);
        ms->cur.frag_num     = te.frag_num;
        ms->cur.frag_total   = te.frag_total;
        ms->cur.shard_id     = ms->shard_id;
        ms->cur.offset_div16 = te.offset_div16;
        ms->cur.byte_len     = te.byte_len;
    }
    ms->remaining--;
    ms->valid = true;
    return true;
}

static void heap_sift_up(MergeHeap* h, size_t idx) {
    while (idx > 0) {
        size_t par = (idx - 1) / 2;
        if (cmp_merge_items(&h->streams[h->heap[idx]], &h->streams[h->heap[par]]) < 0) {
            size_t tmp = h->heap[idx]; h->heap[idx] = h->heap[par]; h->heap[par] = tmp;
            idx = par;
        } else break;
    }
}

static void heap_sift_down(MergeHeap* h, size_t idx) {
    while (true) {
        size_t l = 2 * idx + 1, r = 2 * idx + 2, best = idx;
        if (l < h->heap_size && cmp_merge_items(&h->streams[h->heap[l]], &h->streams[h->heap[best]]) < 0) best = l;
        if (r < h->heap_size && cmp_merge_items(&h->streams[h->heap[r]], &h->streams[h->heap[best]]) < 0) best = r;
        if (best == idx) break;
        size_t tmp = h->heap[idx]; h->heap[idx] = h->heap[best]; h->heap[best] = tmp;
        idx = best;
    }
}

static bool merge_streams_to_fp(MergeStream* streams, size_t count, FILE* out_fp, uint64_t* out_written) {
    MergeHeap h = {
        .streams = streams, .stream_count = count,
        .heap = malloc((count > 0 ? count : 1) * sizeof(size_t)), .heap_size = 0
    };
    if (!h.heap) return false;

    for (size_t i = 0; i < count; i++) {
        if (merge_stream_advance(&streams[i])) {
            h.heap[h.heap_size] = i;
            heap_sift_up(&h, h.heap_size);
            h.heap_size++;
        }
    }

    uint64_t written = 0;
    bool have_prev = false;
    char prev_acc[12] = {0};
    uint16_t prev_frag = 0;
    uint16_t prev_frag_total = 0;

    while (h.heap_size > 0) {
        size_t top_s = h.heap[0];
        BpdbGlobalIdxEntry item = streams[top_s].cur;

        if (merge_stream_advance(&streams[top_s])) {
            heap_sift_down(&h, 0);
        } else {
            h.heap[0] = h.heap[--h.heap_size];
            if (h.heap_size > 0) heap_sift_down(&h, 0);
        }

        // Deduplicate: earlier shard_id (v6 overlays) wins; drop stale trailing frags
        if (have_prev && memcmp(prev_acc, item.accession, 12) == 0) {
            if (item.frag_num == prev_frag || item.frag_num > prev_frag_total) continue;
            prev_frag = item.frag_num;
        } else {
            memcpy(prev_acc, item.accession, 12);
            prev_frag       = item.frag_num;
            prev_frag_total = (item.frag_total > 0) ? item.frag_total : 1;
            have_prev       = true;
        }

        if (fwrite(&item, sizeof(item), 1, out_fp) != 1) { free(h.heap); return false; }
        written++;
    }

    free(h.heap);
    if (out_written) *out_written = written;
    return true;
}

typedef struct {
    char*    rel_path;
    uint32_t priority;
} RankedPack;

static int cmp_ranked_pack(const void* a, const void* b) {
    const RankedPack* pa = (const RankedPack*)a;
    const RankedPack* pb = (const RankedPack*)b;
    if (pa->priority != pb->priority) return (pa->priority < pb->priority) ? -1 : 1;
    return strcmp(pa->rel_path, pb->rel_path);
}

static int cmd_build(int argc, char** argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: bpdb_index build <db_dir> <out.idx>\n");
        return 2;
    }
    const char* db_dir   = argv[0];
    const char* out_path = argv[1];

    PathList pl = {0};
    if (!walk_packs(db_dir, "", &pl) || pl.count == 0) {
        fprintf(stderr, "Error: no .bpdbpack files found under '%s'\n", db_dir);
        pathlist_free(&pl);
        return 1;
    }

    RankedPack* ranked = calloc(pl.count, sizeof(RankedPack));
    uint32_t valid_count = 0;
    for (uint32_t i = 0; i < pl.count; i++) {
        char full[1024];
        snprintf(full, sizeof(full), "%s/%s", db_dir, pl.paths[i]);
        FILE* f = fopen(full, "rb");
        if (!f) continue;
        BpdbPackHeader hdr;
        bool ok = (fread(&hdr, 1, sizeof(hdr), f) == sizeof(hdr));
        fclose(f);
        if (!ok || memcmp(hdr.magic, BPDB_PACK_MAGIC, 4) != 0 || hdr.version != BPDB_PACK_VERSION) {
            fprintf(stderr, "Warning: skipping invalid pack '%s'\n", pl.paths[i]);
            continue;
        }
        ranked[valid_count].rel_path = pl.paths[i];
        ranked[valid_count].priority = (hdr.flags & BPDB_PACK_FLAG_AFDB_V6) ? 0U : 1U;
        valid_count++;
    }
    if (valid_count == 0) {
        fprintf(stderr, "Error: zero valid .bpdbpack files found.\n");
        free(ranked); pathlist_free(&pl);
        return 1;
    }
    qsort(ranked, valid_count, sizeof(RankedPack), cmp_ranked_pack);
    printf("Discovered %u valid .bpdbpack files under '%s'...\n", valid_count, db_dir);

    size_t num_runs = ((size_t)valid_count + MERGE_FANIN - 1U) / MERGE_FANIN;
    char (*run_paths)[512] = calloc(num_runs, sizeof(*run_paths));
    uint64_t* run_counts   = calloc(num_runs, sizeof(uint64_t));

    char tmp_idx[512];
    snprintf(tmp_idx, sizeof(tmp_idx), "%s.tmp", out_path);

    bool ok = true;
    for (size_t r = 0; r < num_runs && ok; r++) {
        size_t start = r * MERGE_FANIN;
        size_t end   = start + MERGE_FANIN;
        if (end > valid_count) end = valid_count;
        size_t k = end - start;

        MergeStream* ms = calloc(k, sizeof(MergeStream));
        for (size_t i = 0; i < k; i++) {
            char full[1024];
            snprintf(full, sizeof(full), "%s/%s", db_dir, ranked[start + i].rel_path);
            FILE* f = fopen(full, "rb");
            if (!f) { ok = false; break; }
            BpdbPackHeader hdr;
            if (fread(&hdr, 1, sizeof(hdr), f) != sizeof(hdr) ||
                fseeko(f, (off_t)hdr.toc_offset, SEEK_SET) != 0) {
                fclose(f); ok = false; break;
            }
            ms[i].fp         = f;
            ms[i].shard_id   = (uint32_t)(start + i);
            ms[i].priority   = ranked[start + i].priority;
            ms[i].is_raw_run = false;
            ms[i].remaining  = hdr.record_count;
        }

        if (ok) {
            snprintf(run_paths[r], sizeof(run_paths[r]), "%s.run%zu.tmp", out_path, r);
            FILE* rfp = fopen(run_paths[r], "wb");
            if (!rfp) ok = false;
            else {
                setvbuf(rfp, NULL, _IOFBF, 1024 * 1024);
                ok = merge_streams_to_fp(ms, k, rfp, &run_counts[r]);
                fclose(rfp);
            }
        }
        for (size_t i = 0; i < k; i++) if (ms[i].fp) fclose(ms[i].fp);
        free(ms);
    }

    uint64_t final_entries = 0;
    if (ok) {
        FILE* out_fp = fopen(tmp_idx, "wb");
        if (!out_fp) {
            ok = false;
        } else {
            setvbuf(out_fp, NULL, _IOFBF, 1024 * 1024);
            BpdbGlobalIdxHeader hdr;
            memset(&hdr, 0, sizeof(hdr));
            if (fwrite(&hdr, 1, sizeof(hdr), out_fp) != sizeof(hdr)) ok = false;

            if (ok) {
                MergeStream* ms = calloc(num_runs, sizeof(MergeStream));
                for (size_t r = 0; r < num_runs; r++) {
                    ms[r].fp         = fopen(run_paths[r], "rb");
                    ms[r].is_raw_run = true;
                    ms[r].remaining  = run_counts[r];
                    if (!ms[r].fp) { ok = false; break; }
                }
                if (ok) ok = merge_streams_to_fp(ms, num_runs, out_fp, &final_entries);
                for (size_t r = 0; r < num_runs; r++) if (ms[r].fp) fclose(ms[r].fp);
                free(ms);
            }

            if (ok) {
                uint64_t shard_tbl_off = sizeof(BpdbGlobalIdxHeader) + final_entries * sizeof(BpdbGlobalIdxEntry);
                for (uint32_t i = 0; i < valid_count; i++) {
                    size_t slen = strlen(ranked[i].rel_path) + 1;
                    if (fwrite(ranked[i].rel_path, 1, slen, out_fp) != slen) { ok = false; break; }
                }
                if (ok) {
                    memcpy(hdr.magic, BPDB_IDX_MAGIC, 4);
                    hdr.version            = BPDB_IDX_VERSION;
                    hdr.entry_count        = final_entries;
                    hdr.shard_count        = valid_count;
                    hdr.shard_table_offset = shard_tbl_off;
                    if (fseeko(out_fp, 0, SEEK_SET) != 0 ||
                        fwrite(&hdr, 1, sizeof(hdr), out_fp) != sizeof(hdr) ||
                        fflush(out_fp) != 0) ok = false;
                }
            }
            fclose(out_fp);
        }
    }

    for (size_t r = 0; r < num_runs; r++) if (run_paths[r][0] != '\0') unlink(run_paths[r]);
    free(run_paths);
    free(run_counts);

    if (ok && rename(tmp_idx, out_path) == 0) {
        printf("Built '%s': %llu entries across %u packs.\n",
               out_path, (unsigned long long)final_entries, valid_count);
    } else {
        unlink(tmp_idx);
        fprintf(stderr, "Error: failed building index '%s'\n", out_path);
        ok = false;
    }

    free(ranked);
    pathlist_free(&pl);
    return ok ? 0 : 1;
}

static const char* find_shard_path(const uint8_t* base, size_t fsize,
                                   uint64_t tbl_off, uint32_t shard_id) {
    size_t pos = (size_t)tbl_off;
    for (uint32_t i = 0; i < shard_id; i++) {
        const uint8_t* nul = (const uint8_t*)memchr(base + pos, 0, fsize - pos);
        if (!nul) return NULL;
        pos = (size_t)(nul - base) + 1;
        if (pos >= fsize) return NULL;
    }
    return (const char*)(base + pos);
}

static int cmd_get(int argc, char** argv) {
    if (argc < 3 || argc > 4) {
        fprintf(stderr, "Usage: bpdb_index get <afdb.idx> <db_dir> <accession> [out_prefix|-]\n");
        return 2;
    }
    const char* idx_path   = argv[0];
    const char* db_dir     = argv[1];
    const char* acc_query  = argv[2];
    const char* out_prefix = (argc == 4) ? argv[3] : NULL;

    char acc_key[12] = {0};
    size_t alen = strlen(acc_query);
    if (alen > 12) alen = 12;
    memcpy(acc_key, acc_query, alen);

    int fd = open(idx_path, O_RDONLY);
    if (fd < 0) { fprintf(stderr, "Error: cannot open '%s': %s\n", idx_path, strerror(errno)); return 1; }
    struct stat st;
    if (fstat(fd, &st) != 0 || (size_t)st.st_size < sizeof(BpdbGlobalIdxHeader)) { close(fd); return 1; }
    size_t fsize = (size_t)st.st_size;
    const uint8_t* map = mmap(NULL, fsize, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (map == MAP_FAILED) return 1;

    BpdbGlobalIdxHeader hdr;
    memcpy(&hdr, map, sizeof(hdr));
    if (memcmp(hdr.magic, BPDB_IDX_MAGIC, 4) != 0 || hdr.version != BPDB_IDX_VERSION) {
        fprintf(stderr, "Error: invalid bIX2 index header.\n");
        munmap((void*)map, fsize);
        return 1;
    }

    const BpdbGlobalIdxEntry* entries = (const BpdbGlobalIdxEntry*)(map + sizeof(BpdbGlobalIdxHeader));
    uint64_t lo = 0, hi = hdr.entry_count;
    while (lo < hi) {
        uint64_t mid = lo + (hi - lo) / 2;
        if (memcmp(entries[mid].accession, acc_key, 12) < 0) lo = mid + 1;
        else                                                 hi = mid;
    }
    if (lo >= hdr.entry_count || memcmp(entries[lo].accession, acc_key, 12) != 0) {
        fprintf(stderr, "Accession '%s' not found in '%s'.\n", acc_query, idx_path);
        munmap((void*)map, fsize);
        return 1;
    }

    uint64_t cur = lo;
    while (cur < hdr.entry_count && memcmp(entries[cur].accession, acc_key, 12) == 0) {
        const BpdbGlobalIdxEntry* e = &entries[cur];
        const char* rel = find_shard_path(map, fsize, hdr.shard_table_offset, e->shard_id);
        if (!rel) break;

        char full[1024];
        snprintf(full, sizeof(full), "%s/%s", db_dir, rel);
        FILE* pf = fopen(full, "rb");
        if (!pf) { fprintf(stderr, "Error: cannot open shard '%s'\n", full); munmap((void*)map, fsize); return 1; }

        uint64_t off = (uint64_t)e->offset_div16 << 4;
        uint8_t* rec = malloc(e->byte_len);
        bool rd_ok = rec && fseeko(pf, (off_t)off, SEEK_SET) == 0 && fread(rec, 1, e->byte_len, pf) == e->byte_len;
        fclose(pf);
        if (!rd_ok) { free(rec); munmap((void*)map, fsize); return 1; }

        BpdbHotView hv;
        if (bpdb_open_hot(rec, e->byte_len, &hv)) {
            fprintf(stderr, "AF-%s-F%u (frag %u/%u) | res=%u atoms=%u chain=%c | %u B @ %s+%llu\n",
                    hv.accession, (unsigned)e->frag_num, (unsigned)e->frag_num, (unsigned)e->frag_total,
                    hv.total_res, hv.total_atoms, hv.chain_id, e->byte_len, rel, (unsigned long long)off);
        }

        if (out_prefix) {
            if (strcmp(out_prefix, "-") == 0) {
                decode_binary_to_afdb(rec, e->byte_len, stdout);
            } else {
                char out_cif[512];
                snprintf(out_cif, sizeof(out_cif), "%s-F%u.cif", out_prefix, (unsigned)e->frag_num);
                FILE* out_f = fopen(out_cif, "w");
                if (out_f) {
                    decode_binary_to_afdb(rec, e->byte_len, out_f);
                    fclose(out_f);
                    fprintf(stderr, "  -> wrote %s\n", out_cif);
                }
            }
        }
        free(rec);
        cur++;
    }

    munmap((void*)map, fsize);
    return 0;
}

static int cmd_info(int argc, char** argv) {
    if (argc != 1) {
        fprintf(stderr, "Usage: bpdb_index info <file.bpdbpack>\n");
        return 2;
    }
    FILE* f = fopen(argv[0], "rb");
    if (!f) { fprintf(stderr, "Error: cannot open '%s'\n", argv[0]); return 1; }

    BpdbPackHeader hdr;
    if (fread(&hdr, 1, sizeof(hdr), f) != sizeof(hdr) ||
        memcmp(hdr.magic, BPDB_PACK_MAGIC, 4) != 0 || hdr.version != BPDB_PACK_VERSION) {
        fprintf(stderr, "Error: invalid bPK2 pack file '%s'\n", argv[0]);
        fclose(f);
        return 1;
    }

    uint64_t multi_frags = 0, united_ok = 0, split_err = 0, order_err = 0;
    if (fseeko(f, (off_t)hdr.toc_offset, SEEK_SET) == 0 && hdr.record_count > 0) {
        BpdbPackTocEntry* toc = malloc((size_t)hdr.record_count * sizeof(BpdbPackTocEntry));
        if (toc && fread(toc, sizeof(BpdbPackTocEntry), hdr.record_count, f) == hdr.record_count) {
            for (uint32_t i = 0; i < hdr.record_count; i++) {
                if (i > 0 && bpdb_acc_frag_cmp(toc[i - 1].accession, toc[i - 1].frag_num,
                                               toc[i].accession, toc[i].frag_num) >= 0) {
                    order_err++;
                }
                if (toc[i].frag_total > 1 && toc[i].frag_num == 1) {
                    multi_frags++;
                    bool contiguous = true;
                    for (uint16_t k = 1; k < toc[i].frag_total; k++) {
                        if (i + k >= hdr.record_count ||
                            memcmp(toc[i + k].accession, toc[i].accession, 12) != 0 ||
                            toc[i + k].frag_num != (uint16_t)(k + 1)) {
                            contiguous = false;
                            break;
                        }
                        uint64_t prev_end = ((uint64_t)toc[i + k - 1].offset_div16 << 4) +
                                            ((toc[i + k - 1].byte_len + 15U) & ~15U);
                        if (((uint64_t)toc[i + k].offset_div16 << 4) != prev_end) contiguous = false;
                    }
                    if (contiguous) united_ok++;
                    else            split_err++;
                }
            }
        }
        free(toc);
    }
    fclose(f);

    printf("Pack:        %s\n", argv[0]);
    printf("Header:      bPK2 v%u | flags=0x%04x | tax_id=%u | part=%u/%u\n",
           hdr.version, hdr.flags, hdr.tax_id, (unsigned)(hdr.part_idx + 1), (unsigned)hdr.part_count);
    printf("Records:     %u (payload=%.2f MiB, toc_offset=%llu)\n",
           hdr.record_count, (double)hdr.payload_bytes / (1024.0 * 1024.0),
           (unsigned long long)hdr.toc_offset);
    printf("Fragments:   %llu multi-fragment proteins (%llu physically contiguous, %llu errors)\n",
           (unsigned long long)multi_frags, (unsigned long long)united_ok, (unsigned long long)split_err);
    printf("TOC Order:   %s\n", (order_err == 0) ? "OK (strictly sorted)" : "ERROR");
    return (order_err == 0 && split_err == 0) ? 0 : 1;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr,
            "Usage:\n"
            "  %s build <db_dir> <out.idx>\n"
            "  %s get   <afdb.idx> <db_dir> <accession> [out_prefix|-]\n"
            "  %s info  <file.bpdbpack>\n",
            argv[0], argv[0], argv[0]);
        return 2;
    }
    if (strcmp(argv[1], "build") == 0) return cmd_build(argc - 2, argv + 2);
    if (strcmp(argv[1], "get") == 0)   return cmd_get(argc - 2, argv + 2);
    if (strcmp(argv[1], "info") == 0)  return cmd_info(argc - 2, argv + 2);
    fprintf(stderr, "Unknown subcommand: %s\n", argv[1]);
    return 2;
}
