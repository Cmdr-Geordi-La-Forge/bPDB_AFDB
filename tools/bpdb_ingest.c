#define _POSIX_C_SOURCE 200809L
#define _FILE_OFFSET_BITS 64

#include "bpdb.h"
#include "gz_inflate.h"
#include "pack_format.h"
#include "pack_writer.h"
#include "gcs_stream.h"

#include <curl/curl.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define MIB (1024ULL * 1024ULL)
#define GIB (1024ULL * MIB)
#define TIB (1024ULL * GIB)

#define PROGRESS_FLUSH_EVERY 256U

typedef struct {
    uint32_t tax_id;
    uint32_t shard_idx;
    uint64_t tar_bytes;
    char     filename[64];
} ManifestShard;

typedef struct {
    uint32_t tax_id;
    uint64_t total_tar_bytes;
    size_t   shard_start;
    size_t   shard_count;
} TaxGroup;

typedef struct {
    uint32_t batch_id;
    size_t   tax_start;
    size_t   tax_count;
    size_t   shard_count;
    uint64_t total_tar_bytes;
} SmallBatch;

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void fmt_bytes(uint64_t n, char out[32]) {
    if (n >= TIB)      snprintf(out, 32, "%.2f TiB", (double)n / (double)TIB);
    else if (n >= GIB) snprintf(out, 32, "%.2f GiB", (double)n / (double)GIB);
    else if (n >= MIB) snprintf(out, 32, "%.2f MiB", (double)n / (double)MIB);
    else               snprintf(out, 32, "%.1f KiB", (double)n / 1024.0);
}

static void fmt_eta(double sec, char out[32]) {
    if (sec <= 0.0 || sec > 365.0 * 86400.0) {
        snprintf(out, 32, "--");
        return;
    }
    uint64_t s = (uint64_t)sec;
    uint64_t d = s / 86400ULL; s %= 86400ULL;
    uint64_t h = s / 3600ULL;  s %= 3600ULL;
    uint64_t m = s / 60ULL;    s %= 60ULL;
    if (d > 0)      snprintf(out, 32, "%llud%02lluh", (unsigned long long)d, (unsigned long long)h);
    else if (h > 0) snprintf(out, 32, "%lluh%02llum", (unsigned long long)h, (unsigned long long)m);
    else if (m > 0) snprintf(out, 32, "%llum%02llus", (unsigned long long)m, (unsigned long long)s);
    else            snprintf(out, 32, "%llus", (unsigned long long)s);
}

static int cmp_manifest_shard(const void* a, const void* b) {
    const ManifestShard* sa = (const ManifestShard*)a;
    const ManifestShard* sb = (const ManifestShard*)b;
    if (sa->tax_id != sb->tax_id) return (sa->tax_id < sb->tax_id) ? -1 : 1;
    if (sa->shard_idx != sb->shard_idx) return (sa->shard_idx < sb->shard_idx) ? -1 : 1;
    return 0;
}

static int cmp_tax_desc(const void* a, const void* b) {
    const TaxGroup* ta = (const TaxGroup*)a;
    const TaxGroup* tb = (const TaxGroup*)b;
    if (tb->total_tar_bytes != ta->total_tar_bytes)
        return (tb->total_tar_bytes > ta->total_tar_bytes) ? 1 : -1;
    return (ta->tax_id < tb->tax_id) ? -1 : (ta->tax_id > tb->tax_id);
}

static bool parse_manifest_v4(const char* path, ManifestShard** out_shards, size_t* out_count) {
    FILE* fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "Error: cannot open manifest '%s'\n", path);
        return false;
    }

    size_t cap = 65536, count = 0;
    ManifestShard* arr = malloc(cap * sizeof(ManifestShard));
    if (!arr) { fclose(fp); return false; }

    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        const char* p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p < '0' || *p > '9') continue;

        char* endp = NULL;
        uint64_t sz = strtoull(p, &endp, 10);
        if (!endp || (*endp != ' ' && *endp != '\t')) continue;

        const char* marker = strstr(endp, "proteome-tax_id-");
        if (!marker) continue;

        char fname[64];
        size_t flen = 0;
        while (marker[flen] && marker[flen] != '\n' && marker[flen] != '\r' &&
               marker[flen] != ' ' && marker[flen] != '\t' && flen + 1 < sizeof(fname)) {
            fname[flen] = marker[flen];
            flen++;
        }
        fname[flen] = '\0';
        if (flen < 7 || memcmp(fname + flen - 7, "_v4.tar", 7) != 0) continue;

        const char* num = fname + 16;
        char* dash = NULL;
        unsigned long tax = strtoul(num, &dash, 10);
        if (!dash || *dash != '-') continue;
        char* under = NULL;
        unsigned long sh = strtoul(dash + 1, &under, 10);
        if (!under || *under != '_') continue;

        if (count >= cap) {
            cap *= 2;
            ManifestShard* na = realloc(arr, cap * sizeof(ManifestShard));
            if (!na) { free(arr); fclose(fp); return false; }
            arr = na;
        }
        arr[count].tax_id    = (uint32_t)tax;
        arr[count].shard_idx = (uint32_t)sh;
        arr[count].tar_bytes = sz;
        memcpy(arr[count].filename, fname, flen + 1);
        count++;
    }
    fclose(fp);

    if (count > 1) qsort(arr, count, sizeof(ManifestShard), cmp_manifest_shard);
    *out_shards = arr;
    *out_count  = count;
    return true;
}

/* ===================== Live Progress & Verbose Logging ===================== */

typedef struct {
    pthread_mutex_t mu;
    uint64_t        packs_done;
    uint64_t        packs_skipped;
    uint64_t        packs_failed;
    uint64_t        packs_total;
    uint64_t        tars_done;
    uint64_t        est_total_bpdb;
    uint64_t        structures;
    uint64_t        verify_ok;
    uint64_t        verify_fail;
    uint64_t        raw_bytes;
    uint64_t        bpdb_bytes;
    double          start_time;
    double          last_print;
    bool            verbose;
    bool            progress_on_line;
} IngestStats;

static void log_line(IngestStats* s, const char* fmt, ...) {
    va_list ap;
    pthread_mutex_lock(&s->mu);
    if (s->progress_on_line) {
        fprintf(stderr, "\r\033[K");
        s->progress_on_line = false;
    }
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    fflush(stderr);
    pthread_mutex_unlock(&s->mu);
}

typedef struct {
    GzInflater*     inf;
    BpdbEncoderCtx* enc;
    BpdbPackWriter* pw;
    IngestStats*    stats;
    const char*     phase_label;
    uint64_t        verify_every;
    uint64_t        verify_seq;
    uint64_t        local_count;
    uint64_t        local_raw;
    uint64_t        local_bpdb;
    uint64_t        local_vok;
    uint64_t        local_vfail;
    uint64_t        pack_structs;
    uint64_t        pack_raw;
    uint64_t        pack_bpdb;
} CifEncodeCtx;

static void stats_flush_local(IngestStats* s, CifEncodeCtx* c, bool force_print) {
    pthread_mutex_lock(&s->mu);
    s->structures  += c->local_count;
    s->raw_bytes   += c->local_raw;
    s->bpdb_bytes  += c->local_bpdb;
    s->verify_ok   += c->local_vok;
    s->verify_fail += c->local_vfail;
    c->local_count = c->local_raw = c->local_bpdb = c->local_vok = c->local_vfail = 0;

    double t = now_sec();
    if (force_print || (t - s->last_print) >= 1.0) {
        double dt = t - s->start_time;
        double rate = (dt > 0.0) ? (double)s->structures / dt : 0.0;
        double raw_mb_s = (dt > 0.0) ? ((double)s->raw_bytes / (double)MIB) / dt : 0.0;
        double ratio = (s->bpdb_bytes > 0) ? ((double)s->raw_bytes / (double)s->bpdb_bytes) : 0.0;
        char bpdb_str[32], raw_str[32];
        fmt_bytes(s->bpdb_bytes, bpdb_str);
        fmt_bytes(s->raw_bytes, raw_str);

        if (s->packs_total > 0) {
            char eta_str[32] = "--";
            double pct = 0.0;
            if (s->est_total_bpdb > 0 && s->bpdb_bytes > 0) {
                pct = 100.0 * (double)s->bpdb_bytes / (double)s->est_total_bpdb;
                if (pct > 99.9) pct = 99.9;
                if (dt >= 5.0) {
                    double bps = (double)s->bpdb_bytes / dt;
                    double rem = (s->est_total_bpdb > s->bpdb_bytes)
                               ? (double)(s->est_total_bpdb - s->bpdb_bytes) / bps : 0.0;
                    fmt_eta(rem, eta_str);
                }
            }
            fprintf(stderr,
                    "\r\033[K[%s] %4.1f%% | packs=%llu/%llu (skip=%llu fail=%llu) | "
                    "structs=%llu (%.0f/s) | %.1f MB/s raw | %s bPDB (%.1fx) | ETA %s",
                    c->phase_label, pct,
                    (unsigned long long)(s->packs_done + s->packs_skipped),
                    (unsigned long long)s->packs_total,
                    (unsigned long long)s->packs_skipped,
                    (unsigned long long)s->packs_failed,
                    (unsigned long long)s->structures, rate, raw_mb_s, bpdb_str, ratio, eta_str);
        } else {
            fprintf(stderr,
                    "\r\033[K[%s] structs=%llu (%.0f/s) | %s raw (%.1f MB/s) | %s bPDB (%.2fx)",
                    c->phase_label,
                    (unsigned long long)s->structures, rate, raw_str, raw_mb_s, bpdb_str, ratio);
        }
        if (c->verify_every > 0) {
            fprintf(stderr, " | vfy=%llu/%llu",
                    (unsigned long long)s->verify_ok,
                    (unsigned long long)(s->verify_ok + s->verify_fail));
        }
        fflush(stderr);
        s->progress_on_line = true;
        s->last_print = t;
    }
    pthread_mutex_unlock(&s->mu);
}

static bool encode_cif_callback(const char* entry_name,
                                const uint8_t* raw_data, size_t raw_len,
                                void* user_data) {
    CifEncodeCtx* ctx = (CifEncodeCtx*)user_data;
    const char* cif_data = NULL;
    size_t      cif_len  = 0;

    if (!gz_inflate_or_passthrough(ctx->inf, raw_data, raw_len, &cif_data, &cif_len)) {
        log_line(ctx->stats, "ERROR: gzip inflate failed for '%s'", entry_name);
        return false;
    }

    const uint8_t* bpdb_ptr = NULL;
    size_t bpdb_len = bpdb_encode_reuse(ctx->enc, cif_data, cif_len, &bpdb_ptr);
    if (bpdb_len == 0 || !bpdb_ptr) {
        log_line(ctx->stats, "ERROR: bPDB encode failed for '%s'", entry_name);
        return false;
    }

    ctx->verify_seq++;
    ctx->local_count++;
    ctx->local_raw   += cif_len;
    ctx->local_bpdb  += bpdb_len;
    ctx->pack_structs++;
    ctx->pack_raw    += cif_len;
    ctx->pack_bpdb   += bpdb_len;

    if (ctx->verify_every > 0 && (ctx->verify_seq % ctx->verify_every) == 0) {
        char*  dec_buf = NULL;
        size_t dec_len = 0;
        FILE*  mf      = open_memstream(&dec_buf, &dec_len);
        if (mf) {
            decode_binary_to_afdb(bpdb_ptr, bpdb_len, mf);
            fclose(mf);
            if (dec_len == cif_len && memcmp(cif_data, dec_buf, cif_len) == 0) {
                ctx->local_vok++;
            } else {
                ctx->local_vfail++;
                log_line(ctx->stats, "WARNING: roundtrip verify mismatch on '%s' (orig=%zu dec=%zu)",
                         entry_name, cif_len, dec_len);
            }
            free(dec_buf);
        }
    }

    if (ctx->local_count >= PROGRESS_FLUSH_EVERY) {
        stats_flush_local(ctx->stats, ctx, false);
    }

    return pack_writer_add_record(ctx->pw, bpdb_ptr, bpdb_len);
}

static void build_gcs_url(char* out, size_t cap, const char* bucket, const char* fname) {
    snprintf(out, cap,
             "https://storage.googleapis.com/download/storage/v1/b/%s/o/proteomes%%2F%s?alt=media",
             bucket, fname);
}

/* ===================== Large Proteome Phase (LPT Scheduled) ===================== */

typedef struct {
    const ManifestShard* shards;
    const TaxGroup*      large_taxa;
    size_t               large_count;
    size_t               next_idx;
    pthread_mutex_t      idx_mu;

    const char*          bucket;
    const char*          out_dir;
    uint64_t             target_pack_bytes;
    uint64_t             verify_every;
    TokenManager*        tm;
    IngestStats*         stats;
} LargePhaseCtx;

static void* large_worker_main(void* arg) {
    LargePhaseCtx*  ctx = (LargePhaseCtx*)arg;
    GcsStreamer*    gs  = gcs_streamer_create();
    GzInflater*     inf = gz_inflater_create();
    BpdbEncoderCtx* enc = bpdb_encoder_create();

    if (!gs || !inf || !enc) {
        gcs_streamer_destroy(gs);
        gz_inflater_destroy(inf);
        bpdb_encoder_destroy(enc);
        return NULL;
    }

    while (true) {
        pthread_mutex_lock(&ctx->idx_mu);
        size_t idx = ctx->next_idx++;
        pthread_mutex_unlock(&ctx->idx_mu);
        if (idx >= ctx->large_count) break;

        const TaxGroup* tg = &ctx->large_taxa[idx];
        char base_path[480];
        snprintf(base_path, sizeof(base_path), "%s/proteomes/%02u/tax_%u",
                 ctx->out_dir, (unsigned)(tg->tax_id % 100U), (unsigned)tg->tax_id);

        uint64_t est_bpdb = (uint64_t)((double)tg->total_tar_bytes * 0.100);
        uint16_t planned  = pack_plan_parts(est_bpdb, ctx->target_pack_bytes);

        if (pack_already_complete(base_path, planned)) {
            pthread_mutex_lock(&ctx->stats->mu);
            ctx->stats->packs_skipped++;
            if (ctx->stats->est_total_bpdb > est_bpdb) ctx->stats->est_total_bpdb -= est_bpdb;
            pthread_mutex_unlock(&ctx->stats->mu);
            if (ctx->stats->verbose) {
                log_line(ctx->stats, "[large] SKIP  TaxID %u (already complete, %u part%s)",
                         (unsigned)tg->tax_id, (unsigned)planned, planned > 1 ? "s" : "");
            }
            continue;
        }

        double t_tax0 = now_sec();
        if (ctx->stats->verbose) {
            char sz_str[32];
            fmt_bytes(tg->total_tar_bytes, sz_str);
            log_line(ctx->stats, "[large] START TaxID %u (%zu shard%s, %s tar -> %u part%s)",
                     (unsigned)tg->tax_id, tg->shard_count, tg->shard_count > 1 ? "s" : "",
                     sz_str, (unsigned)planned, planned > 1 ? "s" : "");
        }

        bool tax_ok = false;
        uint64_t final_structs = 0, final_raw = 0, final_bpdb = 0, final_mfrags = 0;

        for (int retry = 0; retry < 3 && !tax_ok; retry++) {
            BpdbPackWriter pw;
            if (!pack_writer_init(&pw, base_path, tg->tax_id, 0, est_bpdb, ctx->target_pack_bytes)) break;

            CifEncodeCtx ectx = {
                .inf = inf, .enc = enc, .pw = &pw, .stats = ctx->stats,
                .phase_label = "large", .verify_every = ctx->verify_every
            };

            bool all_shards_ok = true;
            for (size_t s = 0; s < tg->shard_count; s++) {
                const ManifestShard* sh = &ctx->shards[tg->shard_start + s];
                char url[512];
                build_gcs_url(url, sizeof(url), ctx->bucket, sh->filename);

                if (ctx->stats->verbose && tg->shard_count > 1) {
                    char sh_sz[32];
                    fmt_bytes(sh->tar_bytes, sh_sz);
                    log_line(ctx->stats, "[large]   TaxID %u shard %zu/%zu: %s (%s)",
                             (unsigned)tg->tax_id, s + 1, tg->shard_count, sh->filename, sh_sz);
                }

                if (!gcs_stream_tar_url(gs, ctx->tm, url, true, sh->tar_bytes,
                                        encode_cif_callback, &ectx, NULL, NULL)) {
                    log_line(ctx->stats, "WARNING: TaxID %u shard %s failed (attempt %d/3)",
                             (unsigned)tg->tax_id, sh->filename, retry + 1);
                    all_shards_ok = false;
                    break;
                }
                stats_flush_local(ctx->stats, &ectx, false);
            }

            if (all_shards_ok) {
                final_structs = ectx.pack_structs;
                final_raw     = ectx.pack_raw;
                final_bpdb    = ectx.pack_bpdb;
                bool fin_ok   = pack_writer_finish(&pw);
                final_mfrags  = pw.multi_frag_proteins;
                if (fin_ok) tax_ok = true;
            }
            if (!tax_ok) {
                pack_writer_abort(&pw);
                if (retry < 2) sleep(2U << retry);
            }
        }

        pthread_mutex_lock(&ctx->stats->mu);
        if (tax_ok) {
            ctx->stats->packs_done++;
            ctx->stats->tars_done += tg->shard_count;
        } else {
            ctx->stats->packs_failed++;
        }
        pthread_mutex_unlock(&ctx->stats->mu);

        if (tax_ok && ctx->stats->verbose) {
            char b_str[32];
            fmt_bytes(final_bpdb, b_str);
            double dt = now_sec() - t_tax0;
            double ratio = (final_bpdb > 0) ? ((double)final_raw / (double)final_bpdb) : 0.0;
            log_line(ctx->stats,
                     "[large] DONE  TaxID %u: %llu structs, %s bPDB (%.2fx), %llu multi-frag united (%.1f s)",
                     (unsigned)tg->tax_id, (unsigned long long)final_structs,
                     b_str, ratio, (unsigned long long)final_mfrags, dt);
        } else if (!tax_ok) {
            log_line(ctx->stats, "ERROR: TaxID %u failed after 3 attempts!", (unsigned)tg->tax_id);
        }
    }

    bpdb_encoder_destroy(enc);
    gz_inflater_destroy(inf);
    gcs_streamer_destroy(gs);
    return NULL;
}

/* ===================== Small Taxa Clustered Phase ===================== */

typedef struct {
    const ManifestShard* shards;
    const size_t*        shard_indices;
    size_t               shard_count;
    size_t               next_pos;
    pthread_mutex_t      pos_mu;
    bool                 any_failed;

    const char*          bucket;
    uint64_t             verify_every;
    TokenManager*        tm;
    BpdbPackWriter*      pw;
    IngestStats*         stats;
} SmallBatchRunCtx;

static void* small_worker_main(void* arg) {
    SmallBatchRunCtx* ctx = (SmallBatchRunCtx*)arg;
    GcsStreamer*      gs  = gcs_streamer_create();
    GzInflater*       inf = gz_inflater_create();
    BpdbEncoderCtx*   enc = bpdb_encoder_create();

    if (!gs || !inf || !enc) {
        gcs_streamer_destroy(gs);
        gz_inflater_destroy(inf);
        bpdb_encoder_destroy(enc);
        return NULL;
    }

    CifEncodeCtx ectx = {
        .inf = inf, .enc = enc, .pw = ctx->pw, .stats = ctx->stats,
        .phase_label = "small", .verify_every = ctx->verify_every
    };

    while (true) {
        pthread_mutex_lock(&ctx->pos_mu);
        if (ctx->any_failed || ctx->next_pos >= ctx->shard_count) {
            pthread_mutex_unlock(&ctx->pos_mu);
            break;
        }
        size_t pos = ctx->next_pos++;
        pthread_mutex_unlock(&ctx->pos_mu);

        const ManifestShard* sh = &ctx->shards[ctx->shard_indices[pos]];
        char url[512];
        build_gcs_url(url, sizeof(url), ctx->bucket, sh->filename);

        if (!gcs_stream_tar_url(gs, ctx->tm, url, true, sh->tar_bytes,
                                encode_cif_callback, &ectx, NULL, NULL)) {
            log_line(ctx->stats, "ERROR: small shard %s failed!", sh->filename);
            pthread_mutex_lock(&ctx->pos_mu);
            ctx->any_failed = true;
            pthread_mutex_unlock(&ctx->pos_mu);
            break;
        }
        if ((pos & 63U) == 0) stats_flush_local(ctx->stats, &ectx, false);
    }

    stats_flush_local(ctx->stats, &ectx, false);
    bpdb_encoder_destroy(enc);
    gz_inflater_destroy(inf);
    gcs_streamer_destroy(gs);
    return NULL;
}

static bool run_small_batch(const SmallBatch* batch, size_t total_batches,
                            const TaxGroup* small_taxa, const ManifestShard* shards,
                            const char* bucket, const char* out_dir, int small_jobs,
                            uint64_t verify_every, TokenManager* tm, IngestStats* stats) {
    char base_path[480];
    snprintf(base_path, sizeof(base_path), "%s/small/small_%04u", out_dir, (unsigned)batch->batch_id);

    if (pack_already_complete(base_path, 1)) {
        uint64_t est = (uint64_t)((double)batch->total_tar_bytes * 0.105);
        pthread_mutex_lock(&stats->mu);
        stats->packs_skipped++;
        if (stats->est_total_bpdb > est) stats->est_total_bpdb -= est;
        pthread_mutex_unlock(&stats->mu);
        if (stats->verbose) {
            log_line(stats, "[small] SKIP  batch %04u/%zu (already complete)",
                     (unsigned)batch->batch_id, total_batches);
        }
        return true;
    }

    double t0 = now_sec();
    if (stats->verbose) {
        char sz_str[32];
        fmt_bytes(batch->total_tar_bytes, sz_str);
        log_line(stats, "[small] START batch %04u/%zu (%zu TaxIDs, %zu tars, %s tar)",
                 (unsigned)batch->batch_id, total_batches, batch->tax_count, batch->shard_count, sz_str);
    }

    size_t* sidx = malloc(batch->shard_count * sizeof(size_t));
    if (!sidx) return false;
    size_t k = 0;
    for (size_t i = 0; i < batch->tax_count; i++) {
        const TaxGroup* tg = &small_taxa[batch->tax_start + i];
        for (size_t s = 0; s < tg->shard_count; s++) sidx[k++] = tg->shard_start + s;
    }

    BpdbPackWriter pw;
    if (!pack_writer_init(&pw, base_path, 0, BPDB_PACK_FLAG_SMALL_CLUSTER, 0, 0)) {
        free(sidx);
        return false;
    }

    SmallBatchRunCtx bctx = {
        .shards = shards, .shard_indices = sidx, .shard_count = k, .next_pos = 0,
        .any_failed = false, .bucket = bucket, .verify_every = verify_every,
        .tm = tm, .pw = &pw, .stats = stats
    };
    pthread_mutex_init(&bctx.pos_mu, NULL);

    int nt = (small_jobs < 1) ? 1 : small_jobs;
    if ((size_t)nt > k && k > 0) nt = (int)k;

    pthread_t* threads = malloc((size_t)nt * sizeof(pthread_t));
    for (int i = 0; i < nt; i++) pthread_create(&threads[i], NULL, small_worker_main, &bctx);
    for (int i = 0; i < nt; i++) pthread_join(threads[i], NULL);
    free(threads);
    pthread_mutex_destroy(&bctx.pos_mu);
    free(sidx);

    uint64_t b_recs = 0, b_bytes = 0, b_mfrags = 0;
    bool ok = !bctx.any_failed && pack_writer_finish(&pw);
    if (ok) {
        b_recs   = pw.total_records;
        b_bytes  = pw.total_payload_bytes;
        b_mfrags = pw.multi_frag_proteins;
    } else {
        pack_writer_abort(&pw);
        pthread_mutex_lock(&stats->mu);
        stats->packs_failed++;
        pthread_mutex_unlock(&stats->mu);
        return false;
    }

    pthread_mutex_lock(&stats->mu);
    stats->packs_done++;
    stats->tars_done += k;
    pthread_mutex_unlock(&stats->mu);

    if (stats->verbose) {
        char b_str[32];
        fmt_bytes(b_bytes, b_str);
        log_line(stats, "[small] DONE  batch %04u/%zu: %llu structs, %s bPDB, %llu multi-frag united (%.1f s)",
                 (unsigned)batch->batch_id, total_batches, (unsigned long long)b_recs,
                 b_str, (unsigned long long)b_mfrags, now_sec() - t0);
    }
    return true;
}

/* ===================== Subcommand: manifest ===================== */

static int cmd_manifest(int argc, char** argv) {
    const char* phase        = "all";
    const char* bucket       = "public-datasets-deepmind-alphafold-v4";
    uint64_t    thresh_mib   = 100;
    uint64_t    target_mib   = 1024;
    int         large_jobs   = 6;
    int         small_jobs   = 48;
    uint64_t    verify_every = 0;
    size_t      limit        = 0;
    uint32_t    only_tax_id  = 0;
    bool        dry_run      = false;
    bool        verbose      = false;

    int i = 0;
    while (i < argc && argv[i][0] == '-') {
        if      (strcmp(argv[i], "--phase") == 0 && i + 1 < argc)         phase        = argv[++i];
        else if (strcmp(argv[i], "--bucket") == 0 && i + 1 < argc)        bucket       = argv[++i];
        else if (strcmp(argv[i], "--threshold-mib") == 0 && i + 1 < argc) thresh_mib   = strtoull(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--target-mib") == 0 && i + 1 < argc)    target_mib   = strtoull(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--large-jobs") == 0 && i + 1 < argc)    large_jobs   = atoi(argv[++i]);
        else if (strcmp(argv[i], "--small-jobs") == 0 && i + 1 < argc)    small_jobs   = atoi(argv[++i]);
        else if (strcmp(argv[i], "--verify-every") == 0 && i + 1 < argc)  verify_every = strtoull(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--limit") == 0 && i + 1 < argc)         limit        = (size_t)strtoull(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--tax-id") == 0 && i + 1 < argc)        only_tax_id  = (uint32_t)strtoul(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--dry-run") == 0)                       dry_run      = true;
        else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) verbose = true;
        else { fprintf(stderr, "Unknown option: %s\n", argv[i]); return 2; }
        i++;
    }

    if (i + 2 > argc) {
        fprintf(stderr, "Usage: bpdb_ingest manifest [options] <proteomes_manifest.txt> <out_dir>\n");
        return 2;
    }
    const char* manifest_path = argv[i];
    const char* out_dir       = argv[i + 1];

    ManifestShard* shards = NULL;
    size_t shard_count = 0;
    if (!parse_manifest_v4(manifest_path, &shards, &shard_count) || shard_count == 0) {
        fprintf(stderr, "Error: no *_v4.tar entries found in '%s'\n", manifest_path);
        free(shards);
        return 1;
    }

    size_t tax_cap = 65536, tax_count = 0;
    TaxGroup* taxa = malloc(tax_cap * sizeof(TaxGroup));
    uint64_t total_tar_bytes = 0;
    for (size_t s = 0; s < shard_count; ) {
        size_t e = s + 1;
        uint64_t sum = shards[s].tar_bytes;
        while (e < shard_count && shards[e].tax_id == shards[s].tax_id) {
            sum += shards[e].tar_bytes;
            e++;
        }
        if (only_tax_id == 0 || shards[s].tax_id == only_tax_id) {
            if (tax_count >= tax_cap) { tax_cap *= 2; taxa = realloc(taxa, tax_cap * sizeof(TaxGroup)); }
            taxa[tax_count++] = (TaxGroup){ shards[s].tax_id, sum, s, e - s };
            total_tar_bytes += sum;
        }
        s = e;
    }
    if (tax_count == 0) {
        fprintf(stderr, "Error: TaxID %u not found in '%s'\n", (unsigned)only_tax_id, manifest_path);
        free(taxa); free(shards);
        return 1;
    }

    uint64_t thresh_bytes = (only_tax_id != 0) ? 0ULL : (thresh_mib * MIB);
    uint64_t target_bytes = target_mib * MIB;

    TaxGroup* large_taxa = malloc(tax_count * sizeof(TaxGroup));
    TaxGroup* small_taxa = malloc(tax_count * sizeof(TaxGroup));
    size_t large_count = 0, small_count = 0;
    uint64_t large_tar = 0, small_tar = 0, large_pack_files = 0;

    for (size_t t = 0; t < tax_count; t++) {
        if (taxa[t].total_tar_bytes >= thresh_bytes) {
            large_taxa[large_count++] = taxa[t];
            large_tar += taxa[t].total_tar_bytes;
            uint64_t est = (uint64_t)((double)taxa[t].total_tar_bytes * 0.100);
            large_pack_files += pack_plan_parts(est, target_bytes);
        } else {
            small_taxa[small_count++] = taxa[t];
            small_tar += taxa[t].total_tar_bytes;
        }
    }

    if (large_count > 1) qsort(large_taxa, large_count, sizeof(TaxGroup), cmp_tax_desc);

    uint64_t small_est_bpdb = (uint64_t)((double)small_tar * 0.105);
    uint16_t n_small_batches = (small_count > 0) ? pack_plan_parts(small_est_bpdb, target_bytes) : 0;
    uint64_t small_batch_tar_target = (n_small_batches > 0) ? (small_tar / n_small_batches) : 0;

    SmallBatch* batches = calloc(n_small_batches + 1, sizeof(SmallBatch));
    size_t batch_count = 0;
    if (small_count > 0) {
        size_t cur_start = 0, cur_shards = 0;
        uint64_t cur_bytes = 0;
        for (size_t t = 0; t < small_count; t++) {
            if (cur_shards > 0 && batch_count + 1 < n_small_batches &&
                cur_bytes + small_taxa[t].total_tar_bytes > small_batch_tar_target) {
                batches[batch_count] = (SmallBatch){ (uint32_t)batch_count, cur_start, t - cur_start, cur_shards, cur_bytes };
                batch_count++;
                cur_start = t; cur_shards = 0; cur_bytes = 0;
            }
            cur_shards += small_taxa[t].shard_count;
            cur_bytes  += small_taxa[t].total_tar_bytes;
        }
        if (cur_shards > 0) {
            batches[batch_count] = (SmallBatch){ (uint32_t)batch_count, cur_start, small_count - cur_start, cur_shards, cur_bytes };
            batch_count++;
        }
    }

    char s_tot[32], s_lrg[32], s_sml[32], s_est[32];
    fmt_bytes(total_tar_bytes, s_tot);
    fmt_bytes(large_tar, s_lrg);
    fmt_bytes(small_tar, s_sml);
    fmt_bytes((uint64_t)((double)total_tar_bytes * 0.105), s_est);

    printf("========================================================================\n");
    printf("bPDB Ingest Plan (%s)\n", manifest_path);
    printf("========================================================================\n");
    printf("Total TaxIDs:     %zu (%s tar -> ~%s bPDB)\n", tax_count, s_tot, s_est);
    printf("Large proteomes:  %zu TaxIDs -> %llu balanced .bpdbpack files (%s tar)\n",
           large_count, (unsigned long long)large_pack_files, s_lrg);
    printf("Small taxa:       %zu TaxIDs -> %zu clustered .bpdbpack files (%s tar)\n",
           small_count, batch_count, s_sml);
    printf("Total pack files: %llu (target ~%llu MiB/pack, -F1 aligned, fragments united)\n",
           (unsigned long long)(large_pack_files + batch_count), (unsigned long long)target_mib);

    if (dry_run) {
        printf("\nTop 10 Largest Proteomes (balanced split plan):\n");
        size_t topn = (large_count < 10) ? large_count : 10;
        for (size_t t = 0; t < topn; t++) {
            uint64_t est = (uint64_t)((double)large_taxa[t].total_tar_bytes * 0.100);
            uint16_t k   = pack_plan_parts(est, target_bytes);
            char sb[32], se[32];
            fmt_bytes(large_taxa[t].total_tar_bytes, sb);
            fmt_bytes(est, se);
            printf("  TaxID %-8u | %3zu shards | %10s tar -> ~%10s bPDB (%2u balanced parts)\n",
                   (unsigned)large_taxa[t].tax_id, large_taxa[t].shard_count, sb, se, (unsigned)k);
        }
        free(batches); free(small_taxa); free(large_taxa); free(taxa); free(shards);
        return 0;
    }

    TokenManager tm;
    if (!token_manager_init(&tm)) {
        fprintf(stderr, "Error: failed to fetch GCS OAuth token via `gcloud auth print-access-token`.\n");
        free(batches); free(small_taxa); free(large_taxa); free(taxa); free(shards);
        return 1;
    }

    bool do_large = (strcmp(phase, "all") == 0 || strcmp(phase, "large") == 0);
    bool do_small = (strcmp(phase, "all") == 0 || strcmp(phase, "small") == 0) && (only_tax_id == 0);

    size_t run_large = (do_large && large_count > 0) ? ((limit > 0 && limit < large_count) ? limit : large_count) : 0;
    size_t run_small = (do_small && batch_count > 0) ? ((limit > 0 && limit < batch_count) ? limit : batch_count) : 0;

    uint64_t planned_tar_bytes = 0;
    for (size_t t = 0; t < run_large; t++) planned_tar_bytes += large_taxa[t].total_tar_bytes;
    for (size_t b = 0; b < run_small; b++) planned_tar_bytes += batches[b].total_tar_bytes;

    IngestStats stats;
    memset(&stats, 0, sizeof(stats));
    pthread_mutex_init(&stats.mu, NULL);
    stats.packs_total    = (uint64_t)(run_large + run_small);
    stats.est_total_bpdb = (uint64_t)((double)planned_tar_bytes * 0.105);
    stats.start_time     = now_sec();
    stats.last_print     = stats.start_time;
    stats.verbose        = verbose;

    if (run_large > 0) {
        LargePhaseCtx lctx = {
            .shards = shards, .large_taxa = large_taxa, .large_count = run_large, .next_idx = 0,
            .bucket = bucket, .out_dir = out_dir, .target_pack_bytes = target_bytes,
            .verify_every = verify_every, .tm = &tm, .stats = &stats
        };
        pthread_mutex_init(&lctx.idx_mu, NULL);

        int nt = (large_jobs < 1) ? 1 : large_jobs;
        if ((size_t)nt > run_large) nt = (int)run_large;
        pthread_t* threads = malloc((size_t)nt * sizeof(pthread_t));
        for (int t = 0; t < nt; t++) pthread_create(&threads[t], NULL, large_worker_main, &lctx);
        for (int t = 0; t < nt; t++) pthread_join(threads[t], NULL);
        free(threads);
        pthread_mutex_destroy(&lctx.idx_mu);
    }

    if (run_small > 0) {
        for (size_t b = 0; b < run_small; b++) {
            run_small_batch(&batches[b], run_small, small_taxa, shards, bucket, out_dir,
                            small_jobs, verify_every, &tm, &stats);
        }
    }

    if (stats.progress_on_line) fprintf(stderr, "\n");

    double elapsed = now_sec() - stats.start_time;
    printf("========================================================================\n");
    printf("Ingest Complete in %.1f s\n", elapsed);
    printf("Packs:       %llu written, %llu skipped (resumed), %llu failed\n",
           (unsigned long long)stats.packs_done, (unsigned long long)stats.packs_skipped,
           (unsigned long long)stats.packs_failed);
    printf("Structures:  %llu encoded (%.1f/s)\n",
           (unsigned long long)stats.structures,
           (elapsed > 0.0) ? (double)stats.structures / elapsed : 0.0);
    if (stats.bpdb_bytes > 0) {
        char rb[32], bb[32];
        fmt_bytes(stats.raw_bytes, rb);
        fmt_bytes(stats.bpdb_bytes, bb);
        printf("Volume:      %s raw mmCIF -> %s bPDB (%.2fx)\n",
               rb, bb, (double)stats.raw_bytes / (double)stats.bpdb_bytes);
    }
    if (verify_every > 0) {
        printf("Spot Verify: %llu passed, %llu failed\n",
               (unsigned long long)stats.verify_ok, (unsigned long long)stats.verify_fail);
    }
    printf("========================================================================\n");

    pthread_mutex_destroy(&stats.mu);
    token_manager_destroy(&tm);
    free(batches); free(small_taxa); free(large_taxa); free(taxa); free(shards);
    return (stats.packs_failed == 0 && stats.verify_fail == 0) ? 0 : 1;
}

/* ===================== Subcommand: pack (local tars or direct URLs) ===================== */

static bool is_http_url(const char* s) {
    return (strncmp(s, "http://", 7) == 0 || strncmp(s, "https://", 8) == 0);
}

static int cmd_pack(int argc, char** argv) {
    uint32_t tax_id       = 0;
    uint16_t pack_flags   = 0;
    uint64_t est_mib      = 0;
    uint64_t target_mib   = 1024;
    uint64_t verify_every = 0;
    bool     use_auth     = false;
    bool     verbose      = false;

    int i = 0;
    while (i < argc && argv[i][0] == '-') {
        if      (strcmp(argv[i], "--tax-id") == 0 && i + 1 < argc)       tax_id       = (uint32_t)strtoul(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--v6") == 0)                           pack_flags  |= BPDB_PACK_FLAG_AFDB_V6;
        else if (strcmp(argv[i], "--auth") == 0)                         use_auth     = true;
        else if (strcmp(argv[i], "--est-mib") == 0 && i + 1 < argc)      est_mib      = strtoull(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--target-mib") == 0 && i + 1 < argc)   target_mib   = strtoull(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--verify-every") == 0 && i + 1 < argc) verify_every = strtoull(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) verbose = true;
        else { fprintf(stderr, "Unknown option: %s\n", argv[i]); return 2; }
        i++;
    }

    if (i + 2 > argc) {
        fprintf(stderr, "Usage: bpdb_ingest pack [options] <out_base_path> <tar_or_url> [...]\n");
        return 2;
    }
    const char* out_base = argv[i++];

    TokenManager tm;
    bool have_tm = false;
    if (use_auth) {
        if (!token_manager_init(&tm)) {
            fprintf(stderr, "Error: failed to fetch GCS OAuth token.\n");
            return 1;
        }
        have_tm = true;
    }

    BpdbPackWriter pw;
    if (!pack_writer_init(&pw, out_base, tax_id, pack_flags, est_mib * MIB, target_mib * MIB)) {
        fprintf(stderr, "Error: failed to initialize pack writer for '%s'\n", out_base);
        if (have_tm) token_manager_destroy(&tm);
        return 1;
    }

    GcsStreamer*    gs  = gcs_streamer_create();
    GzInflater*     inf = gz_inflater_create();
    BpdbEncoderCtx* enc = bpdb_encoder_create();

    IngestStats stats;
    memset(&stats, 0, sizeof(stats));
    pthread_mutex_init(&stats.mu, NULL);
    stats.start_time = now_sec();
    stats.last_print = stats.start_time;
    stats.verbose    = verbose;

    CifEncodeCtx ectx = {
        .inf = inf, .enc = enc, .pw = &pw, .stats = &stats,
        .phase_label = "pack", .verify_every = verify_every
    };

    int total_inputs = argc - i;
    bool ok = true;
    for (int idx = 0; i < argc; i++, idx++) {
        const char* src = argv[i];
        if (verbose) {
            log_line(&stats, "[pack] Streaming input %d/%d: %s", idx + 1, total_inputs, src);
        }
        uint64_t cifs = 0, skipped = 0;
        double t0 = now_sec();
        bool s_ok = is_http_url(src)
            ? gcs_stream_tar_url(gs, have_tm ? &tm : NULL, src, use_auth, 0, encode_cif_callback, &ectx, &cifs, &skipped)
            : gcs_stream_tar_file(gs, src, encode_cif_callback, &ectx, &cifs, &skipped);
        stats_flush_local(&stats, &ectx, true);
        if (!s_ok) {
            log_line(&stats, "ERROR: failed processing '%s'", src);
            ok = false;
            break;
        }
        if (verbose) {
            log_line(&stats, "[pack] Finished input %d/%d: %llu .cif encoded, %llu non-cif skipped (%.1f s)",
                     idx + 1, total_inputs, (unsigned long long)cifs, (unsigned long long)skipped, now_sec() - t0);
        }
    }

    if (ok && verbose) {
        log_line(&stats, "[pack] Finalizing %zu held fragment candidates & sorting TOC...", pw.frag_count);
    }

    uint64_t multi_frags = 0;
    if (ok) {
        ok = pack_writer_finish(&pw);
        multi_frags = pw.multi_frag_proteins;
    } else {
        pack_writer_abort(&pw);
    }

    if (stats.progress_on_line) fprintf(stderr, "\n");
    double elapsed = now_sec() - stats.start_time;
    if (ok) {
        printf("Packed %llu structures (%llu multi-fragment proteins united) in %.2f s\n",
               (unsigned long long)stats.structures, (unsigned long long)multi_frags, elapsed);
        if (stats.bpdb_bytes > 0) {
            char rb[32], bb[32];
            fmt_bytes(stats.raw_bytes, rb);
            fmt_bytes(stats.bpdb_bytes, bb);
            printf("Volume: %s raw mmCIF -> %s bPDB (%.2fx)\n",
                   rb, bb, (double)stats.raw_bytes / (double)stats.bpdb_bytes);
        }
        if (verify_every > 0) {
            printf("Spot Verify: %llu passed, %llu failed\n",
                   (unsigned long long)stats.verify_ok, (unsigned long long)stats.verify_fail);
        }
    }

    pthread_mutex_destroy(&stats.mu);
    bpdb_encoder_destroy(enc);
    gz_inflater_destroy(inf);
    gcs_streamer_destroy(gs);
    if (have_tm) token_manager_destroy(&tm);
    return (ok && stats.verify_fail == 0) ? 0 : 1;
}

static void print_usage(const char* prog) {
    fprintf(stderr,
        "Usage:\n"
        "  %s manifest [options] <proteomes_manifest.txt> <out_dir>\n"
        "      --phase all|large|small   Which tier to ingest (default: all)\n"
        "      --threshold-mib <N>       Min TaxID tar MiB for standalone proteome (default: 100)\n"
        "      --target-mib <N>          Target bPDB MiB per pack/split part (default: 1024)\n"
        "      --large-jobs <N>          Parallel streaming workers for large proteomes (default: 6)\n"
        "      --small-jobs <N>          Parallel HTTP/2 workers for small taxa (default: 48)\n"
        "      --bucket <name>           GCS bucket (default: public-datasets-deepmind-alphafold-v4)\n"
        "      --tax-id <N>              Only ingest this single TaxID from manifest (great for tests)\n"
        "      --verify-every <N>        Spot-check 1 in N structures with lossless roundtrip\n"
        "      --limit <N>               Process at most N proteomes/batches (for testing)\n"
        "      -v, --verbose             Print per-proteome, per-shard, and finalize logs\n"
        "      --dry-run                 Print pack & split plan without downloading\n"
        "\n"
        "  %s pack [options] <out_base_path> <tar_or_url> [...]\n"
        "      --tax-id <N>              Tag pack with NCBI TaxID\n"
        "      --v6                      Set AFDB v6 overlay flag on pack\n"
        "      --auth                    Attach GCS OAuth Bearer token\n"
        "      --est-mib <N>             Estimated bPDB MiB (enables balanced multi-part split)\n"
        "      --target-mib <N>          Target bPDB MiB per part (default: 1024)\n"
        "      --verify-every <N>        Spot-check 1 in N structures\n"
        "      -v, --verbose             Print per-input and finalize logs\n",
        prog, prog);
}

int main(int argc, char** argv) {
    if (argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        print_usage(argv[0]);
        return (argc < 2) ? 2 : 0;
    }

    curl_global_init(CURL_GLOBAL_DEFAULT);
    int rc = 2;
    if      (strcmp(argv[1], "manifest") == 0) rc = cmd_manifest(argc - 2, argv + 2);
    else if (strcmp(argv[1], "pack") == 0)     rc = cmd_pack(argc - 2, argv + 2);
    else { fprintf(stderr, "Unknown subcommand: %s\n", argv[1]); print_usage(argv[0]); }
    curl_global_cleanup();
    return rc;
}
