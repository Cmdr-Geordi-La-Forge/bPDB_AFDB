#define _POSIX_C_SOURCE 200809L

#include "bpdb.h"
#include "tar_stream.h"
#include "gz_inflate.h"
#include "diff_bucket.h"

#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef struct {
    char*    name;
    uint8_t* data;
    size_t   size;
} WorkItem;

typedef struct {
    WorkItem*       ring;
    size_t          cap;
    size_t          head;
    size_t          tail;
    size_t          count;
    bool            closed;
    pthread_mutex_t mu;
    pthread_cond_t  not_empty;
    pthread_cond_t  not_full;
} WorkQueue;

typedef struct {
    pthread_mutex_t mu;
    uint64_t        processed;
    uint64_t        passed;
    uint64_t        failed;
    uint64_t        gz_bytes;
    uint64_t        raw_bytes;
    uint64_t        bpdb_bytes;
    double          start_time;
    double          last_report_time;
} RunStats;

typedef struct {
    WorkQueue*   queue;
    DiffTracker* tracker;
    RunStats*    stats;
    bool         bench_only;
    bool         quiet;
} WorkerCtx;

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void wq_init(WorkQueue* q, size_t cap) {
    memset(q, 0, sizeof(*q));
    q->ring = calloc(cap, sizeof(WorkItem));
    q->cap  = cap;
    pthread_mutex_init(&q->mu, NULL);
    pthread_cond_init(&q->not_empty, NULL);
    pthread_cond_init(&q->not_full, NULL);
}

static void wq_destroy(WorkQueue* q) {
    pthread_cond_destroy(&q->not_full);
    pthread_cond_destroy(&q->not_empty);
    pthread_mutex_destroy(&q->mu);
    free(q->ring);
}

static void wq_push(WorkQueue* q, WorkItem item) {
    pthread_mutex_lock(&q->mu);
    while (q->count == q->cap) {
        pthread_cond_wait(&q->not_full, &q->mu);
    }
    q->ring[q->tail] = item;
    q->tail = (q->tail + 1) % q->cap;
    q->count++;
    pthread_cond_signal(&q->not_empty);
    pthread_mutex_unlock(&q->mu);
}

static bool wq_pop(WorkQueue* q, WorkItem* out) {
    pthread_mutex_lock(&q->mu);
    while (q->count == 0 && !q->closed) {
        pthread_cond_wait(&q->not_empty, &q->mu);
    }
    if (q->count == 0 && q->closed) {
        pthread_mutex_unlock(&q->mu);
        return false;
    }
    *out = q->ring[q->head];
    q->head = (q->head + 1) % q->cap;
    q->count--;
    pthread_cond_signal(&q->not_full);
    pthread_mutex_unlock(&q->mu);
    return true;
}

static void wq_close(WorkQueue* q) {
    pthread_mutex_lock(&q->mu);
    q->closed = true;
    pthread_cond_broadcast(&q->not_empty);
    pthread_mutex_unlock(&q->mu);
}

static void record_progress(WorkerCtx* ctx, bool ok, size_t gz_sz, size_t raw_sz, size_t bpdb_sz) {
    RunStats* s = ctx->stats;
    pthread_mutex_lock(&s->mu);
    s->processed++;
    if (ok) s->passed++;
    else    s->failed++;
    s->gz_bytes   += gz_sz;
    s->raw_bytes  += raw_sz;
    s->bpdb_bytes += bpdb_sz;

    if (!ctx->quiet) {
        double t = now_sec();
        if (t - s->last_report_time >= 1.0) {
            double elapsed = t - s->start_time;
            double mb_s    = (elapsed > 0.0) ? ((double)s->raw_bytes / (1024.0 * 1024.0)) / elapsed : 0.0;
            double rate    = (elapsed > 0.0) ? (double)s->processed / elapsed : 0.0;
            double ratio   = (s->bpdb_bytes > 0) ? (double)s->raw_bytes / (double)s->bpdb_bytes : 0.0;
            fprintf(stderr, "\r[progress] %llu files | ok=%llu fail=%llu | %.1f files/s | %.1f MB/s raw | %.2fx",
                    (unsigned long long)s->processed,
                    (unsigned long long)s->passed,
                    (unsigned long long)s->failed,
                    rate, mb_s, ratio);
            s->last_report_time = t;
        }
    }
    pthread_mutex_unlock(&s->mu);
}

static void* worker_main(void* arg) {
    WorkerCtx*  ctx = (WorkerCtx*)arg;
    GzInflater* inf = gz_inflater_create();
    WorkItem    item;

    while (wq_pop(ctx->queue, &item)) {
        const char* cif_data = NULL;
        size_t      cif_len  = 0;

        if (!inf || !gz_inflate_or_passthrough(inf, item.data, item.size, &cif_data, &cif_len)) {
            diff_tracker_record_error(ctx->tracker, item.name, "gzip_inflate_error", NULL, 0);
            record_progress(ctx, false, item.size, 0, 0);
            free(item.data);
            free(item.name);
            continue;
        }

        uint8_t* bpdb_buf = NULL;
        size_t   bpdb_len = encode_afdb_to_binary(cif_data, cif_len, &bpdb_buf);

        bool ok = false;
        if (bpdb_len == 0 || !bpdb_buf) {
            diff_tracker_record_error(ctx->tracker, item.name, "encode_returned_0", cif_data, cif_len);
        } else if (ctx->bench_only) {
            ok = true;
        } else {
            char*  dec_buf = NULL;
            size_t dec_len = 0;
            FILE*  mem_f   = open_memstream(&dec_buf, &dec_len);
            if (!mem_f) {
                diff_tracker_record_error(ctx->tracker, item.name, "open_memstream_failed", cif_data, cif_len);
            } else {
                decode_binary_to_afdb(bpdb_buf, bpdb_len, mem_f);
                fclose(mem_f);

                if (dec_len == cif_len && memcmp(cif_data, dec_buf, cif_len) == 0) {
                    ok = true;
                } else {
                    diff_tracker_record_mismatch(ctx->tracker, item.name,
                                                 cif_data, cif_len, dec_buf, dec_len);
                }
                free(dec_buf);
            }
        }

        record_progress(ctx, ok, item.size, cif_len, bpdb_len);
        free(bpdb_buf);
        free(item.data);
        free(item.name);
    }

    gz_inflater_destroy(inf);
    return NULL;
}

static bool has_suffix(const char* s, const char* suf) {
    size_t sl = strlen(s), sufl = strlen(suf);
    return (sl >= sufl) && (memcmp(s + sl - sufl, suf, sufl) == 0);
}

static bool enqueue_single_file(WorkQueue* q, const char* path) {
    FILE* f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "Error: cannot open '%s'\n", path);
        return false;
    }
    fseek(f, 0, SEEK_END);
    long fsz = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (fsz <= 0) { fclose(f); return false; }

    uint8_t* data = malloc((size_t)fsz);
    if (!data || fread(data, 1, (size_t)fsz, f) != (size_t)fsz) {
        free(data);
        fclose(f);
        return false;
    }
    fclose(f);

    WorkItem item = { .name = strdup(path), .data = data, .size = (size_t)fsz };
    wq_push(q, item);
    return true;
}

static void print_usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s [options] <input.tar | - | file.cif[.gz]> [...]\n"
        "\n"
        "Options:\n"
        "  -j <N>    Number of worker threads (default: online CPUs)\n"
        "  -o <dir>  Directory to dump failing sample .cif & .diff.txt files\n"
        "  -k <N>    Max failing samples dumped per diff bucket (default: 5)\n"
        "  -n <N>    Stop after enqueuing N structures (0 = unlimited)\n"
        "  -b        Benchmark mode (encode-only; skips decode & round-trip diff)\n"
        "  -q        Quiet mode (suppress periodic progress lines)\n"
        "  -h        Show this help\n",
        prog);
}

int main(int argc, char** argv) {
    long default_cpus = sysconf(_SC_NPROCESSORS_ONLN);
    int         num_threads = (default_cpus > 0) ? (int)default_cpus : 4;
    const char* dump_dir    = NULL;
    size_t      max_samples = 5;
    uint64_t    max_entries = 0;
    bool        bench_only  = false;
    bool        quiet       = false;

    int first_input = 1;
    while (first_input < argc && argv[first_input][0] == '-' && argv[first_input][1] != '\0') {
        const char* opt = argv[first_input];
        if (strcmp(opt, "-h") == 0 || strcmp(opt, "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        } else if (strcmp(opt, "-b") == 0) {
            bench_only = true;
        } else if (strcmp(opt, "-q") == 0) {
            quiet = true;
        } else if (strcmp(opt, "-j") == 0 && first_input + 1 < argc) {
            num_threads = atoi(argv[++first_input]);
            if (num_threads < 1) num_threads = 1;
        } else if (strcmp(opt, "-o") == 0 && first_input + 1 < argc) {
            dump_dir = argv[++first_input];
        } else if (strcmp(opt, "-k") == 0 && first_input + 1 < argc) {
            max_samples = (size_t)strtoull(argv[++first_input], NULL, 10);
        } else if (strcmp(opt, "-n") == 0 && first_input + 1 < argc) {
            max_entries = strtoull(argv[++first_input], NULL, 10);
        } else {
            fprintf(stderr, "Unknown option: %s\n", opt);
            print_usage(argv[0]);
            return 2;
        }
        first_input++;
    }

    if (first_input >= argc) {
        print_usage(argv[0]);
        return 2;
    }

    size_t q_cap = (size_t)num_threads * 4;
    if (q_cap < 16)  q_cap = 16;
    if (q_cap > 256) q_cap = 256;

    WorkQueue   queue;
    DiffTracker tracker;
    RunStats    stats;

    wq_init(&queue, q_cap);
    diff_tracker_init(&tracker, dump_dir, max_samples);
    memset(&stats, 0, sizeof(stats));
    pthread_mutex_init(&stats.mu, NULL);
    stats.start_time       = now_sec();
    stats.last_report_time = stats.start_time;

    WorkerCtx ctx = {
        .queue      = &queue,
        .tracker    = &tracker,
        .stats      = &stats,
        .bench_only = bench_only,
        .quiet      = quiet
    };

    pthread_t* threads = malloc((size_t)num_threads * sizeof(pthread_t));
    for (int i = 0; i < num_threads; i++) {
        pthread_create(&threads[i], NULL, worker_main, &ctx);
    }

    uint64_t enqueued = 0;
    uint64_t skipped  = 0;

    for (int a = first_input; a < argc; a++) {
        if (max_entries > 0 && enqueued >= max_entries) break;
        const char* path = argv[a];

        if (strcmp(path, "-") == 0 || has_suffix(path, ".tar")) {
            TarStream ts;
            if (!tar_stream_open(&ts, path)) {
                fprintf(stderr, "Error: cannot open tar '%s'\n", path);
                continue;
            }
            TarEntry entry;
            bool tar_err = false;
            while (tar_stream_next_cif(&ts, &entry, &tar_err)) {
                WorkItem item = {
                    .name = strdup(entry.name),
                    .data = entry.data,
                    .size = entry.size
                };
                wq_push(&queue, item);
                enqueued++;
                if (max_entries > 0 && enqueued >= max_entries) break;
            }
            if (tar_err) {
                fprintf(stderr, "\nWarning: tar stream '%s' ended with read/format error\n", path);
            }
            skipped += ts.skipped_entries;
            tar_stream_close(&ts);
        } else {
            if (enqueue_single_file(&queue, path)) {
                enqueued++;
            }
        }
    }

    wq_close(&queue);
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }
    free(threads);

    double elapsed = now_sec() - stats.start_time;
    if (!quiet) fprintf(stderr, "\n");

    double raw_mb  = (double)stats.raw_bytes  / (1024.0 * 1024.0);
    double gz_mb   = (double)stats.gz_bytes   / (1024.0 * 1024.0);
    double bpdb_mb = (double)stats.bpdb_bytes / (1024.0 * 1024.0);

    printf("============================================================\n");
    printf("Mode:        %s (%d worker threads)\n", bench_only ? "encode-only benchmark" : "full round-trip verify", num_threads);
    printf("Elapsed:     %.2f s\n", elapsed);
    printf("Structures:  %llu processed (%llu passed, %llu failed, %llu non-cif skipped)\n",
           (unsigned long long)stats.processed,
           (unsigned long long)stats.passed,
           (unsigned long long)stats.failed,
           (unsigned long long)skipped);
    printf("Volume:      %.2f MB input -> %.2f MB raw mmCIF -> %.2f MB bPDB\n", gz_mb, raw_mb, bpdb_mb);
    if (stats.bpdb_bytes > 0) {
        printf("Compression: %.2fx vs raw mmCIF | %.2fx vs .cif.gz\n",
               (double)stats.raw_bytes / (double)stats.bpdb_bytes,
               (double)stats.gz_bytes  / (double)stats.bpdb_bytes);
    }
    if (elapsed > 0.0) {
        printf("Throughput:  %.1f structures/s | %.1f MB/s raw mmCIF\n",
               (double)stats.processed / elapsed,
               raw_mb / elapsed);
    }
    printf("============================================================\n");

    diff_tracker_print_report(&tracker);

    pthread_mutex_destroy(&stats.mu);
    diff_tracker_destroy(&tracker);
    wq_destroy(&queue);

    return (stats.failed == 0) ? 0 : 1;
}
