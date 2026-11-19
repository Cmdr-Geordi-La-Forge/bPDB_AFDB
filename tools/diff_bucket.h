#ifndef BPDB_DIFF_BUCKET_H
#define BPDB_DIFF_BUCKET_H

#include <stddef.h>
#include <stdint.h>
#include <pthread.h>

typedef struct {
    char     category[64];
    uint64_t count;
    size_t   saved_samples;
    char     first_example[128];
} DiffBucket;

typedef struct {
    pthread_mutex_t mu;
    DiffBucket*     items;
    size_t          count;
    size_t          cap;
    const char*     dump_dir;
    size_t          max_samples_per_bucket;
} DiffTracker;

void diff_tracker_init(DiffTracker* dt, const char* dump_dir, size_t max_samples);
void diff_tracker_destroy(DiffTracker* dt);

void diff_tracker_record_error(DiffTracker* dt, const char* entry_name,
                               const char* category,
                               const char* orig, size_t orig_len);

void diff_tracker_record_mismatch(DiffTracker* dt, const char* entry_name,
                                  const char* orig, size_t orig_len,
                                  const char* dec, size_t dec_len);

void diff_tracker_print_report(DiffTracker* dt);

#endif // BPDB_DIFF_BUCKET_H
