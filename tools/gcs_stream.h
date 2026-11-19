#ifndef BPDB_GCS_STREAM_H
#define BPDB_GCS_STREAM_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <pthread.h>
#include <time.h>

typedef struct {
    pthread_mutex_t mu;
    char            token[2048];
    time_t          fetched_at;
    uint64_t        generation;
    bool            disabled;
} TokenManager;

bool     token_manager_init(TokenManager* tm);
void     token_manager_destroy(TokenManager* tm);
uint64_t token_manager_get(TokenManager* tm, char* dst, size_t dst_cap);
bool     token_manager_refresh(TokenManager* tm, uint64_t failed_gen);

typedef struct GcsStreamer GcsStreamer;

GcsStreamer* gcs_streamer_create(void);
void         gcs_streamer_destroy(GcsStreamer* gs);

typedef bool (*TarCifCallback)(const char* entry_name,
                               const uint8_t* raw_data, size_t raw_len,
                               void* user_data);

bool gcs_stream_tar_url(GcsStreamer* gs, TokenManager* tm,
                        const char* url, bool use_auth, uint64_t known_size,
                        TarCifCallback cb, void* user_data,
                        uint64_t* out_cif_count, uint64_t* out_skipped);

bool gcs_stream_tar_file(GcsStreamer* gs, const char* path,
                         TarCifCallback cb, void* user_data,
                         uint64_t* out_cif_count, uint64_t* out_skipped);

#endif // BPDB_GCS_STREAM_H
