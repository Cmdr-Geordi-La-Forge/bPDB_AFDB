#ifndef BPDB_TAR_STREAM_H
#define BPDB_TAR_STREAM_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

typedef struct {
    FILE*    fp;
    bool     close_on_destroy;
    bool     is_seekable;
    char     pending_name[4096];
    bool     has_pending_name;
    uint64_t skipped_entries;
} TarStream;

typedef struct {
    char     name[4096];
    uint8_t* data;
    size_t   size;
} TarEntry;

bool tar_stream_open(TarStream* ts, const char* path);
void tar_stream_close(TarStream* ts);

// Returns true if a matching .cif / .cif.gz entry was loaded into out_entry.
// Caller owns out_entry->data and must free() it.
// Returns false on clean EOF; sets *out_error = true on truncated/corrupt archive.
bool tar_stream_next_cif(TarStream* ts, TarEntry* out_entry, bool* out_error);

#endif // BPDB_TAR_STREAM_H
