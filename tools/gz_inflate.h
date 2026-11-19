#ifndef BPDB_GZ_INFLATE_H
#define BPDB_GZ_INFLATE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct GzInflater GzInflater;

GzInflater* gz_inflater_create(void);
void        gz_inflater_destroy(GzInflater* inf);

// If src is gzipped (0x1f 0x8b), inflates into inf's reusable buffer and sets
// *out_data / *out_len (valid until the next call on the same GzInflater).
// If src is plain text, points *out_data directly to src (zero-copy).
bool gz_inflate_or_passthrough(GzInflater* inf, const uint8_t* src, size_t src_len,
                               const char** out_data, size_t* out_len);

#endif // BPDB_GZ_INFLATE_H
