#ifndef BPDB_COORD_CODEC_H
#define BPDB_COORD_CODEC_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "bpdb_internal.h"
#include "binary/byte_buf.h"
#include "simd/simd.h"

typedef struct {
    simd_vec_t     acc;
    const uint8_t* ptr;
    const uint8_t* end;
} CoordStream;

size_t coord_max_chain_size(const ParsedChain* chain);
void   coord_write_chain(ByteWriter* w, const ParsedChain* chain);
bool   coord_read_chain_view(ByteReader* r, ChainView* out_view);

void   coord_stream_init(CoordStream* cs, const ChainView* view);
bool   coord_stream_step(CoordStream* cs, int32_t* b, int32_t* x, int32_t* y, int32_t* z);

#endif // BPDB_COORD_CODEC_H
