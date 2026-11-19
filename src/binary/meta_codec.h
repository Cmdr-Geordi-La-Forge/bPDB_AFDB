#ifndef BPDB_META_CODEC_H
#define BPDB_META_CODEC_H

#include <stdbool.h>
#include <stddef.h>
#include "bpdb_internal.h"
#include "binary/byte_buf.h"

size_t meta_cold_encoded_size(const AfdbMetadata* meta, const ParsedChain* chain);
void   meta_write_cold(ByteWriter* w, const AfdbMetadata* meta, const ParsedChain* chain);
bool   meta_read_cold(ByteReader* r, ChainView* chain, AfdbMetadata* out_meta);

#endif // BPDB_META_CODEC_H
