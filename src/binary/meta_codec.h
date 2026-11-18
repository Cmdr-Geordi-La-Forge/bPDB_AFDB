#ifndef BPDB_META_CODEC_H
#define BPDB_META_CODEC_H

#include <stdbool.h>
#include <stddef.h>
#include "bpdb_internal.h"
#include "binary/byte_buf.h"

size_t meta_encoded_size(const AfdbMetadata* meta);
void   meta_write_binary(ByteWriter* w, const AfdbMetadata* meta);
bool   meta_read_binary(ByteReader* r, AfdbMetadata* out_meta);

#endif // BPDB_META_CODEC_H
