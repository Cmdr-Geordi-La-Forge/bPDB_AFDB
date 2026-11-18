#include "bpdb.h"
#include "bpdb_internal.h"
#include "binary/byte_buf.h"
#include "binary/meta_codec.h"
#include "binary/coord_codec.h"
#include "cif/cif_parser.h"
#include <stdlib.h>

size_t encode_afdb_to_binary(const char* cif_data, size_t cif_len, uint8_t** out_buf) {
    AfdbMetadata meta  = {0};
    ParsedChain  chain = {0};

    if (!cif_parse_afdb(cif_data, cif_len, &meta, &chain)) {
        afdb_meta_free(&meta);
        parsed_chain_free(&chain);
        *out_buf = NULL;
        return 0;
    }

    size_t max_out_size = meta_encoded_size(&meta) + coord_max_chain_size(&chain);
    uint8_t* buf = calloc(1, max_out_size);
    if (!buf) {
        afdb_meta_free(&meta);
        parsed_chain_free(&chain);
        *out_buf = NULL;
        return 0;
    }

    ByteWriter w = bw_init(buf, max_out_size);
    meta_write_binary(&w, &meta);
    coord_write_chain(&w, &chain);

    afdb_meta_free(&meta);
    parsed_chain_free(&chain);

    *out_buf = buf;
    return bw_len(&w);
}
