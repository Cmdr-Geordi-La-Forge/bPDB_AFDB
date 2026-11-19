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

    HotPrep hot = {0};
    if (!coord_prepare_v2(&chain, &meta, &hot)) {
        coord_hot_prep_free(&hot);
        afdb_meta_free(&meta);
        parsed_chain_free(&chain);
        *out_buf = NULL;
        return 0;
    }

    size_t hot_size  = coord_hot_encoded_size(&hot);
    size_t cold_size = meta_cold_encoded_size(&meta);
    size_t total     = hot_size + cold_size;

    uint8_t* buf = calloc(1, total);
    if (!buf) {
        coord_hot_prep_free(&hot);
        afdb_meta_free(&meta);
        parsed_chain_free(&chain);
        *out_buf = NULL;
        return 0;
    }

    ByteWriter w = bw_init(buf, total);
    coord_write_v2_hot(&w, &meta, &chain, &hot, (uint32_t)hot_size);
    meta_write_cold(&w, &meta, &chain);

    coord_hot_prep_free(&hot);
    afdb_meta_free(&meta);
    parsed_chain_free(&chain);

    *out_buf = buf;
    return bw_len(&w);
}
