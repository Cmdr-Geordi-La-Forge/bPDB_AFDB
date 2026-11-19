#include "bpdb.h"
#include "bpdb_internal.h"
#include "binary/byte_buf.h"
#include "binary/meta_codec.h"
#include "binary/coord_codec.h"
#include "cif/cif_parser.h"
#include <stdlib.h>

struct BpdbEncoderCtx {
    AfdbMetadata meta;
    ParsedChain  chain;
    HotPrep      hot;
    uint8_t*     out_buf;
    size_t       out_cap;
};

BpdbEncoderCtx* bpdb_encoder_create(void) {
    return calloc(1, sizeof(BpdbEncoderCtx));
}

void bpdb_encoder_destroy(BpdbEncoderCtx* ctx) {
    if (!ctx) return;
    coord_hot_prep_free(&ctx->hot);
    afdb_meta_free(&ctx->meta);
    parsed_chain_free(&ctx->chain);
    free(ctx->out_buf);
    free(ctx);
}

size_t bpdb_encode_reuse(BpdbEncoderCtx* ctx, const char* cif_data, size_t cif_len,
                         const uint8_t** out_ptr) {
    if (out_ptr) *out_ptr = NULL;
    if (!ctx || !cif_data || cif_len == 0 || !out_ptr) return 0;

    afdb_meta_reset_reuse(&ctx->meta);

    if (!cif_parse_afdb(cif_data, cif_len, &ctx->meta, &ctx->chain)) {
        afdb_meta_reset_reuse(&ctx->meta);
        return 0;
    }

    if (!coord_prepare_v2(&ctx->chain, &ctx->meta, &ctx->hot)) {
        afdb_meta_reset_reuse(&ctx->meta);
        return 0;
    }

    size_t hot_size  = coord_hot_encoded_size(&ctx->hot);
    size_t cold_size = meta_cold_encoded_size(&ctx->meta, &ctx->chain);
    size_t total     = hot_size + cold_size;

    if (ctx->out_cap < total) {
        size_t ncap = (ctx->out_cap == 0) ? 32768U : ctx->out_cap;
        while (ncap < total) ncap *= 2U;
        uint8_t* nb = realloc(ctx->out_buf, ncap);
        if (!nb) {
            afdb_meta_reset_reuse(&ctx->meta);
            return 0;
        }
        ctx->out_buf = nb;
        ctx->out_cap = ncap;
    }

    ByteWriter w = bw_init(ctx->out_buf, total);
    coord_write_v2_hot(&w, &ctx->meta, &ctx->chain, &ctx->hot, (uint32_t)hot_size);
    meta_write_cold(&w, &ctx->meta, &ctx->chain);

    afdb_meta_reset_reuse(&ctx->meta);

    *out_ptr = ctx->out_buf;
    return bw_len(&w);
}

size_t encode_afdb_to_binary(const char* cif_data, size_t cif_len, uint8_t** out_buf) {
    if (!out_buf) return 0;
    *out_buf = NULL;
    if (!cif_data || cif_len == 0) return 0;

    AfdbMetadata meta  = {0};
    ParsedChain  chain = {0};

    if (!cif_parse_afdb(cif_data, cif_len, &meta, &chain)) {
        afdb_meta_free(&meta);
        parsed_chain_free(&chain);
        return 0;
    }

    HotPrep hot = {0};
    if (!coord_prepare_v2(&chain, &meta, &hot)) {
        coord_hot_prep_free(&hot);
        afdb_meta_free(&meta);
        parsed_chain_free(&chain);
        return 0;
    }

    size_t hot_size  = coord_hot_encoded_size(&hot);
    size_t cold_size = meta_cold_encoded_size(&meta, &chain);
    size_t total     = hot_size + cold_size;

    uint8_t* buf = malloc(total);
    if (!buf) {
        coord_hot_prep_free(&hot);
        afdb_meta_free(&meta);
        parsed_chain_free(&chain);
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
