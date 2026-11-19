#include "gz_inflate.h"
#include <stdlib.h>
#include <string.h>

#ifdef BPDB_USE_LIBDEFLATE
    #include <libdeflate.h>
#else
    #include <zlib.h>
#endif

struct GzInflater {
#ifdef BPDB_USE_LIBDEFLATE
    struct libdeflate_decompressor* decomp;
#else
    z_stream zs;
    bool     zs_ready;
#endif
    uint8_t* buf;
    size_t   cap;
};

static bool ensure_cap(GzInflater* inf, size_t need) {
    if (need <= inf->cap) return true;
    size_t next = (inf->cap == 0) ? (256 * 1024) : inf->cap;
    while (next < need) next *= 2;
    uint8_t* grown = realloc(inf->buf, next);
    if (!grown) return false;
    inf->buf = grown;
    inf->cap = next;
    return true;
}

GzInflater* gz_inflater_create(void) {
    GzInflater* inf = calloc(1, sizeof(*inf));
    if (!inf) return NULL;

#ifdef BPDB_USE_LIBDEFLATE
    inf->decomp = libdeflate_alloc_decompressor();
    if (!inf->decomp) { free(inf); return NULL; }
#else
    if (inflateInit2(&inf->zs, 16 + MAX_WBITS) != Z_OK) {
        free(inf);
        return NULL;
    }
    inf->zs_ready = true;
#endif
    ensure_cap(inf, 512 * 1024);
    return inf;
}

void gz_inflater_destroy(GzInflater* inf) {
    if (!inf) return;
#ifdef BPDB_USE_LIBDEFLATE
    if (inf->decomp) libdeflate_free_decompressor(inf->decomp);
#else
    if (inf->zs_ready) inflateEnd(&inf->zs);
#endif
    free(inf->buf);
    free(inf);
}

bool gz_inflate_or_passthrough(GzInflater* inf, const uint8_t* src, size_t src_len,
                               const char** out_data, size_t* out_len) {
    // Plain uncompressed mmCIF? Zero-copy passthrough.
    if (src_len < 2 || src[0] != 0x1F || src[1] != 0x8B) {
        *out_data = (const char*)src;
        *out_len  = src_len;
        return true;
    }

#ifdef BPDB_USE_LIBDEFLATE
    if (src_len < 18) return false;
    uint32_t isize = 0;
    memcpy(&isize, src + src_len - 4, 4);
    size_t hint = (isize > 0 && isize <= (64U * 1024U * 1024U)) ? (size_t)isize : (src_len * 6);
    if (!ensure_cap(inf, hint + 1)) return false;

    while (true) {
        size_t actual = 0;
        enum libdeflate_result res = libdeflate_gzip_decompress(
            inf->decomp, src, src_len, inf->buf, inf->cap - 1, &actual
        );
        if (res == LIBDEFLATE_SUCCESS) {
            inf->buf[actual] = '\0';
            *out_data = (const char*)inf->buf;
            *out_len  = actual;
            return true;
        }
        if (res == LIBDEFLATE_INSUFFICIENT_SPACE) {
            if (inf->cap >= (128U * 1024U * 1024U)) return false;
            if (!ensure_cap(inf, inf->cap * 2)) return false;
            continue;
        }
        return false;
    }
#else
    if (inflateReset(&inf->zs) != Z_OK) return false;
    inf->zs.next_in  = (Bytef*)src;
    inf->zs.avail_in = (uInt)src_len;

    size_t total = 0;
    while (true) {
        if (!ensure_cap(inf, total + (256 * 1024) + 1)) return false;
        inf->zs.next_out  = inf->buf + total;
        inf->zs.avail_out = (uInt)(inf->cap - 1 - total);

        int ret = inflate(&inf->zs, Z_NO_FLUSH);
        total = (inf->cap - 1) - inf->zs.avail_out;

        if (ret == Z_STREAM_END) {
            // Support multi-member concatenated gzip streams
            if (inf->zs.avail_in >= 2 &&
                inf->zs.next_in[0] == 0x1F && inf->zs.next_in[1] == 0x8B) {
                if (inflateReset(&inf->zs) != Z_OK) return false;
                continue;
            }
            inf->buf[total] = '\0';
            *out_data = (const char*)inf->buf;
            *out_len  = total;
            return true;
        }
        if (ret != Z_OK && ret != Z_BUF_ERROR) {
            return false;
        }
        if (ret == Z_BUF_ERROR && inf->zs.avail_out > 0) {
            return false; // truncated input
        }
    }
#endif
}
