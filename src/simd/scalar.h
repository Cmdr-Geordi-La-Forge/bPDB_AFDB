#ifndef BPDB_SIMD_SCALAR_H
#define BPDB_SIMD_SCALAR_H

#include <stdint.h>
#include <string.h>

typedef struct {
    int32_t b, x, y, z;
} simd_vec_t;

static inline simd_vec_t simd_init(void) {
    simd_vec_t v = {0, 0, 0, 0};
    return v;
}

static inline int32_t scalar_unpack_anchor(int16_t msb_raw, int16_t lsb_raw) {
    uint32_t hi = (uint32_t)(int32_t)(msb_raw >> 1) << 15;
    uint32_t lo = (uint32_t)(int32_t)(lsb_raw >> 1) & 0x7FFFU;
    return (int32_t)(hi | lo);
}

static inline void simd_decode_step(simd_vec_t* acc, const uint8_t** stream_ptr,
                                    int32_t* b, int32_t* x, int32_t* y, int32_t* z) {
    const uint8_t* p = *stream_ptr;

    if ((p[0] & 1U) == 1U) {
        int16_t s[8];
        memcpy(s, p, sizeof(s));
        acc->b = scalar_unpack_anchor(s[0], s[4]);
        acc->x = scalar_unpack_anchor(s[1], s[5]);
        acc->y = scalar_unpack_anchor(s[2], s[6]);
        acc->z = scalar_unpack_anchor(s[3], s[7]);
        *stream_ptr = p + 16;
    } else {
        int16_t s[4];
        memcpy(s, p, sizeof(s));
        acc->b += (s[0] >> 1);
        acc->x += (s[1] >> 1);
        acc->y += (s[2] >> 1);
        acc->z += (s[3] >> 1);
        *stream_ptr = p + 8;
    }

    *b = acc->b;
    *x = acc->x;
    *y = acc->y;
    *z = acc->z;
}

#endif // BPDB_SIMD_SCALAR_H
