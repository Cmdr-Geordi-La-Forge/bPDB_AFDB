#ifndef SIMD_SCALAR_H
#define SIMD_SCALAR_H

#include <stdint.h>

typedef struct {
    int32_t b, x, y, z;
} simd_vec_t;

static inline simd_vec_t simd_init(void) {
    simd_vec_t v = {0, 0, 0, 0};
    return v;
}

static inline void simd_decode_step(simd_vec_t* acc, const int16_t** stream_ptr, int32_t* b, int32_t* x, int32_t* y, int32_t* z) {
    const int16_t* stream = *stream_ptr;
    
    if ((stream[0] & 1) == 1) { 
        acc->b = ((int32_t)(stream[0] >> 1) << 15) | ((stream[4] >> 1) & 0x7FFF);
        acc->x = ((int32_t)(stream[1] >> 1) << 15) | ((stream[5] >> 1) & 0x7FFF);
        acc->y = ((int32_t)(stream[2] >> 1) << 15) | ((stream[6] >> 1) & 0x7FFF);
        acc->z = ((int32_t)(stream[3] >> 1) << 15) | ((stream[7] >> 1) & 0x7FFF);
        *stream_ptr += 8;
    } else { 
        acc->b += (stream[0] >> 1);
        acc->x += (stream[1] >> 1);
        acc->y += (stream[2] >> 1);
        acc->z += (stream[3] >> 1);
        *stream_ptr += 4;
    }
    
    *b = acc->b; *x = acc->x; *y = acc->y; *z = acc->z;
}

#endif // SIMD_SCALAR_H
