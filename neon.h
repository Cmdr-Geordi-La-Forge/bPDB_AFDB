#ifndef SIMD_NEON_H
#define SIMD_NEON_H

#include <arm_neon.h>
#include <stdint.h>

typedef int32x4_t simd_vec_t;

static inline simd_vec_t simd_init(void) {
    return vdupq_n_s32(0);
}

static inline void simd_decode_step(simd_vec_t* acc, const int16_t** stream_ptr, int32_t* b, int32_t* x, int32_t* y, int32_t* z) {
    const int16_t* stream = *stream_ptr;
    
    if ((stream[0] & 1) == 1) { 
        // Load 8x16-bit block and arithmetic shift right by 1
        int16x8_t v = vld1q_s16(stream);
        int16x8_t v_shifted = vshrq_n_s16(v, 1);
        
        // Extract lower 4x16-bit (MSB), sign-extend to 32-bit, and shift left 15
        int32x4_t msb32 = vshlq_n_s32(vmovl_s16(vget_low_s16(v_shifted)), 15);
        // Extract upper 4x16-bit (LSB), sign-extend to 32-bit, and mask 15-bit payload
        int32x4_t lsb32 = vandq_s32(vmovl_s16(vget_high_s16(v_shifted)), vdupq_n_s32(0x7FFF));
        
        *acc = vorrq_s32(msb32, lsb32);
        *stream_ptr += 8;
    } else { 
        // Load 4x16-bit block, shift, sign-extend, and accumulate
        int16x4_t v = vld1_s16(stream);
        int32x4_t v_delta32 = vmovl_s16(vshr_n_s16(v, 1));
        *acc = vaddq_s32(*acc, v_delta32);
        *stream_ptr += 4;
    }
    
    *b = vgetq_lane_s32(*acc, 0);
    *x = vgetq_lane_s32(*acc, 1);
    *y = vgetq_lane_s32(*acc, 2);
    *z = vgetq_lane_s32(*acc, 3);
}

#endif // SIMD_NEON_H
