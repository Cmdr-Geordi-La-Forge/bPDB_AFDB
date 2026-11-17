#ifndef SIMD_SSE_H
#define SIMD_SSE_H

#include <smmintrin.h>
#include <stdint.h>

typedef __m128i simd_vec_t;

static inline simd_vec_t simd_init(void) {
    return _mm_setzero_si128();
}

static inline void simd_decode_step(simd_vec_t* acc, const int16_t** stream_ptr, int32_t* b, int32_t* x, int32_t* y, int32_t* z) {
    const int16_t* stream = *stream_ptr;
    
    if ((stream[0] & 1) == 1) { 
        // Use unaligned load because 8-byte deltas shift the 16-byte alignment boundary
        __m128i v = _mm_loadu_si128((const __m128i*)stream); 
        __m128i v_shifted = _mm_srai_epi16(v, 1);
        
        // Extract lower 64 bits (MSB), sign-extend to 32x4, shift left 15
        __m128i msb32 = _mm_slli_epi32(_mm_cvtepi16_epi32(v_shifted), 15);
        // Shift register right 64 bits to access LSB block, sign-extend, and mask
        __m128i lsb32 = _mm_and_si128(_mm_cvtepi16_epi32(_mm_srli_si128(v_shifted, 8)), _mm_set1_epi32(0x7FFF));
        
        *acc = _mm_or_si128(msb32, lsb32);
        *stream_ptr += 8;
    } else { 
        // _mm_loadl_epi64 naturally handles unaligned 8-byte reads
        __m128i v = _mm_loadl_epi64((const __m128i*)stream);
        __m128i v_delta32 = _mm_cvtepi16_epi32(_mm_srai_epi16(v, 1));
        *acc = _mm_add_epi32(*acc, v_delta32);
        *stream_ptr += 4;
    }
    
    *b = _mm_extract_epi32(*acc, 0);
    *x = _mm_extract_epi32(*acc, 1);
    *y = _mm_extract_epi32(*acc, 2);
    *z = _mm_extract_epi32(*acc, 3);
}

#endif // SIMD_SSE_H
