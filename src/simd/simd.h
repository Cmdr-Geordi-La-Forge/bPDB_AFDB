#ifndef BPDB_SIMD_H
#define BPDB_SIMD_H

#if defined(__SSE4_1__)
    #include "simd/sse.h"
#elif defined(__aarch64__) || defined(__ARM_NEON)
    #include "simd/neon.h"
#else
    #include "simd/scalar.h"
#endif

#endif // BPDB_SIMD_H
