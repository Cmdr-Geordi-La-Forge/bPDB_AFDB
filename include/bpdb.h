#ifndef BPDB_H
#define BPDB_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

// Lossless mmCIF <-> bPDB v2
size_t encode_afdb_to_binary(const char* cif_data, size_t cif_len, uint8_t** out_buf);
void   decode_binary_to_afdb(const uint8_t* buf, size_t buf_len, FILE* out_f);

// Zero-allocation Hot Geometry API for single-core geometric analysis
typedef struct {
    uint32_t        total_res;
    uint32_t        total_atoms;
    int32_t         first_seq_id;
    char            chain_id;
    uint8_t         c_term_oxt;
    const char*     seq;           // direct pointer to ASCII 'A'..'Z' (length = total_res)
    const uint16_t* plddt;         // direct pointer to uint16_t[total_res] (scaled x100)
    char            accession[13]; // NUL-terminated UniProt accession
    uint16_t        frag_num;
} BpdbHotView;

bool bpdb_open_hot(const uint8_t* buf, size_t buf_len, BpdbHotView* out_view);

// Decodes all heavy atoms into float4 [x, y, z, 0] (out_xyzw must hold total_atoms * 4 floats).
// Uses AArch64 NEON / x86-64 SSE4.1 with a branchless fast path when escape_count == 0.
bool bpdb_decode_coords_f32(const uint8_t* buf, size_t buf_len, float* out_xyzw);

#endif // BPDB_H
