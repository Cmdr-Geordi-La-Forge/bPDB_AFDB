#ifndef BPDB_COORD_CODEC_H
#define BPDB_COORD_CODEC_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "bpdb_internal.h"
#include "binary/byte_buf.h"

typedef struct {
    uint32_t    total_res;
    uint32_t    total_atoms;
    uint8_t     c_term_oxt;
    int32_t     origin[3];
    uint16_t*   plddt;
    uint8_t*    packed_coords;
    size_t      packed_len;
    EscapeAtom* escapes;
    uint32_t    escape_count;
} HotPrep;

typedef struct {
    const uint8_t*    ptr;
    const uint8_t*    end;
    const EscapeAtom* escapes;
    uint32_t          esc_pos;
    uint32_t          esc_count;
    int32_t           origin[3];
    int32_t           prev_c[3];
    int32_t           pair_d[2][3];
    uint32_t          parity;
    bool              first_res;
} CoordStream;

bool   coord_prepare_v2(const ParsedChain* chain, AfdbMetadata* meta, HotPrep* out);
void   coord_hot_prep_free(HotPrep* prep);
size_t coord_hot_encoded_size(const HotPrep* prep);
void   coord_write_v2_hot(ByteWriter* w, const AfdbMetadata* meta, const ParsedChain* chain,
                          const HotPrep* prep, uint32_t cold_offset);
bool   coord_read_v2_hot(const uint8_t* buf, size_t buf_len, BpdbV2Header* out_hdr, ChainView* out_view);

void   coord_stream_init(CoordStream* cs, const ChainView* view);
bool   coord_stream_next_residue(CoordStream* cs, const ResidueDict* rd, int num_atoms,
                                 int32_t out_xyz[30][3]);

#endif // BPDB_COORD_CODEC_H
