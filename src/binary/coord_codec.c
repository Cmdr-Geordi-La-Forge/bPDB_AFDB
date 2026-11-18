#include "binary/coord_codec.h"
#include "chem/dict.h"
#include <stdlib.h>
#include <string.h>

void parsed_chain_free(ParsedChain* chain) {
    if (!chain) return;
    free(chain->residues);
    memset(chain, 0, sizeof(*chain));
}

static void push_anchor(ByteWriter* w, int32_t b, int32_t x, int32_t y, int32_t z) {
    int16_t blk[8];
    blk[0] = (int16_t)(((uint32_t)(b >> 15) << 1) | 1U);
    blk[1] = (int16_t)((uint32_t)(x >> 15) << 1);
    blk[2] = (int16_t)((uint32_t)(y >> 15) << 1);
    blk[3] = (int16_t)((uint32_t)(z >> 15) << 1);
    blk[4] = (int16_t)((uint32_t)(b & 0x7FFF) << 1);
    blk[5] = (int16_t)((uint32_t)(x & 0x7FFF) << 1);
    blk[6] = (int16_t)((uint32_t)(y & 0x7FFF) << 1);
    blk[7] = (int16_t)((uint32_t)(z & 0x7FFF) << 1);
    bw_write(w, blk, sizeof(blk));
}

static void push_delta(ByteWriter* w, int32_t db, int32_t dx, int32_t dy, int32_t dz) {
    int16_t blk[4];
    blk[0] = (int16_t)(((uint32_t)db << 1) & ~1U);
    blk[1] = (int16_t)((uint32_t)dx << 1);
    blk[2] = (int16_t)((uint32_t)dy << 1);
    blk[3] = (int16_t)((uint32_t)dz << 1);
    bw_write(w, blk, sizeof(blk));
}

static bool delta_in_range(int32_t d) {
    return d >= -16384 && d <= 16383;
}

size_t coord_max_chain_size(const ParsedChain* chain) {
    size_t total_atoms = 0;
    for (int i = 0; i < chain->res_count; i++) {
        const ResidueDict* rd = dict_find_aa(chain->residues[i].code);
        if (rd) {
            total_atoms += (i == chain->res_count - 1) ? rd->numHeavyAtomsTerm : rd->numHeavyAtoms;
        }
    }
    size_t pad = (16 - ((size_t)chain->res_count % 16)) % 16;
    return 16 + (size_t)chain->res_count + pad + (total_atoms * 16);
}

static void serialize_coordinates(ByteWriter* w, const ParsedChain* chain) {
    int32_t prev_b = 0, prev_x = 0, prev_y = 0, prev_z = 0;
    bool first_atom = true;

    for (int i = 0; i < chain->res_count; i++) {
        const ResidueDict* rd = dict_find_aa(chain->residues[i].code);
        if (!rd) continue;

        int num_atoms = (i == chain->res_count - 1) ? rd->numHeavyAtomsTerm : rd->numHeavyAtoms;
        const ParsedResidue* res = &chain->residues[i];

        for (int a = 0; a < num_atoms; a++) {
            int32_t b = res->present[a] ? res->coords[a][0] : 0;
            int32_t x = res->present[a] ? res->coords[a][1] : 0;
            int32_t y = res->present[a] ? res->coords[a][2] : 0;
            int32_t z = res->present[a] ? res->coords[a][3] : 0;

            if (first_atom) {
                push_anchor(w, b, x, y, z);
                first_atom = false;
            } else {
                int32_t db = b - prev_b;
                int32_t dx = x - prev_x;
                int32_t dy = y - prev_y;
                int32_t dz = z - prev_z;

                if (delta_in_range(db) && delta_in_range(dx) &&
                    delta_in_range(dy) && delta_in_range(dz)) {
                    push_delta(w, db, dx, dy, dz);
                } else {
                    push_anchor(w, b, x, y, z);
                }
            }

            prev_b = b; prev_x = x; prev_y = y; prev_z = z;
        }
    }
}

void coord_write_chain(ByteWriter* w, const ParsedChain* chain) {
    uint32_t tot_res = (uint32_t)chain->res_count;
    char cid[2] = { chain->chain_id, '\0' };
    int16_t res_seq_shifted = (int16_t)(((uint32_t)chain->first_seq_id << 1) | 1U);

    bw_write(w, &tot_res, 4);
    bw_write(w, cid, 2);
    bw_write(w, &res_seq_shifted, 2);
    bw_write(w, chain->fmt_widths, 8);

    for (int i = 0; i < chain->res_count; i++) {
        bw_write_u8(w, chain->residues[i].code);
    }

    size_t pad = (16 - ((size_t)chain->res_count % 16)) % 16;
    bw_write_zeros(w, pad);

    serialize_coordinates(w, chain);
}

bool coord_read_chain_view(ByteReader* r, ChainView* out_view) {
    uint8_t hdr[16];
    if (!br_read(r, hdr, sizeof(hdr))) return false;

    uint32_t total_res = 0;
    int16_t  res_seq_shifted = 0;
    memcpy(&total_res, hdr, 4);
    memcpy(&res_seq_shifted, hdr + 6, 2);

    size_t pad = (16 - (total_res % 16)) % 16;
    if (br_remaining(r) < (size_t)total_res + pad) {
        r->error = true;
        return false;
    }

    out_view->total_res    = total_res;
    out_view->chain_id     = (char)hdr[4];
    out_view->first_seq_id = res_seq_shifted >> 1;
    memcpy(out_view->fmt_widths, hdr + 8, 8);

    out_view->seq = r->ptr;
    r->ptr += (size_t)total_res + pad;

    out_view->coord_stream = r->ptr;
    out_view->stream_end   = r->end;
    return true;
}

void coord_stream_init(CoordStream* cs, const ChainView* view) {
    cs->acc = simd_init();
    cs->ptr = view->coord_stream;
    cs->end = view->stream_end;
}

bool coord_stream_step(CoordStream* cs, int32_t* b, int32_t* x, int32_t* y, int32_t* z) {
    if (cs->ptr >= cs->end) return false;
    size_t need = (cs->ptr[0] & 1U) ? 16 : 8;
    if ((size_t)(cs->end - cs->ptr) < need) return false;
    simd_decode_step(&cs->acc, &cs->ptr, b, x, y, z);
    return true;
}
