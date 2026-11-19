#include "binary/coord_codec.h"
#include "bpdb.h"
#include "chem/dict.h"
#include <stdlib.h>
#include <string.h>

#if defined(__SSE4_1__)
#include <smmintrin.h>
#define BPDB_HAVE_SSE41 1
#define BPDB_HAVE_NEON  0
#elif defined(__aarch64__)
#include <arm_neon.h>
#define BPDB_HAVE_SSE41 0
#define BPDB_HAVE_NEON  1
#else
#define BPDB_HAVE_SSE41 0
#define BPDB_HAVE_NEON  0
#endif

void parsed_chain_free(ParsedChain* chain) {
    if (!chain) return;
    free(chain->residues);
    memset(chain, 0, sizeof(*chain));
}

void coord_hot_prep_free(HotPrep* prep) {
    if (!prep) return;
    free(prep->plddt);
    free(prep->packed_coords);
    free(prep->escapes);
    memset(prep, 0, sizeof(*prep));
}

static inline size_t align4(size_t n) {
    return (n + 3U) & ~(size_t)3U;
}

static inline size_t coord_padded_size(uint32_t total_atoms) {
    if (total_atoms <= 1) return 0;
    size_t n_pairs = ((size_t)(total_atoms - 1U) + 1U) / 2U;
    // Guarantee >= 16 readable bytes from the start of the last 9-byte pair
    return align4(n_pairs * 9U + 7U);
}

static inline void pack_pair_12(uint8_t* dst, const int32_t d0[3], const int32_t d1[3]) {
    uint64_t a0 = ((uint64_t)(d0[0] & 0xFFF)) |
                  ((uint64_t)(d0[1] & 0xFFF) << 12) |
                  ((uint64_t)(d0[2] & 0xFFF) << 24);
    uint64_t a1 = ((uint64_t)(d1[0] & 0xFFF)) |
                  ((uint64_t)(d1[1] & 0xFFF) << 12) |
                  ((uint64_t)(d1[2] & 0xFFF) << 24);
    uint64_t lo = a0 | ((a1 & 0x0FFFFFFFULL) << 36);
    memcpy(dst, &lo, 8);
    dst[8] = (uint8_t)(a1 >> 28);
}

static inline int32_t sext12(uint32_t u) {
    return (int32_t)(u << 20) >> 20;
}

static inline void unpack_pair_12_scalar(const uint8_t* p, int32_t d0[3], int32_t d1[3]) {
    uint64_t v0, v1;
    memcpy(&v0, p,     8);
    memcpy(&v1, p + 4, 8);
    v1 >>= 4;
    d0[0] = sext12((uint32_t)(v0));
    d0[1] = sext12((uint32_t)(v0 >> 12));
    d0[2] = sext12((uint32_t)(v0 >> 24));
    d1[0] = sext12((uint32_t)(v1));
    d1[1] = sext12((uint32_t)(v1 >> 12));
    d1[2] = sext12((uint32_t)(v1 >> 24));
}

#if BPDB_HAVE_SSE41
static inline void unpack_pair_12_sse(const uint8_t* p, __m128i* out0, __m128i* out1) {
    const __m128i shuf = _mm_setr_epi8(
        0, 1,  1, 2,  3, 4,  -1, -1,
        4, 5,  6, 7,  7, 8,  -1, -1
    );
    const __m128i mul = _mm_setr_epi16(16, 1, 16, 0,  1, 16, 1, 0);
    __m128i raw = _mm_loadu_si128((const __m128i*)p);
    __m128i v   = _mm_srai_epi16(_mm_mullo_epi16(_mm_shuffle_epi8(raw, shuf), mul), 4);
    *out0 = _mm_cvtepi16_epi32(v);
    *out1 = _mm_cvtepi16_epi32(_mm_srli_si128(v, 8));
}
#elif BPDB_HAVE_NEON
static inline void unpack_pair_12_neon(const uint8_t* p, int32x4_t* out0, int32x4_t* out1) {
    static const uint8_t shuf_data[16] = {
        0, 1,  1, 2,  3, 4,  0xFF, 0xFF,
        4, 5,  6, 7,  7, 8,  0xFF, 0xFF
    };
    static const int16_t shift_data[8] = { 4, 0, 4, 0,  0, 4, 0, 0 };
    uint8x16_t raw = vld1q_u8(p);
    int16x8_t  v   = vreinterpretq_s16_u8(vqtbl1q_u8(raw, vld1q_u8(shuf_data)));
    v = vshrq_n_s16(vshlq_s16(v, vld1q_s16(shift_data)), 4);
    *out0 = vmovl_s16(vget_low_s16(v));
    *out1 = vmovl_high_s16(v);
}
#endif

static bool push_neg_zero(AfdbMetadata* meta, uint32_t entry) {
    if (meta->neg_zero_count >= meta->neg_zero_cap) {
        uint32_t ncap = (meta->neg_zero_cap == 0) ? 8U : (meta->neg_zero_cap * 2U);
        uint32_t* nr = realloc(meta->neg_zeros, (size_t)ncap * sizeof(uint32_t));
        if (!nr) return false;
        meta->neg_zeros = nr;
        meta->neg_zero_cap = ncap;
    }
    meta->neg_zeros[meta->neg_zero_count++] = entry;
    return true;
}

bool coord_prepare_v2(const ParsedChain* chain, AfdbMetadata* meta, HotPrep* out) {
    if (!chain || chain->res_count <= 0) return false;

    uint32_t n_res = (uint32_t)chain->res_count;
    const ParsedResidue* last_res = &chain->residues[n_res - 1];
    const ResidueDict*   last_rd  = dict_find_aa(last_res->code);
    if (!last_rd) return false;

    uint8_t c_term_oxt = (last_res->present_mask & (1U << last_rd->numHeavyAtoms)) ? 1 : 0;

    uint32_t total_atoms = 0;
    for (uint32_t r = 0; r < n_res; r++) {
        const ResidueDict* rd = dict_find_aa(chain->residues[r].code);
        if (!rd) return false;
        total_atoms += (r == n_res - 1 && c_term_oxt) ? rd->numHeavyAtomsTerm : rd->numHeavyAtoms;
    }
    if (total_atoms == 0) return false;

    size_t packed_len = coord_padded_size(total_atoms);

    if (out->plddt_cap < n_res) {
        uint16_t* np = realloc(out->plddt, (size_t)n_res * sizeof(uint16_t));
        if (!np) return false;
        out->plddt = np;
        out->plddt_cap = n_res;
    }
    if (packed_len > out->packed_cap) {
        uint8_t* nc = realloc(out->packed_coords, packed_len);
        if (!nc) return false;
        out->packed_coords = nc;
        out->packed_cap = packed_len;
    }

    out->total_res       = n_res;
    out->total_atoms     = total_atoms;
    out->c_term_oxt      = c_term_oxt;
    out->packed_len      = packed_len;
    out->escape_count    = 0;
    meta->neg_zero_count = 0;

    int32_t  prev_c[3] = {0, 0, 0};
    int32_t  pair_d[2][3];
    uint32_t parity   = 0;
    uint32_t atom_idx = 0;
    uint8_t* dst      = out->packed_coords;

    for (uint32_t r = 0; r < n_res; r++) {
        const ParsedResidue* res = &chain->residues[r];
        const ResidueDict*   rd  = dict_find_aa(res->code);
        int num_atoms = (r == n_res - 1 && c_term_oxt) ? rd->numHeavyAtomsTerm : rd->numHeavyAtoms;
        uint32_t pmask = res->present_mask;

        out->plddt[r] = res->plddt;

        int32_t cur_xyz[25][3];
        for (int a = 0; a < num_atoms; a++) {
            int32_t x, y, z;
            uint8_t lsb;
            if (__builtin_expect((pmask & (1U << a)) != 0, 1)) {
                x   = res->xyz[a][0];
                y   = res->xyz[a][1];
                z   = res->xyz[a][2];
                lsb = res->xyz_lsb[a];
            } else {
                x = y = z = 0;
                lsb = 0;
            }
            cur_xyz[a][0] = x;
            cur_xyz[a][1] = y;
            cur_xyz[a][2] = z;

            if (__builtin_expect(lsb != 0, 0)) {
                for (uint32_t axis = 0; axis < 3; axis++) {
                    if (lsb & (1U << axis)) {
                        if (!push_neg_zero(meta, (atom_idx << 2) | (axis + 1U))) {
                            return false;
                        }
                    }
                }
            }

            if (__builtin_expect(atom_idx == 0, 0)) {
                out->origin[0] = x;
                out->origin[1] = y;
                out->origin[2] = z;
            } else {
                uint8_t par = rd->parent[a];
                const int32_t* pxyz = (par == BPDB_PARENT_PREV_C) ? prev_c : cur_xyz[par];
                int32_t dx = x - pxyz[0];
                int32_t dy = y - pxyz[1];
                int32_t dz = z - pxyz[2];

                if (__builtin_expect(
                        (uint32_t)(dx + 2047) > 4094U ||
                        (uint32_t)(dy + 2047) > 4094U ||
                        (uint32_t)(dz + 2047) > 4094U, 0)) {
                    pair_d[parity][0] = -2048;
                    pair_d[parity][1] = 0;
                    pair_d[parity][2] = 0;
                    if (out->escape_count >= out->esc_cap) {
                        uint32_t ncap = (out->esc_cap == 0) ? 4U : (out->esc_cap * 2U);
                        EscapeAtom* ne = realloc(out->escapes, (size_t)ncap * sizeof(EscapeAtom));
                        if (!ne) return false;
                        out->escapes = ne;
                        out->esc_cap = ncap;
                    }
                    out->escapes[out->escape_count++] = (EscapeAtom){ x, y, z };
                } else {
                    pair_d[parity][0] = dx;
                    pair_d[parity][1] = dy;
                    pair_d[parity][2] = dz;
                }

                if (parity == 1) {
                    pack_pair_12(dst, pair_d[0], pair_d[1]);
                    dst += 9;
                    parity = 0;
                } else {
                    parity = 1;
                }
            }
            atom_idx++;
        }
        prev_c[0] = cur_xyz[2][0];
        prev_c[1] = cur_xyz[2][1];
        prev_c[2] = cur_xyz[2][2];
    }

    if (parity == 1) {
        pair_d[1][0] = 0;
        pair_d[1][1] = 0;
        pair_d[1][2] = 0;
        pack_pair_12(dst, pair_d[0], pair_d[1]);
        dst += 9;
    }
    if (packed_len > 0) {
        size_t written = (size_t)(dst - out->packed_coords);
        if (written < packed_len) {
            memset(dst, 0, packed_len - written);
        }
    }

    if (meta->neg_zero_count > 0) {
        meta->flags |= FLAG_NEG_ZERO;
    }
    return true;
}

size_t coord_hot_encoded_size(const HotPrep* prep) {
    return sizeof(BpdbV2Header)
         + align4(prep->total_res)
         + align4((size_t)prep->total_res * 2U)
         + prep->packed_len
         + (size_t)prep->escape_count * sizeof(EscapeAtom);
}

void coord_write_v2_hot(ByteWriter* w, const AfdbMetadata* meta, const ParsedChain* chain,
                        const HotPrep* prep, uint32_t cold_offset) {
    BpdbV2Header hdr;
    memset(&hdr, 0, sizeof(hdr));
    memcpy(hdr.magic, "bPDB", 4);
    hdr.version      = 2;
    hdr.chain_id     = (uint8_t)chain->chain_id;
    hdr.c_term_oxt   = prep->c_term_oxt;
    hdr.cold_codec   = 0;
    hdr.flags        = meta->flags;
    hdr.total_res    = prep->total_res;
    hdr.total_atoms  = prep->total_atoms;
    hdr.first_seq_id = chain->first_seq_id;
    hdr.escape_count = prep->escape_count;
    hdr.origin[0]    = prep->origin[0];
    hdr.origin[1]    = prep->origin[1];
    hdr.origin[2]    = prep->origin[2];
    hdr.cold_offset  = cold_offset;
    memcpy(hdr.accession, meta->accession, 12);
    hdr.frag_num     = (uint16_t)meta->frag_num;
    bw_write(w, &hdr, sizeof(hdr));

    for (uint32_t r = 0; r < prep->total_res; r++) {
        bw_write_u8(w, chain->residues[r].code);
    }
    size_t seq_pad = align4(prep->total_res) - prep->total_res;
    bw_write_zeros(w, seq_pad);

    bw_write(w, prep->plddt, (size_t)prep->total_res * 2U);
    size_t plddt_pad = align4((size_t)prep->total_res * 2U) - (size_t)prep->total_res * 2U;
    bw_write_zeros(w, plddt_pad);

    if (prep->packed_len > 0) {
        bw_write(w, prep->packed_coords, prep->packed_len);
    }
    if (prep->escape_count > 0) {
        bw_write(w, prep->escapes, (size_t)prep->escape_count * sizeof(EscapeAtom));
    }
}

bool coord_read_v2_hot(const uint8_t* buf, size_t buf_len, BpdbV2Header* out_hdr, ChainView* out_view) {
    if (!buf || buf_len < sizeof(BpdbV2Header)) return false;

    BpdbV2Header hdr;
    memcpy(&hdr, buf, sizeof(hdr));
    if (memcmp(hdr.magic, "bPDB", 4) != 0 || hdr.version != 2) return false;
    if (hdr.cold_offset < sizeof(BpdbV2Header) || hdr.cold_offset > buf_len) return false;

    size_t seq_off   = sizeof(BpdbV2Header);
    size_t seq_len   = align4(hdr.total_res);
    size_t plddt_off = seq_off + seq_len;
    size_t plddt_len = align4((size_t)hdr.total_res * 2U);
    size_t coord_off = plddt_off + plddt_len;
    size_t coord_len = coord_padded_size(hdr.total_atoms);
    size_t esc_off   = coord_off + coord_len;
    size_t esc_len   = (size_t)hdr.escape_count * sizeof(EscapeAtom);

    if (esc_off < coord_off || esc_off + esc_len != hdr.cold_offset) return false;

    if (out_hdr) *out_hdr = hdr;
    if (out_view) {
        memset(out_view, 0, sizeof(*out_view));
        out_view->total_res    = hdr.total_res;
        out_view->total_atoms  = hdr.total_atoms;
        out_view->escape_count = hdr.escape_count;
        out_view->first_seq_id = hdr.first_seq_id;
        out_view->origin[0]    = hdr.origin[0];
        out_view->origin[1]    = hdr.origin[1];
        out_view->origin[2]    = hdr.origin[2];
        out_view->chain_id     = (char)hdr.chain_id;
        out_view->c_term_oxt   = hdr.c_term_oxt;
        out_view->seq          = buf + seq_off;
        out_view->plddt_bytes  = buf + plddt_off;
        out_view->coord_stream = buf + coord_off;
        out_view->coord_len    = coord_len;
        out_view->escapes      = hdr.escape_count ? (const EscapeAtom*)(buf + esc_off) : NULL;
    }
    return true;
}

/* ===================== Streaming Residue Decoder (for mmCIF Writer) ===================== */

void coord_stream_init(CoordStream* cs, const ChainView* view) {
    memset(cs, 0, sizeof(*cs));
    cs->ptr       = view->coord_stream;
    cs->end       = view->coord_stream + view->coord_len;
    cs->escapes   = view->escapes;
    cs->esc_count = view->escape_count;
    cs->origin[0] = view->origin[0];
    cs->origin[1] = view->origin[1];
    cs->origin[2] = view->origin[2];
    cs->first_res = true;
}

bool coord_stream_next_residue(CoordStream* cs, const ResidueDict* rd, int num_atoms,
                               int32_t out_xyz[30][3]) {
    int a_start = 0;
    if (cs->first_res) {
        out_xyz[0][0] = cs->origin[0];
        out_xyz[0][1] = cs->origin[1];
        out_xyz[0][2] = cs->origin[2];
        cs->first_res = false;
        a_start = 1;
    }

    for (int a = a_start; a < num_atoms; a++) {
        if (cs->parity == 0) {
            if ((size_t)(cs->end - cs->ptr) < 16) return false;
            unpack_pair_12_scalar(cs->ptr, cs->pair_d[0], cs->pair_d[1]);
            cs->ptr += 9;
        }
        const int32_t* d = cs->pair_d[cs->parity];
        cs->parity ^= 1U;

        if (__builtin_expect(d[0] == -2048, 0)) {
            if (cs->esc_pos >= cs->esc_count) return false;
            EscapeAtom e;
            memcpy(&e, cs->escapes + cs->esc_pos++, sizeof(e));
            out_xyz[a][0] = e.x;
            out_xyz[a][1] = e.y;
            out_xyz[a][2] = e.z;
        } else {
            uint8_t par = rd->parent[a];
            const int32_t* pv = (par == BPDB_PARENT_PREV_C) ? cs->prev_c : out_xyz[par];
            out_xyz[a][0] = pv[0] + d[0];
            out_xyz[a][1] = pv[1] + d[1];
            out_xyz[a][2] = pv[2] + d[2];
        }
    }

    cs->prev_c[0] = out_xyz[2][0];
    cs->prev_c[1] = out_xyz[2][1];
    cs->prev_c[2] = out_xyz[2][2];
    return true;
}

/* ===================== Public Zero-Allocation Hot Geometry API ===================== */

bool bpdb_open_hot(const uint8_t* buf, size_t buf_len, BpdbHotView* out_view) {
    if (!out_view) return false;
    BpdbV2Header hdr;
    ChainView    cv;
    if (!coord_read_v2_hot(buf, buf_len, &hdr, &cv)) return false;

    out_view->total_res    = hdr.total_res;
    out_view->total_atoms  = hdr.total_atoms;
    out_view->first_seq_id = hdr.first_seq_id;
    out_view->chain_id     = (char)hdr.chain_id;
    out_view->c_term_oxt   = hdr.c_term_oxt;
    out_view->seq          = (const char*)cv.seq;
    out_view->plddt        = (const uint16_t*)cv.plddt_bytes;
    memcpy(out_view->accession, hdr.accession, 12);
    out_view->accession[12] = '\0';
    out_view->frag_num     = hdr.frag_num;
    return true;
}

bool bpdb_decode_coords_f32(const uint8_t* buf, size_t buf_len, float* out_xyzw) {
    if (!out_xyzw) return false;
    BpdbV2Header hdr;
    ChainView    cv;
    if (!coord_read_v2_hot(buf, buf_len, &hdr, &cv)) return false;
    if (cv.total_res == 0 || cv.total_atoms == 0) return false;

    const uint8_t* p   = cv.coord_stream;
    const uint8_t* seq = cv.seq;
    float*         out = out_xyzw;
    uint32_t       rem = cv.total_atoms;

#if BPDB_HAVE_SSE41
    const __m128 scale = _mm_set1_ps(0.001f);
    __m128i res_v[30];
    __m128i pair_v[2];
    uint32_t parity = 0;

    #define NEXT_SSE(dst) do { \
        if (parity == 0) { unpack_pair_12_sse(p, &pair_v[0], &pair_v[1]); p += 9; } \
        (dst) = pair_v[parity]; parity ^= 1U; \
    } while (0)
    #define STORE_SSE(ptr, vec) _mm_storeu_ps((ptr), _mm_mul_ps(_mm_cvtepi32_ps(vec), scale))

    if (__builtin_expect(cv.escape_count == 0, 1)) {
        // Branchless fast path (99.95% of AFDB structures)
        const ResidueDict* rd0 = dict_find_aa(seq[0]);
        if (!rd0) return false;
        int na0 = (cv.total_res == 1 && cv.c_term_oxt) ? rd0->numHeavyAtomsTerm : rd0->numHeavyAtoms;
        if (na0 < 4 || (uint32_t)na0 > rem) return false;
        rem -= (uint32_t)na0;

        res_v[0] = _mm_setr_epi32(cv.origin[0], cv.origin[1], cv.origin[2], 0);
        STORE_SSE(out, res_v[0]); out += 4;
        for (int a = 1; a < na0; a++) {
            __m128i dv; NEXT_SSE(dv);
            res_v[a] = _mm_add_epi32(res_v[rd0->parent[a]], dv);
            STORE_SSE(out, res_v[a]); out += 4;
        }
        __m128i v_prev_c = res_v[2];

        for (uint32_t r = 1; r < cv.total_res; r++) {
            const ResidueDict* rd = dict_find_aa(seq[r]);
            if (!rd) return false;
            int na = (r == cv.total_res - 1 && cv.c_term_oxt) ? rd->numHeavyAtomsTerm : rd->numHeavyAtoms;
            if (na < 4 || (uint32_t)na > rem) return false;
            rem -= (uint32_t)na;

            __m128i d0, d1, d2, d3;
            NEXT_SSE(d0); __m128i v_n  = _mm_add_epi32(v_prev_c, d0); STORE_SSE(out + 0,  v_n);
            NEXT_SSE(d1); __m128i v_ca = _mm_add_epi32(v_n,      d1); STORE_SSE(out + 4,  v_ca);
            NEXT_SSE(d2); __m128i v_c  = _mm_add_epi32(v_ca,     d2); STORE_SSE(out + 8,  v_c);
            NEXT_SSE(d3); __m128i v_o  = _mm_add_epi32(v_c,      d3); STORE_SSE(out + 12, v_o);
            out += 16;
            res_v[1] = v_ca;
            res_v[2] = v_c;
            v_prev_c = v_c;

            const uint8_t* par = rd->parent;
            for (int a = 4; a < na; a++) {
                __m128i dv; NEXT_SSE(dv);
                res_v[a] = _mm_add_epi32(res_v[par[a]], dv);
                STORE_SSE(out, res_v[a]); out += 4;
            }
        }
        #undef NEXT_SSE
        #undef STORE_SSE
        return rem == 0;
    }
    #undef NEXT_SSE
    #undef STORE_SSE
#elif BPDB_HAVE_NEON
    const float32x4_t scale = vdupq_n_f32(0.001f);
    int32x4_t res_v[30];
    int32x4_t pair_v[2];
    uint32_t  parity = 0;

    #define NEXT_NEON(dst) do { \
        if (parity == 0) { unpack_pair_12_neon(p, &pair_v[0], &pair_v[1]); p += 9; } \
        (dst) = pair_v[parity]; parity ^= 1U; \
    } while (0)
    #define STORE_NEON(ptr, vec) vst1q_f32((ptr), vmulq_f32(vcvtq_f32_s32(vec), scale))

    if (__builtin_expect(cv.escape_count == 0, 1)) {
        const ResidueDict* rd0 = dict_find_aa(seq[0]);
        if (!rd0) return false;
        int na0 = (cv.total_res == 1 && cv.c_term_oxt) ? rd0->numHeavyAtomsTerm : rd0->numHeavyAtoms;
        if (na0 < 4 || (uint32_t)na0 > rem) return false;
        rem -= (uint32_t)na0;

        int32_t init[4] = { cv.origin[0], cv.origin[1], cv.origin[2], 0 };
        res_v[0] = vld1q_s32(init);
        STORE_NEON(out, res_v[0]); out += 4;
        for (int a = 1; a < na0; a++) {
            int32x4_t dv; NEXT_NEON(dv);
            res_v[a] = vaddq_s32(res_v[rd0->parent[a]], dv);
            STORE_NEON(out, res_v[a]); out += 4;
        }
        int32x4_t v_prev_c = res_v[2];

        for (uint32_t r = 1; r < cv.total_res; r++) {
            const ResidueDict* rd = dict_find_aa(seq[r]);
            if (!rd) return false;
            int na = (r == cv.total_res - 1 && cv.c_term_oxt) ? rd->numHeavyAtomsTerm : rd->numHeavyAtoms;
            if (na < 4 || (uint32_t)na > rem) return false;
            rem -= (uint32_t)na;

            int32x4_t d0, d1, d2, d3;
            NEXT_NEON(d0); int32x4_t v_n  = vaddq_s32(v_prev_c, d0); STORE_NEON(out + 0,  v_n);
            NEXT_NEON(d1); int32x4_t v_ca = vaddq_s32(v_n,      d1); STORE_NEON(out + 4,  v_ca);
            NEXT_NEON(d2); int32x4_t v_c  = vaddq_s32(v_ca,     d2); STORE_NEON(out + 8,  v_c);
            NEXT_NEON(d3); int32x4_t v_o  = vaddq_s32(v_c,      d3); STORE_NEON(out + 12, v_o);
            out += 16;
            res_v[1] = v_ca;
            res_v[2] = v_c;
            v_prev_c = v_c;

            const uint8_t* par = rd->parent;
            for (int a = 4; a < na; a++) {
                int32x4_t dv; NEXT_NEON(dv);
                res_v[a] = vaddq_s32(res_v[par[a]], dv);
                STORE_NEON(out, res_v[a]); out += 4;
            }
        }
        #undef NEXT_NEON
        #undef STORE_NEON
        return rem == 0;
    }
    #undef NEXT_NEON
    #undef STORE_NEON
#endif

    // Escaped-structure or scalar fallback path
    CoordStream cs;
    coord_stream_init(&cs, &cv);
    for (uint32_t r = 0; r < cv.total_res; r++) {
        const ResidueDict* rd = dict_find_aa(seq[r]);
        if (!rd) return false;
        int na = (r == cv.total_res - 1 && cv.c_term_oxt) ? rd->numHeavyAtomsTerm : rd->numHeavyAtoms;
        if ((uint32_t)na > rem) return false;
        rem -= (uint32_t)na;

        int32_t xyz[30][3];
        if (!coord_stream_next_residue(&cs, rd, na, xyz)) return false;
        for (int a = 0; a < na; a++) {
            out[0] = xyz[a][0] * 0.001f;
            out[1] = xyz[a][1] * 0.001f;
            out[2] = xyz[a][2] * 0.001f;
            out[3] = 0.0f;
            out += 4;
        }
    }
    return rem == 0;
}
