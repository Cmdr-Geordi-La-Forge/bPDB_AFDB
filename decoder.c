#include "bpdb.h"
#include "pdb_utils.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#if defined(__SSE4_1__)
    #include "sse.h"
#elif defined(__aarch64__) || defined(__ARM_NEON)
    #include "neon.h"
#else
    #include "scalar.h"
#endif

extern const ResidueDict DICT[];
extern const size_t NUM_DICT_RES;

static const char* ATOM37_ORDER[37] = {
    "N", "CA", "C", "CB", "O", "CG", "CG1", "CG2", "OG", "OG1", "SG", "CD",
    "CD1", "CD2", "ND1", "ND2", "OD1", "OD2", "SD", "CE", "CE1", "CE2", "CE3",
    "NE", "NE1", "NE2", "OE1", "OE2", "CH2", "NH1", "NH2", "OH", "CZ", "CZ2",
    "CZ3", "NZ", "OXT"
};

// --- 1. Coordinate Formatters ---
static void trim_string(const char* src, char* dst, int max_len) {
    int d = 0;
    for (int i = 0; i < max_len && src[i] != '\0'; i++) {
        if (src[i] != ' ') dst[d++] = src[i];
    }
    dst[d] = '\0';
}

static void format_coord(char* buf, int32_t val, int decimals) {
    uint32_t uval = (val < 0) ? -val : val; 
    if (decimals == 3) {
        sprintf(buf, "%s%u.%03u", (val < 0) ? "-" : "", uval / 1000, uval % 1000);
    } else {
        sprintf(buf, "%s%u.%02u", (val < 0) ? "-" : "", uval / 100, uval % 100);
    }
}

// --- 2. Meta Deserialization ---
extern void write_meta_struct_conf(FILE* out_f, const MetaStructConf* confs, size_t count);

static const uint8_t* deserialize_meta_section(const uint8_t* ptr, FILE* out_f) {
    uint64_t bitmask;
    memcpy(&bitmask, ptr, 8); 
    ptr += 8;

    MetaStructConf* decoded_confs = NULL;
    size_t decoded_conf_count = 0;
    
    if (bitmask & (1ULL << 40)) { 
        uint64_t struct_len;
        memcpy(&struct_len, ptr, 8); 
        ptr += 8;
        
        decoded_conf_count = (size_t)struct_len;
        size_t payload_bytes = decoded_conf_count * sizeof(MetaStructConf);
        
        decoded_confs = malloc(payload_bytes);
        memcpy(decoded_confs, ptr, payload_bytes);
        ptr += payload_bytes;
    }

    if (decoded_conf_count > 0 && decoded_confs != NULL) {
        write_meta_struct_conf(out_f, decoded_confs, decoded_conf_count);
        free(decoded_confs);
    }
    
    return ptr;
}

// --- 3. Atom Loop Helpers ---
static void print_atom_site_header(FILE* out_f) {
    fprintf(out_f, "loop_\n"
                   "_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n"
                   "_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n"
                   "_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n"
                   "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n"
                   "_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n"
                   "_atom_site.pdbx_formal_charge\n_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n"
                   "_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n_atom_site.pdbx_PDB_model_num\n"
                   "_atom_site.pdbx_sifts_xref_db_acc\n_atom_site.pdbx_sifts_xref_db_name\n"
                   "_atom_site.pdbx_sifts_xref_db_num\n_atom_site.pdbx_sifts_xref_db_res\n");
}

static void decode_atom_coordinates(const uint8_t* chain_hdr, FILE* out_f, const char* accession) {
    uint32_t total_res;
    memcpy(&total_res, chain_hdr, 4);
    char chain_id = chain_hdr[4];
    
    int16_t resSeqShifted;
    memcpy(&resSeqShifted, chain_hdr + 6, 2);
    int first_seq_id = resSeqShifted >> 1;

    uint8_t wid_atom_id = chain_hdr[8], wid_atom_name = chain_hdr[9], wid_res_seq = chain_hdr[10];
    uint8_t wid_x = chain_hdr[11], wid_y = chain_hdr[12], wid_z = chain_hdr[13];
    uint8_t wid_b = chain_hdr[14], wid_db_num = chain_hdr[15];

    uint8_t pad = (16 - (total_res % 16)) % 16;
    const uint8_t* seq = chain_hdr + 16;
    const int16_t* stream = (const int16_t*)(chain_hdr + 16 + total_res + pad);

    print_atom_site_header(out_f);

    int atom_id = 1;
    simd_vec_t simd_acc = simd_init();

    for (uint32_t i = 0; i < total_res; i++) {
        uint8_t code = seq[i];
        int dict_idx = -1;
        for (size_t d = 0; d < NUM_DICT_RES; d++) {
            if (DICT[d].oneLetter == code) { dict_idx = d; break; }
        }
        if (dict_idx == -1) continue;

        int num_atoms = (i == total_res - 1) ? DICT[dict_idx].numHeavyAtomsTerm : DICT[dict_idx].numHeavyAtoms;
        int seq_id = first_seq_id + i;

        int32_t dec_b[30] = {0}, dec_x[30] = {0}, dec_y[30] = {0}, dec_z[30] = {0};
        for (int a = 0; a < num_atoms; a++) {
            simd_decode_step(&simd_acc, &stream, &dec_b[a], &dec_x[a], &dec_y[a], &dec_z[a]);
        }

        for (int j = 0; j < 37; j++) {
            for (int a = 0; a < num_atoms; a++) {
                char clean_dict[5];
                trim_string(DICT[dict_idx].atomNames[a], clean_dict, 4);
                
                if (strcmp(clean_dict, ATOM37_ORDER[j]) == 0) {
                    char elem[3];
                    trim_string(DICT[dict_idx].elements[a], elem, 2);
                    
                    char sx[16], sy[16], sz[16], sb[16];
                    format_coord(sx, dec_x[a], 3); format_coord(sy, dec_y[a], 3);
                    format_coord(sz, dec_z[a], 3); format_coord(sb, dec_b[a], 2);
                    
                    fprintf(out_f, "ATOM %-*d%-2s%-*s. %-4s%c 1 %-*d? %-*s%-*s%-*s1.00 %-*s? %-*d%-4s%c %-*s1 %s UNP %-*d%c \n",
                            wid_atom_id, atom_id, elem, wid_atom_name, clean_dict, DICT[dict_idx].name, chain_id,
                            wid_res_seq, seq_id, wid_x, sx, wid_y, sy, wid_z, sz, wid_b, sb,
                            wid_res_seq, seq_id, DICT[dict_idx].name, chain_id, wid_atom_name, clean_dict,
                            accession, wid_db_num, seq_id, code);
                    
                    atom_id++;
                    break; 
                }
            }
        }
    }
}

// --- 4. Main Decoder Pipeline ---
void decode_binary_to_afdb(const uint8_t* buf, size_t buf_len, FILE* out_f) {
    if (buf_len < 80) return; 

    // 1. Global Header parsing
    char accession[13] = {0};
    memcpy(accession, buf + 8, 12); 
    
    uint16_t cur_mod;
    memcpy(&cur_mod, buf + 22, 2);
    int frag_num = cur_mod + 1; // Reconstruct original fragment number

    const uint8_t* ptr = buf + 64; 

    // 2. Deserialize Meta-Bitmask
    uint64_t bitmask;
    memcpy(&bitmask, ptr, 8); 
    ptr += 8;
    uint8_t elements_mask[16] = {0};
    if (bitmask & FLAG_ELEMENTS) {
        memcpy(elements_mask, ptr, 16);
        ptr += 16;
    }

    MetaStructConf* decoded_confs = NULL;
    size_t decoded_conf_count = 0;
    if (bitmask & FLAG_STRUCT) { 
        uint64_t struct_len;
        memcpy(&struct_len, ptr, 8); 
        ptr += 8;
        decoded_conf_count = (size_t)struct_len;
        size_t payload_bytes = decoded_conf_count * sizeof(MetaStructConf);
        decoded_confs = malloc(payload_bytes);
        memcpy(decoded_confs, ptr, payload_bytes);
        ptr += payload_bytes;
    }

    // 3. Write AlphaFold Headers and Elements 
    if ((bitmask & FLAG_AFDB) && (bitmask & FLAG_CIF)) {
        // Reconstruct the full AFDB ID using the extracted UniProt ID and fragment number
        fprintf(out_f, "data_AF-%s-F%d\n#\n_entry.id AF-%s-F%d\n#\n", accession, frag_num, accession, frag_num);
        
        if (bitmask & FLAG_ELEMENTS) {
            fprintf(out_f, "loop_\n_atom_type.symbol\nC \nN \nO \n");
            if (elements_mask[1] & (1 << 7)) { // Check Bit 15 of mask (Sulfur)
                fprintf(out_f, "S \n");
            }
            fprintf(out_f, "#\n");
        }

        // Print the statically defined audit/author block cleanly
        fprintf(out_f, "%s", AFDB_STATIC_META_STRING);
    }

    // 4. Write Secondary Structure
    if (decoded_conf_count > 0 && decoded_confs != NULL) {
        write_meta_struct_conf(out_f, decoded_confs, decoded_conf_count);
        free(decoded_confs);
    }

    // 5. Decode Coordinates
    decode_atom_coordinates(ptr, out_f, accession);
}
