#include "bpdb.h"
#include "pdb_utils.h"
#include "utils.h"
#include "afdb_strings.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

#if defined(__SSE4_1__)
    #include "sse.h"
#elif defined(__aarch64__) || defined(__ARM_NEON)
    #include "neon.h"
#else
    #include "scalar.h"
#endif

extern const ResidueDict DICT[];
extern const size_t NUM_DICT_RES;

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

// --- 2. Atom Loop Helpers ---
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

// --- 3. Main Decoder Pipeline ---
void decode_binary_to_afdb(const uint8_t* buf, size_t buf_len, FILE* out_f) {
    if (buf_len < 80) return; 

    char accession[13] = {0};
    memcpy(accession, buf + 8, 12); 
    
    uint16_t cur_mod;
    memcpy(&cur_mod, buf + 22, 2);
    int frag_num = cur_mod + 1;

    const uint8_t* ptr = buf + 64; 

    uint64_t bitmask;
    memcpy(&bitmask, ptr, 8); 
    ptr += 8;

    float global_plddt = 0.0f;
    if (bitmask & FLAG_RESOLUTION) {
        memcpy(&global_plddt, ptr, 4);
        ptr += 4;
    }

    char* description = NULL;
    if (bitmask & FLAG_ENTITY) {
        uint32_t dlen;
        memcpy(&dlen, ptr, 4); 
        ptr += 4;
        if (dlen > 0) {
            description = malloc(dlen + 1);
            memcpy(description, ptr, dlen);
            description[dlen] = '\0';
            ptr += dlen;
        }
    }
    if (!description) {
        description = malloc(8);
        memcpy(description, "Unknown", 8);
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

    MetaTemplate* decoded_temps = NULL;
    size_t decoded_temp_count = 0;
    if (bitmask & FLAG_TEMPLATES) {
        uint64_t t_len;
        memcpy(&t_len, ptr, 8);
        ptr += 8;
        decoded_temp_count = (size_t)t_len;
        size_t t_bytes = decoded_temp_count * sizeof(MetaTemplate);
        decoded_temps = malloc(t_bytes);
        memcpy(decoded_temps, ptr, t_bytes);
        ptr += t_bytes;
    }

    if ((bitmask & FLAG_AFDB) && (bitmask & FLAG_CIF)) {
        fprintf(out_f, "data_AF-%s-F%d\n#\n_entry.id AF-%s-F%d\n#\n", accession, frag_num, accession, frag_num);
        
        uint32_t total_res;
        char chain_id;
        const uint8_t* seq = NULL;
        int first_seq_id = 1;
        {
            memcpy(&total_res, ptr, 4);
            chain_id = ptr[4];
            
            int16_t resSeqShifted;
            memcpy(&resSeqShifted, ptr + 6, 2);
            first_seq_id = resSeqShifted >> 1;
            
            seq = ptr + 16;
        }
        
        bool present[26] = {false};
        for(uint32_t i = 0; i < total_res; i++) {
            uint8_t c = seq[i];
            if(c >= 'A' && c <= 'Z') present[c - 'A'] = true;
        }

        fprintf(out_f, "loop_\n_atom_type.symbol\nC \nN \nO \n");
        if (present['C' - 'A'] || present['M' - 'A']) {
            fprintf(out_f, "S \n");
        }
        fprintf(out_f, "#\n");

        fprintf(out_f, "%s", AFDB_AUDIT_BLOCK);
        
        fprintf(out_f, "loop_\n_chem_comp.formula\n_chem_comp.formula_weight\n_chem_comp.id\n_chem_comp.mon_nstd_flag\n_chem_comp.name\n_chem_comp.pdbx_synonyms\n_chem_comp.type\n");
        for (size_t k = 0; k < NUM_CHEM_COMP_ROWS; k++) {
            if (present[CHEM_COMP_ROWS[k].code - 'A']) {
                fprintf(out_f, "%s\n", CHEM_COMP_ROWS[k].str);
            }
        }
        fprintf(out_f, "#\n");

        fprintf(out_f, "%s", AFDB_CITATION_BLOCK);

        fprintf(out_f, "_database_2.database_code AF-%s-F%d\n", accession, frag_num);
        fprintf(out_f, "_database_2.database_id   AlphaFoldDB\n#\n");

        fprintf(out_f, "_entity.details                  ?\n");
        fprintf(out_f, "_entity.formula_weight           ?\n");
        fprintf(out_f, "_entity.id                       1\n");
        fprintf(out_f, "_entity.pdbx_description         \"%s\"\n", description);
        fprintf(out_f, "_entity.pdbx_ec                  ?\n");
        fprintf(out_f, "_entity.pdbx_fragment            ?\n");
        fprintf(out_f, "_entity.pdbx_mutation            ?\n");
        fprintf(out_f, "_entity.pdbx_number_of_molecules 1\n");
        fprintf(out_f, "_entity.src_method               man\n");
        fprintf(out_f, "_entity.type                     polymer\n#\n");
        
        fprintf(out_f, "_entity_poly.entity_id                    1\n");
        fprintf(out_f, "_entity_poly.nstd_linkage                 no\n");
        fprintf(out_f, "_entity_poly.nstd_monomer                 no\n");

        if (total_res <= 80) {
            fprintf(out_f, "_entity_poly.pdbx_seq_one_letter_code     %.*s\n", total_res, seq);
            fprintf(out_f, "_entity_poly.pdbx_seq_one_letter_code_can %.*s\n", total_res, seq);
        } else {
            fprintf(out_f, "_entity_poly.pdbx_seq_one_letter_code     \n;");
            for (uint32_t i = 0; i < total_res; i++) {
                fputc(seq[i], out_f);
                if ((i + 1) % 80 == 0 && i != total_res - 1) fputc('\n', out_f);
            }
            fprintf(out_f, "\n;\n_entity_poly.pdbx_seq_one_letter_code_can \n;");
            for (uint32_t i = 0; i < total_res; i++) {
                fputc(seq[i], out_f);
                if ((i + 1) % 80 == 0 && i != total_res - 1) fputc('\n', out_f);
            }
            fprintf(out_f, "\n;\n");
        }
        
        fprintf(out_f, "_entity_poly.pdbx_strand_id               %c\n", chain_id);
        fprintf(out_f, "_entity_poly.type                         polypeptide(L)\n#\n");
        
        // --- 1. _entity_poly_seq block ---
        fprintf(out_f, "loop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.hetero\n_entity_poly_seq.mon_id\n_entity_poly_seq.num\n");
        int num_width = snprintf(NULL, 0, "%u", total_res);
        for (uint32_t i = 0; i < total_res; i++) {
            char res3[4] = "UNK";
            for (size_t d = 0; d < NUM_DICT_RES; d++) {
                if (DICT[d].oneLetter == seq[i]) { trim_string(DICT[d].name, res3, 4); break; }
            }
            fprintf(out_f, "1 n %-3s %-*d \n", res3, num_width, i + 1);
        }
        fprintf(out_f, "#\n");

        // --- 2. Dynamic _ma_data & Static Models ---
        if (decoded_temp_count > 0) {
            fprintf(out_f, "loop_\n_ma_data.content_type\n_ma_data.id\n_ma_data.name\n");
            fprintf(out_f, "\"model coordinates\" 1 Model             \n");
            fprintf(out_f, "\"input structure\"   2 \"Input structure\" \n#\n");
        } else {
            fprintf(out_f, "_ma_data.content_type \"model coordinates\"\n");
            fprintf(out_f, "_ma_data.id           1\n");
            fprintf(out_f, "_ma_data.name         Model\n#\n");
        }
        
        fprintf(out_f, "%s", AFDB_MA_MODEL_LIST);

        // --- 3. _ma_template_ref_db_details ---
        if (decoded_temp_count > 0) {
            fprintf(out_f, "loop_\n"
                           "_ma_template_ref_db_details.db_accession_code\n"
                           "_ma_template_ref_db_details.db_name\n"
                           "_ma_template_ref_db_details.template_id\n");
            for (size_t i = 0; i < decoded_temp_count; i++) {
                fprintf(out_f, "%s %s %u \n", decoded_temps[i].accession, decoded_temps[i].db_name, decoded_temps[i].template_id);
            }
            fprintf(out_f, "#\n");
        }

        // --- 4. _ma_qa_metric_global ---
        if (bitmask & FLAG_RESOLUTION) {
            fprintf(out_f, "_ma_qa_metric_global.metric_id    1\n");
            fprintf(out_f, "_ma_qa_metric_global.metric_value %.2f\n", global_plddt);
            fprintf(out_f, "_ma_qa_metric_global.model_id     1\n");
            fprintf(out_f, "_ma_qa_metric_global.ordinal_id   1\n#\n");
        }

        // --- 5. _ma_qa_metric_local ---
        fprintf(out_f, "loop_\n"
                       "_ma_qa_metric_local.label_asym_id\n"
                       "_ma_qa_metric_local.label_comp_id\n"
                       "_ma_qa_metric_local.label_seq_id\n"
                       "_ma_qa_metric_local.metric_id\n"
                       "_ma_qa_metric_local.metric_value\n"
                       "_ma_qa_metric_local.model_id\n"
                       "_ma_qa_metric_local.ordinal_id\n");
                       
        const int16_t* c_stream = (const int16_t*)(ptr + 16 + total_res + ((16 - (total_res % 16)) % 16));
        simd_vec_t simd_acc_copy = simd_init();
        
        for (uint32_t i = 0; i < total_res; i++) {
            uint8_t code = seq[i];
            int dict_idx = -1;
            for (size_t d = 0; d < NUM_DICT_RES; d++) {
                if (DICT[d].oneLetter == code) { dict_idx = d; break; }
            }
            if (dict_idx == -1) continue;

            char res3[4] = "UNK";
            trim_string(DICT[dict_idx].name, res3, 4);
            int num_atoms = (i == total_res - 1) ? DICT[dict_idx].numHeavyAtomsTerm : DICT[dict_idx].numHeavyAtoms;
            int seq_id = first_seq_id + i;
            
            int32_t b = 0;
            for (int a = 0; a < num_atoms; a++) {
                int32_t tb, tx, ty, tz;
                simd_decode_step(&simd_acc_copy, &c_stream, &tb, &tx, &ty, &tz);
                if (a == 0) { b = tb; }
            }
            
            char sb[16];
            format_coord(sb, b, 2);
            fprintf(out_f, "%c %-3s %-*d 2 %-5s 1 %-*d \n", chain_id, res3, num_width, seq_id, sb, num_width, i + 1);
        }
        fprintf(out_f, "#\n");

        // --- 6. Static Target Entity Instance ---
        fprintf(out_f, "%s", AFDB_MA_POST_BLOCK);
    }
    
    if (decoded_temps) free(decoded_temps);
    free(description);

    if (decoded_conf_count > 0 && decoded_confs != NULL) {
        write_meta_struct_conf(out_f, decoded_confs, decoded_conf_count);
        free(decoded_confs);
    }

    decode_atom_coordinates(ptr, out_f, accession);
}
