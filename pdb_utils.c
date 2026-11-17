#include "pdb_utils.h"
#include "utils.h"
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

#define MAX_LOOP_HEADERS 64

// --- Amino Acid Code Converters ---

uint8_t res3_to_type_id(const char* src, int width) {
    char buf[3] = {' ', ' ', ' '};
    int d = 0;
    for (int i = 0; i < width && d < 3; i++) {
        if (src[i] != ' ' && src[i] != '\n' && src[i] != '\r' && src[i] != '\0') {
            buf[d++] = src[i];
        }
    }
    
    uint32_t val = (uint8_t)buf[0] | ((uint8_t)buf[1] << 8) | ((uint8_t)buf[2] << 16);
    char one_letter;
    
    switch(val) {
        case 'A' | ('L'<<8) | ('A'<<16): one_letter = 'A'; break;
        case 'A' | ('R'<<8) | ('G'<<16): one_letter = 'R'; break;
        case 'A' | ('S'<<8) | ('N'<<16): one_letter = 'N'; break;
        case 'A' | ('S'<<8) | ('P'<<16): one_letter = 'D'; break;
        case 'C' | ('Y'<<8) | ('S'<<16): one_letter = 'C'; break;
        case 'G' | ('L'<<8) | ('N'<<16): one_letter = 'Q'; break;
        case 'G' | ('L'<<8) | ('U'<<16): one_letter = 'E'; break;
        case 'G' | ('L'<<8) | ('Y'<<16): one_letter = 'G'; break;
        case 'H' | ('I'<<8) | ('S'<<16): one_letter = 'H'; break;
        case 'I' | ('L'<<8) | ('E'<<16): one_letter = 'I'; break;
        case 'L' | ('E'<<8) | ('U'<<16): one_letter = 'L'; break;
        case 'L' | ('Y'<<8) | ('S'<<16): one_letter = 'K'; break;
        case 'M' | ('E'<<8) | ('T'<<16): one_letter = 'M'; break;
        case 'P' | ('H'<<8) | ('E'<<16): one_letter = 'F'; break;
        case 'P' | ('R'<<8) | ('O'<<16): one_letter = 'P'; break;
        case 'S' | ('E'<<8) | ('R'<<16): one_letter = 'S'; break;
        case 'T' | ('H'<<8) | ('R'<<16): one_letter = 'T'; break;
        case 'T' | ('R'<<8) | ('P'<<16): one_letter = 'W'; break;
        case 'T' | ('Y'<<8) | ('R'<<16): one_letter = 'Y'; break;
        case 'V' | ('A'<<8) | ('L'<<16): one_letter = 'V'; break;
        default: one_letter = 'X'; break;
    }
    return (uint8_t)(one_letter - 0x40);
}

const char* type_id_to_res3(uint8_t resType) {
    char c = (char)(resType + 0x40); 
    switch(c) {
        case 'A': return "ALA";
        case 'R': return "ARG";
        case 'N': return "ASN";
        case 'D': return "ASP";
        case 'C': return "CYS";
        case 'Q': return "GLN";
        case 'E': return "GLU";
        case 'G': return "GLY";
        case 'H': return "HIS";
        case 'I': return "ILE";
        case 'L': return "LEU";
        case 'K': return "LYS";
        case 'M': return "MET";
        case 'F': return "PHE";
        case 'P': return "PRO";
        case 'S': return "SER";
        case 'T': return "THR";
        case 'W': return "TRP";
        case 'Y': return "TYR";
        case 'V': return "VAL";
        default:  return "UNK";
    }
}

// --- Structural Conformation Encoders & Decoders ---

uint8_t encode_struct_conf_type(const char* type_str) {
    if (type_str[0] == 'O') return TYPE_OTHER; // "OTHER"
    if (type_str[0] == 'B') return CAT_BEND;   // "BEND"
    if (type_str[0] == 'S') return CAT_STRAND; // "STRN"
    
    if (type_str[0] == 'T') {                  // "TURN"
        if (strcmp(type_str, "TURN_TY1_P") == 0)  return CAT_TURN | 1;
        if (strcmp(type_str, "TURN_TY1P_P") == 0) return CAT_TURN | 2;
        if (strcmp(type_str, "TURN_TY2_P") == 0)  return CAT_TURN | 3;
        if (strcmp(type_str, "TURN_TY2P_P") == 0) return CAT_TURN | 4;
        if (strcmp(type_str, "TURN_TY3_P") == 0)  return CAT_TURN | 5;
        if (strcmp(type_str, "TURN_TY3P_P") == 0) return CAT_TURN | 6;
        return CAT_TURN | 7;
    }
    
    if (type_str[0] == 'H') {                  // "HELX"
        if (strcmp(type_str, "HELX_RH_AL_P") == 0) return CAT_HELIX | (1 << 1) | 0;
        if (strcmp(type_str, "HELX_LH_AL_P") == 0) return CAT_HELIX | (1 << 1) | 1;
        if (strcmp(type_str, "HELX_RH_OM_P") == 0) return CAT_HELIX | (2 << 1) | 0;
        if (strcmp(type_str, "HELX_LH_OM_P") == 0) return CAT_HELIX | (2 << 1) | 1;
        if (strcmp(type_str, "HELX_RH_PI_P") == 0) return CAT_HELIX | (3 << 1) | 0;
        if (strcmp(type_str, "HELX_LH_PI_P") == 0) return CAT_HELIX | (3 << 1) | 1;
        if (strcmp(type_str, "HELX_RH_GA_P") == 0) return CAT_HELIX | (4 << 1) | 0;
        if (strcmp(type_str, "HELX_LH_GA_P") == 0) return CAT_HELIX | (4 << 1) | 1;
        if (strcmp(type_str, "HELX_RH_3T_P") == 0) return CAT_HELIX | (5 << 1) | 0;
        if (strcmp(type_str, "HELX_LH_3T_P") == 0) return CAT_HELIX | (5 << 1) | 1;
        if (strcmp(type_str, "HELX_RH_27_P") == 0) return CAT_HELIX | (9 << 1) | 0;
        if (strcmp(type_str, "HELX_LH_27_P") == 0) return CAT_HELIX | (9 << 1) | 1;
        if (strcmp(type_str, "HELX_RH_PP_P") == 0) return CAT_HELIX | (10 << 1) | 0;
        if (strcmp(type_str, "HELX_LH_PP_P") == 0) return CAT_HELIX | (10 << 1) | 1;
        return CAT_HELIX | (15 << 1); // Fallback for unspecified HELX_P
    }
    
    return TYPE_OTHER;
}

const char* decode_struct_conf_type(uint8_t type_id) {
    if (type_id == TYPE_OTHER) return "OTHER";
    
    uint8_t cat = type_id & CONF_CAT_MASK;
    uint8_t sub = type_id & CONF_SUB_MASK;
    
    if (cat == CAT_BEND) return "BEND";
    if (cat == CAT_STRAND) return "STRN";
    
    if (cat == CAT_TURN) {
        switch(sub) {
            case 1: return "TURN_TY1_P";
            case 2: return "TURN_TY1P_P";
            case 3: return "TURN_TY2_P";
            case 4: return "TURN_TY2P_P";
            case 5: return "TURN_TY3_P";
            case 6: return "TURN_TY3P_P";
            default: return "TURN_P";
        }
    }
    
    if (cat == CAT_HELIX) {
        uint8_t base = sub >> 1;
        uint8_t lh = sub & 1;
        
        switch(base) {
            case 1: return lh ? "HELX_LH_AL_P" : "HELX_RH_AL_P";
            case 2: return lh ? "HELX_LH_OM_P" : "HELX_RH_OM_P";
            case 3: return lh ? "HELX_LH_PI_P" : "HELX_RH_PI_P";
            case 4: return lh ? "HELX_LH_GA_P" : "HELX_RH_GA_P";
            case 5: return lh ? "HELX_LH_3T_P" : "HELX_RH_3T_P";
            case 9: return lh ? "HELX_LH_27_P" : "HELX_RH_27_P";
            case 10: return lh ? "HELX_LH_PP_P" : "HELX_RH_PP_P";
            default: return "HELX_P";
        }
    }
    
    return "OTHER";
}

// --- Coordinate Dictionary Lookups ---

int get_atom_idx(uint8_t one_letter_code, const char* atom_name, bool is_terminal) {
    int dict_idx = -1;
    for (size_t i = 0; i < NUM_DICT_RES; i++) {
        if (DICT[i].oneLetter == one_letter_code) { dict_idx = i; break; }
    }
    if (dict_idx == -1) return -1;

    const ResidueDict* rd = &DICT[dict_idx];
    int limit = is_terminal ? rd->numHeavyAtomsTerm : rd->numHeavyAtoms;
    
    for (int a = 0; a < limit; a++) {
        char d_clean[5];
        int dj = 0;
        for (int k = 0; k < 4; ++k) {
            if (rd->atomNames[a][k] != ' ' && rd->atomNames[a][k] != '\0') d_clean[dj++] = rd->atomNames[a][k];
        }
        d_clean[dj] = '\0';
        
        if (strcmp(atom_name, d_clean) == 0) return a;
    }
    return -1;
}

// --- Parser and Writer logic for MetaStructConf ---

size_t parse_meta_struct_conf(const char* cif_data, const char* end, MetaStructConf** out_confs) {
    const char* ptr = cif_data;
    MetaStructConf* confs = NULL;
    size_t conf_count = 0;
    size_t conf_capacity = 0;

    while (ptr < end) {
        const char* loop_start = find_token(ptr, end, "loop_");
        if (!loop_start) break;
        ptr = loop_start + 5;
        
        const char* headers[MAX_LOOP_HEADERS];
        int header_lens[MAX_LOOP_HEADERS];
        int num_headers = 0;
        
        while (ptr < end) {
            while (ptr < end && (*ptr == ' ' || *ptr == '\r' || *ptr == '\n')) ptr++;
            if (ptr >= end || *ptr != '_') break;
            
            headers[num_headers] = ptr;
            while (ptr < end && *ptr != ' ' && *ptr != '\n' && *ptr != '\r') ptr++;
            header_lens[num_headers] = ptr - headers[num_headers];
            num_headers++;
        }
        
        if (num_headers == 0 || ptr >= end) continue;

        if (header_lens[0] > 13 && strncmp(headers[0], "_struct_conf.", 13) == 0) {
            int idx_beg_comp = -1, idx_beg_seq = -1, idx_beg_chain = -1, idx_beg_icode = -1;
            int idx_end_comp = -1, idx_end_seq = -1, idx_end_chain = -1, idx_end_icode = -1;
            int idx_conf_type = -1;
            
            for (int i = 0; i < num_headers; i++) {
                const char* h = headers[i] + 13;
                int hl = header_lens[i] - 13;
                if (hl == 16 && strncmp(h, "beg_auth_comp_id", 16) == 0) idx_beg_comp = i;
                else if (hl == 15 && strncmp(h, "beg_auth_seq_id", 15) == 0) idx_beg_seq = i;
                else if (hl == 16 && strncmp(h, "beg_auth_asym_id", 16) == 0) idx_beg_chain = i;
                else if (hl == 21 && strncmp(h, "pdbx_beg_PDB_ins_code", 21) == 0) idx_beg_icode = i;
                else if (hl == 16 && strncmp(h, "end_auth_comp_id", 16) == 0) idx_end_comp = i;
                else if (hl == 15 && strncmp(h, "end_auth_seq_id", 15) == 0) idx_end_seq = i;
                else if (hl == 16 && strncmp(h, "end_auth_asym_id", 16) == 0) idx_end_chain = i;
                else if (hl == 21 && strncmp(h, "pdbx_end_PDB_ins_code", 21) == 0) idx_end_icode = i;
                else if (hl == 12 && strncmp(h, "conf_type_id", 12) == 0) idx_conf_type = i;
            }

            const char* line_end = ptr;
            while (line_end < end && *line_end != '\n') line_end++;
            int line_len = line_end - ptr;
            
            int col_off[MAX_LOOP_HEADERS] = {0}, col_wid[MAX_LOOP_HEADERS] = {0};
            int c_idx = 0, i = 0;
            
            while (i < line_len && c_idx < num_headers) {
                while (i < line_len && ptr[i] == ' ') i++; 
                if (i == line_len) break;
                
                col_off[c_idx] = i; 
                while (i < line_len && ptr[i] != ' ') i++; 
                
                if (c_idx > 0) col_wid[c_idx - 1] = col_off[c_idx] - col_off[c_idx - 1];
                c_idx++;
            }
            if (c_idx > 0) col_wid[c_idx - 1] = line_len - col_off[c_idx - 1];

            while (ptr < end && *ptr != '#' && *ptr != '_' && strncmp(ptr, "loop_", 5) != 0) {
                if (*ptr == '\n' || *ptr == '\r' || *ptr == ' ') { ptr++; continue; }
                
                if (conf_count >= conf_capacity) {
                    conf_capacity = conf_capacity == 0 ? 32 : conf_capacity * 2;
                    confs = realloc(confs, conf_capacity * sizeof(MetaStructConf));
                }
                
                MetaStructConf* conf = &confs[conf_count++];
                memset(conf, 0, sizeof(MetaStructConf));
                
                if (idx_beg_comp != -1) conf->start.resType = res3_to_type_id(ptr + col_off[idx_beg_comp], col_wid[idx_beg_comp]);
                if (idx_beg_chain != -1) extract_string_fixed(ptr + col_off[idx_beg_chain], col_wid[idx_beg_chain], conf->start.chain, 2);
                if (idx_beg_seq != -1) conf->start.resSeq = parse_number(ptr + col_off[idx_beg_seq], col_wid[idx_beg_seq]);
                if (idx_beg_icode != -1) {
                    char ic = ptr[col_off[idx_beg_icode]];
                    conf->start.iCode = (ic == '?') ? '\0' : ic;
                }
                
                if (idx_end_comp != -1) conf->end.resType = res3_to_type_id(ptr + col_off[idx_end_comp], col_wid[idx_end_comp]);
                if (idx_end_chain != -1) extract_string_fixed(ptr + col_off[idx_end_chain], col_wid[idx_end_chain], conf->end.chain, 2);
                if (idx_end_seq != -1) conf->end.resSeq = parse_number(ptr + col_off[idx_end_seq], col_wid[idx_end_seq]);
                if (idx_end_icode != -1) {
                    char ic = ptr[col_off[idx_end_icode]];
                    conf->end.iCode = (ic == '?') ? '\0' : ic;
                }
                
                conf->length = conf->end.resSeq - conf->start.resSeq + 1;
                
                if (idx_conf_type != -1) {
                    char type_buf[16] = {0};
                    extract_string_fixed(ptr + col_off[idx_conf_type], col_wid[idx_conf_type], type_buf, 15);
                    conf->type_id = encode_struct_conf_type(type_buf);
                }
                conf->parent_id = 0;
                
                while (ptr < end && *ptr != '\n') ptr++;
                if (ptr < end) ptr++;
            }
        } else {
            while (ptr < end && strncmp(ptr, "loop_", 5) != 0 && *ptr != '#') ptr++;
        }
    }
    
    *out_confs = confs;
    return conf_count;
}

void write_meta_struct_conf(FILE* out_f, const MetaStructConf* confs, size_t count) {
    if (count == 0) return;

    fprintf(out_f, "#\n");
    if (count == 1) {
        const MetaStructConf* c = &confs[0];
        const char* type_str = decode_struct_conf_type(c->type_id);
        const char* r3_s = type_id_to_res3(c->start.resType);
        const char* r3_e = type_id_to_res3(c->end.resType);
        
        fprintf(out_f, "_struct_conf.beg_auth_asym_id      %s\n", c->start.chain);
        fprintf(out_f, "_struct_conf.beg_auth_comp_id      %s\n", r3_s);
        fprintf(out_f, "_struct_conf.beg_auth_seq_id       %d\n", c->start.resSeq);
        fprintf(out_f, "_struct_conf.beg_label_asym_id     %s\n", c->start.chain);
        fprintf(out_f, "_struct_conf.beg_label_comp_id     %s\n", r3_s);
        fprintf(out_f, "_struct_conf.beg_label_seq_id      %d\n", c->start.resSeq);
        fprintf(out_f, "_struct_conf.conf_type_id          %s\n", type_str);
        fprintf(out_f, "_struct_conf.end_auth_asym_id      %s\n", c->end.chain);
        fprintf(out_f, "_struct_conf.end_auth_comp_id      %s\n", r3_e);
        fprintf(out_f, "_struct_conf.end_auth_seq_id       %d\n", c->end.resSeq);
        fprintf(out_f, "_struct_conf.end_label_asym_id     %s\n", c->end.chain);
        fprintf(out_f, "_struct_conf.end_label_comp_id     %s\n", r3_e);
        fprintf(out_f, "_struct_conf.end_label_seq_id      %d\n", c->end.resSeq);
        fprintf(out_f, "_struct_conf.id                    %s1\n", type_str);
        fprintf(out_f, "_struct_conf.pdbx_beg_PDB_ins_code %c\n", c->start.iCode ? c->start.iCode : '?');
        fprintf(out_f, "_struct_conf.pdbx_end_PDB_ins_code %c\n", c->end.iCode ? c->end.iCode : '?');
    } else {
        fprintf(out_f, "loop_\n"
                       "_struct_conf.beg_auth_asym_id\n"
                       "_struct_conf.beg_auth_comp_id\n"
                       "_struct_conf.beg_auth_seq_id\n"
                       "_struct_conf.beg_label_asym_id\n"
                       "_struct_conf.beg_label_comp_id\n"
                       "_struct_conf.beg_label_seq_id\n"
                       "_struct_conf.conf_type_id\n"
                       "_struct_conf.end_auth_asym_id\n"
                       "_struct_conf.end_auth_comp_id\n"
                       "_struct_conf.end_auth_seq_id\n"
                       "_struct_conf.end_label_asym_id\n"
                       "_struct_conf.end_label_comp_id\n"
                       "_struct_conf.end_label_seq_id\n"
                       "_struct_conf.id\n"
                       "_struct_conf.pdbx_beg_PDB_ins_code\n"
                       "_struct_conf.pdbx_end_PDB_ins_code\n");
        
        int counts[256] = {0}; 
        int w_b_chain = 1, w_b_comp = 1, w_b_seq = 1;
        int w_e_chain = 1, w_e_comp = 1, w_e_seq = 1;
        int w_type = 1, w_id = 1;

        for (size_t i = 0; i < count; i++) {
            const char* type_str = decode_struct_conf_type(confs[i].type_id);
            counts[confs[i].type_id]++;
            
            int l;
            l = strlen(confs[i].start.chain); if (l > w_b_chain) w_b_chain = l;
            l = strlen(type_id_to_res3(confs[i].start.resType)); if (l > w_b_comp) w_b_comp = l;
            l = snprintf(NULL, 0, "%d", confs[i].start.resSeq); if (l > w_b_seq) w_b_seq = l;

            l = strlen(confs[i].end.chain); if (l > w_e_chain) w_e_chain = l;
            l = strlen(type_id_to_res3(confs[i].end.resType)); if (l > w_e_comp) w_e_comp = l;
            l = snprintf(NULL, 0, "%d", confs[i].end.resSeq); if (l > w_e_seq) w_e_seq = l;

            l = strlen(type_str); if (l > w_type) w_type = l;
            l = snprintf(NULL, 0, "%s%d", type_str, counts[confs[i].type_id]); if (l > w_id) w_id = l;
        }

        memset(counts, 0, sizeof(counts)); 

        for (size_t i = 0; i < count; i++) {
            const char* type_str = decode_struct_conf_type(confs[i].type_id);
            int id_num = ++counts[confs[i].type_id];
            
            char id_str[32];
            snprintf(id_str, sizeof(id_str), "%s%d", type_str, id_num);
            
            char i_s[2] = { confs[i].start.iCode ? confs[i].start.iCode : '?', 0 };
            char i_e[2] = { confs[i].end.iCode   ? confs[i].end.iCode   : '?', 0 };
            
            const char* r3_s = type_id_to_res3(confs[i].start.resType);
            const char* r3_e = type_id_to_res3(confs[i].end.resType);

            fprintf(out_f, "%-*s %-*s %-*d %-*s %-*s %-*d %-*s %-*s %-*s %-*d %-*s %-*s %-*d %-*s %s %s \n",
                w_b_chain, confs[i].start.chain,
                w_b_comp,  r3_s,
                w_b_seq,   confs[i].start.resSeq,
                w_b_chain, confs[i].start.chain,
                w_b_comp,  r3_s,
                w_b_seq,   confs[i].start.resSeq,
                w_type,    type_str,
                w_e_chain, confs[i].end.chain,
                w_e_comp,  r3_e,
                w_e_seq,   confs[i].end.resSeq,
                w_e_chain, confs[i].end.chain,
                w_e_comp,  r3_e,
                w_e_seq,   confs[i].end.resSeq,
                w_id,      id_str,
                i_s, i_e
            );
        }
    }
    fprintf(out_f, "#\n");
}
