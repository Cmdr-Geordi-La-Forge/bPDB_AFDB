#include "bpdb.h"
#include "utils.h"
#include "pdb_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define MAX_LOOP_HEADERS 64
#define IDX_ATOM_ID   1
#define IDX_ATOM_NAME 3
#define IDX_CHAIN     6
#define IDX_RES_SEQ   8
#define IDX_X         10
#define IDX_Y         11
#define IDX_Z         12
#define IDX_B         14
#define IDX_DB_NUM    23

typedef struct {
    uint8_t code;
    int32_t coords[30][4];
    uint8_t present[30];
} ParsedResidue;

// --- 1. Coordinate Stream Writers ---
static void push_anchor(int16_t** stream_ptr, int32_t b, int32_t x, int32_t y, int32_t z) {
    int16_t* p = *stream_ptr;
    p[0] = (int16_t)(((b >> 15) << 1) | 1);
    p[1] = (int16_t)((x >> 15) << 1);
    p[2] = (int16_t)((y >> 15) << 1);
    p[3] = (int16_t)((z >> 15) << 1);
    p[4] = (int16_t)((b & 0x7FFF) << 1);
    p[5] = (int16_t)((x & 0x7FFF) << 1);
    p[6] = (int16_t)((y & 0x7FFF) << 1);
    p[7] = (int16_t)((z & 0x7FFF) << 1);
    *stream_ptr = p + 8;
}

static void push_delta(int16_t** stream_ptr, int32_t db, int32_t dx, int32_t dy, int32_t dz) {
    int16_t* p = *stream_ptr;
    p[0] = (int16_t)((db << 1) & ~1);
    p[1] = (int16_t)(dx << 1);
    p[2] = (int16_t)(dy << 1);
    p[3] = (int16_t)(dz << 1);
    *stream_ptr = p + 4;
}

// --- 2. Extraction Helpers ---
static void parse_global_ids(const char* cif_data, const char* end, char* accession, int* frag_num) {
    const char* data_hdr = find_token(cif_data, end, "data_AF-");
    if (data_hdr) {
        data_hdr += 8;
        const char* dash = memchr(data_hdr, '-', end - data_hdr);
        if (dash) {
            int acc_len = dash - data_hdr;
            if (acc_len > 12) acc_len = 12;
            memcpy(accession, data_hdr, acc_len);
            accession[acc_len] = '\0';
            if (dash[1] == 'F') {
                *frag_num = atoi(dash + 2);
            }
        }
    }
}

static void parse_entity_desc(const char* cif_data, const char* end, char** out_desc) {
    const char* ptr = find_token(cif_data, end, "_entity.pdbx_description");
    if (!ptr) {
        *out_desc = NULL;
        return;
    }
    ptr += strlen("_entity.pdbx_description");
    while (ptr < end && (*ptr == ' ' || *ptr == '\t')) ptr++;
    
    char quote = 0;
    if (*ptr == '"' || *ptr == '\'') {
        quote = *ptr++;
    }
    
    const char* start = ptr;
    while (ptr < end) {
        if (quote) {
            if (*ptr == quote) break;
        } else {
            if (*ptr == ' ' || *ptr == '\n' || *ptr == '\r') break;
        }
        ptr++;
    }
    
    int len = ptr - start;
    if (len > 0) {
        *out_desc = malloc(len + 1);
        memcpy(*out_desc, start, len);
        (*out_desc)[len] = '\0';
    } else {
        *out_desc = NULL;
    }
}

static int parse_entity_seq(const char* cif_data, const char* end, ParsedResidue** out_residues) {
    const char* seq_ptr = find_token(cif_data, end, "_entity_poly.pdbx_seq_one_letter_code");
    if (!seq_ptr) return 0;
    
    seq_ptr += strlen("_entity_poly.pdbx_seq_one_letter_code");
    while (seq_ptr < end && (*seq_ptr == ' ' || *seq_ptr == '\t')) seq_ptr++;
    
    bool is_semi = false;
    if (seq_ptr < end && *seq_ptr == '\n') {
        seq_ptr++;
        if (seq_ptr < end && *seq_ptr == ';') {
            is_semi = true;
            seq_ptr++;
        }
    } else if (seq_ptr < end && *seq_ptr == ';') {
        is_semi = true;
        seq_ptr++;
    }

    int capacity = 1024;
    int res_count = 0;
    ParsedResidue* residues = calloc(capacity, sizeof(ParsedResidue));
    
    while (seq_ptr < end) {
        if (is_semi) {
            if (*seq_ptr == ';' && *(seq_ptr - 1) == '\n') break; 
        } else {
            if (*seq_ptr == ' ' || *seq_ptr == '\n' || *seq_ptr == '\r') break;
        }
        
        char c = *seq_ptr++;
        if (c >= 'A' && c <= 'Z') {
            if (c == 'U') c = 'C';
            if (c == 'O') c = 'X';
            
            if (res_count >= capacity) {
                capacity *= 2;
                residues = realloc(residues, capacity * sizeof(ParsedResidue));
                memset(residues + res_count, 0, (capacity - res_count) * sizeof(ParsedResidue));
            }
            residues[res_count].code = c;
            res_count++;
        }
    }
    *out_residues = residues;
    return res_count;
}

static void parse_cif_atoms(const char* cif_data, const char* end, ParsedResidue* residues, int res_count,
                            int* out_first_seq_id, char* out_chain_id, uint8_t* fmt_widths) {
    bool offsets_determined = false;
    int off[30] = {0}, wid[30] = {0};
    
    const char* line = find_token(cif_data, end, "ATOM ");
    while (line && line < end) {
        const char* eol = memchr(line, '\n', end - line);
        if (!eol) eol = end;
        int len = eol - line;

        if (len >= 5 && strncmp(line, "ATOM ", 5) == 0) {
            if (!offsets_determined) {
                int tok_starts[30] = {0};
                int tok_count = 0, i = 0;
                
                while (i < len && tok_count < 30) {
                    while (i < len && line[i] == ' ') i++; 
                    if (i == len) break;
                    tok_starts[tok_count++] = i;
                    while (i < len && line[i] != ' ') i++; 
                }

                if (tok_count > IDX_DB_NUM) {
                    wid[IDX_ATOM_ID]   = tok_starts[IDX_ATOM_ID+1]   - tok_starts[IDX_ATOM_ID];
                    off[IDX_ATOM_NAME] = tok_starts[IDX_ATOM_NAME]; wid[IDX_ATOM_NAME] = tok_starts[IDX_ATOM_NAME+1] - tok_starts[IDX_ATOM_NAME];
                    off[IDX_CHAIN]     = tok_starts[IDX_CHAIN];     wid[IDX_CHAIN]     = tok_starts[IDX_CHAIN+1]     - tok_starts[IDX_CHAIN];
                    off[IDX_RES_SEQ]   = tok_starts[IDX_RES_SEQ];   wid[IDX_RES_SEQ]   = tok_starts[IDX_RES_SEQ+1]   - tok_starts[IDX_RES_SEQ];
                    off[IDX_X]         = tok_starts[IDX_X];         wid[IDX_X]         = tok_starts[IDX_X+1]         - tok_starts[IDX_X];
                    off[IDX_Y]         = tok_starts[IDX_Y];         wid[IDX_Y]         = tok_starts[IDX_Y+1]         - tok_starts[IDX_Y];
                    off[IDX_Z]         = tok_starts[IDX_Z];         wid[IDX_Z]         = tok_starts[IDX_Z+1]         - tok_starts[IDX_Z];
                    off[IDX_B]         = tok_starts[IDX_B];         wid[IDX_B]         = tok_starts[IDX_B+1]         - tok_starts[IDX_B];
                    wid[IDX_DB_NUM]    = tok_starts[IDX_DB_NUM+1]   - tok_starts[IDX_DB_NUM];
                    
                    fmt_widths[0] = wid[IDX_ATOM_ID]; fmt_widths[1] = wid[IDX_ATOM_NAME]; fmt_widths[2] = wid[IDX_RES_SEQ];
                    fmt_widths[3] = wid[IDX_X];       fmt_widths[4] = wid[IDX_Y];         fmt_widths[5] = wid[IDX_Z];
                    fmt_widths[6] = wid[IDX_B];       fmt_widths[7] = wid[IDX_DB_NUM];
                    offsets_determined = true;
                }
            }

            if (offsets_determined) {
                int seq_id = parse_number(line + off[IDX_RES_SEQ], wid[IDX_RES_SEQ]);
                if (*out_first_seq_id == -1) {
                    *out_first_seq_id = seq_id;
                    char chain_str[3] = {0};
                    extract_string_fixed(line + off[IDX_CHAIN], wid[IDX_CHAIN], chain_str, 2);
                    *out_chain_id = chain_str[0];
                }

                int res_idx = seq_id - 1; 
                if (res_idx >= 0 && res_idx < res_count) {
                    char atom_name[5] = {0};
                    extract_string_fixed(line + off[IDX_ATOM_NAME], wid[IDX_ATOM_NAME], atom_name, 4);

                    int atom_idx = get_atom_idx(residues[res_idx].code, atom_name, res_idx == res_count - 1);
                    if (atom_idx != -1) {
                        residues[res_idx].coords[atom_idx][0] = parse_number(line + off[IDX_B], wid[IDX_B]);
                        residues[res_idx].coords[atom_idx][1] = parse_number(line + off[IDX_X], wid[IDX_X]);
                        residues[res_idx].coords[atom_idx][2] = parse_number(line + off[IDX_Y], wid[IDX_Y]);
                        residues[res_idx].coords[atom_idx][3] = parse_number(line + off[IDX_Z], wid[IDX_Z]);
                        residues[res_idx].present[atom_idx] = 1;
                    }
                }
            }
        }
        line = eol + 1;
    }
}

static size_t serialize_coordinates(ParsedResidue* residues, int res_count, int16_t* stream) {
    int32_t prev_b = 0, prev_x = 0, prev_y = 0, prev_z = 0;
    int first_atom = 1;
    int16_t* start_stream = stream;

    for (int i = 0; i < res_count; i++) {
        int dict_idx = -1;
        for (size_t d = 0; d < NUM_DICT_RES; d++) {
            if (DICT[d].oneLetter == residues[i].code) { dict_idx = d; break; }
        }
        if (dict_idx == -1) continue;

        int num_atoms = (i == res_count - 1) ? DICT[dict_idx].numHeavyAtomsTerm : DICT[dict_idx].numHeavyAtoms;
        
        for (int a = 0; a < num_atoms; a++) {
            int32_t b = residues[i].present[a] ? residues[i].coords[a][0] : 0;
            int32_t x = residues[i].present[a] ? residues[i].coords[a][1] : 0;
            int32_t y = residues[i].present[a] ? residues[i].coords[a][2] : 0;
            int32_t z = residues[i].present[a] ? residues[i].coords[a][3] : 0;

            if (first_atom) {
                push_anchor(&stream, b, x, y, z);
                prev_b = b; prev_x = x; prev_y = y; prev_z = z;
                first_atom = 0;
                continue;
            }

            int32_t db = b - prev_b, dx = x - prev_x, dy = y - prev_y, dz = z - prev_z;

            if (db < -16384 || db > 16383 || dx < -16384 || dx > 16383 ||
                dy < -16384 || dy > 16383 || dz < -16384 || dz > 16383) {
                push_anchor(&stream, b, x, y, z);
            } else {
                push_delta(&stream, db, dx, dy, dz);
            }

            prev_b = b; prev_x = x; prev_y = y; prev_z = z;
        }
    }
    return (size_t)((uint8_t*)stream - (uint8_t*)start_stream);
}

// --- Template Parser Helper ---
static size_t parse_meta_templates(const char* cif_data, const char* end, MetaTemplate** out_temps) {
    const char* ptr = find_token(cif_data, end, "_ma_template_ref_db_details.template_id");
    if (!ptr) return 0;
    
    while (ptr < end && *ptr != '\n') ptr++;
    if (ptr < end) ptr++;
    
    size_t cap = 4, count = 0;
    MetaTemplate* temps = malloc(cap * sizeof(MetaTemplate));
    
    while (ptr < end && *ptr != '#' && *ptr != '_' && strncmp(ptr, "loop_", 5) != 0) {
        if (*ptr == ' ' || *ptr == '\n' || *ptr == '\r') { ptr++; continue; }
        
        if (count >= cap) { cap *= 2; temps = realloc(temps, cap * sizeof(MetaTemplate)); }
        memset(&temps[count], 0, sizeof(MetaTemplate));
        
        const char* start = ptr;
        while (ptr < end && *ptr != ' ') ptr++;
        int len = ptr - start;
        if (len > 11) len = 11;
        memcpy(temps[count].accession, start, len);
        
        while (ptr < end && *ptr == ' ') ptr++;
        
        start = ptr;
        while (ptr < end && *ptr != ' ') ptr++;
        len = ptr - start;
        if (len > 7) len = 7;
        memcpy(temps[count].db_name, start, len);
        
        while (ptr < end && *ptr == ' ') ptr++;
        
        start = ptr;
        while (ptr < end && *ptr != ' ' && *ptr != '\n' && *ptr != '\r') ptr++;
        temps[count].template_id = parse_number(start, ptr - start);
        
        count++;
        while (ptr < end && *ptr != '\n') ptr++;
    }
    *out_temps = temps;
    return count;
}

extern size_t parse_meta_struct_conf(const char* cif_data, const char* end, MetaStructConf** out_confs);

// --- 3. Main Encoder Pipeline ---
size_t encode_afdb_to_binary(const char* cif_data, size_t cif_len, uint8_t** out_buf) {
    const char* end = cif_data + cif_len;
    
    uint64_t bitmask = FLAG_AF_MODE | FLAG_AFDB | FLAG_CIF;
    size_t meta_size = 8; 

    // 1. Base CIF info
    char accession[13] = {0};
    int frag_num = 1;
    parse_global_ids(cif_data, end, accession, &frag_num);

    char* desc = NULL;
    parse_entity_desc(cif_data, end, &desc);

    ParsedResidue* residues = NULL;
    int res_count = parse_entity_seq(cif_data, end, &residues);
    if (res_count <= 0) {
        if (desc) free(desc);
        *out_buf = NULL; return 0;
    }

    // 2. Entity Desc Bitmask
    bitmask |= FLAG_ENTITY;
    size_t desc_len = desc ? strlen(desc) : 0;
    meta_size += 4 + desc_len; 

    // 3. Meta Extraction (Struct Confs & Templates)
    MetaStructConf* struct_confs = NULL;
    size_t struct_conf_count = parse_meta_struct_conf(cif_data, end, &struct_confs);
    if (struct_conf_count > 0) {
        bitmask |= FLAG_STRUCT; 
        meta_size += 8 + (struct_conf_count * sizeof(MetaStructConf)); 
    }

    MetaTemplate* struct_temps = NULL;
    size_t struct_temp_count = parse_meta_templates(cif_data, end, &struct_temps);
    if (struct_temp_count > 0) {
        bitmask |= FLAG_TEMPLATES;
        meta_size += 8 + (struct_temp_count * sizeof(MetaTemplate));
    }

    // 4. Atom Extraction
    int first_seq_id = -1;
    char chain_id = 'A';
    uint8_t fmt_widths[8] = {0};
    parse_cif_atoms(cif_data, end, residues, res_count, &first_seq_id, &chain_id, fmt_widths);

    int total_atoms = 0;
    for (int i = 0; i < res_count; i++) {
        for (size_t d = 0; d < NUM_DICT_RES; d++) {
            if (DICT[d].oneLetter == residues[i].code) {
                total_atoms += (i == res_count - 1) ? DICT[d].numHeavyAtomsTerm : DICT[d].numHeavyAtoms;
                break;
            }
        }
    }

    // 5. Memory Allocation
    uint8_t pad = (16 - (res_count % 16)) % 16;
    size_t max_out_size = 64 + meta_size + 16 + res_count + pad + (total_atoms * 16); 
    *out_buf = calloc(1, max_out_size);
    uint8_t* out = *out_buf;

    // 6. Write Global Header
    memset(out, 0, 64);
    memcpy(out, "bPDB", 4); 
    out[4] = 0x00; 
    out[5] = 1;    
    out[6] = 1;    
    out[7] = 0;    
    
    memcpy(out + 8, accession, 12); 
    
    uint16_t tot_mod = 0; 
    uint16_t cur_mod = (frag_num > 0) ? (frag_num - 1) : 0;
    
    memcpy(out + 20, &tot_mod, 2);
    memcpy(out + 22, &cur_mod, 2);
    
    out += 64;

    // 7. Write Meta Payloads
    memcpy(out, &bitmask, 8); 
    out += 8;

    if (bitmask & FLAG_ENTITY) {
        uint32_t dlen = desc_len;
        memcpy(out, &dlen, 4); out += 4;
        if (dlen > 0) {
            memcpy(out, desc, dlen); out += dlen;
        }
    }
    if (desc) free(desc);

    if (bitmask & FLAG_STRUCT) {
        uint64_t struct_len = struct_conf_count;
        memcpy(out, &struct_len, 8); 
        out += 8;
        size_t payload_bytes = struct_conf_count * sizeof(MetaStructConf);
        memcpy(out, struct_confs, payload_bytes);
        out += payload_bytes;
    }
    if (struct_confs) free(struct_confs);

    if (bitmask & FLAG_TEMPLATES) {
        uint64_t t_len = struct_temp_count;
        memcpy(out, &t_len, 8); 
        out += 8;
        size_t t_bytes = struct_temp_count * sizeof(MetaTemplate);
        memcpy(out, struct_temps, t_bytes);
        out += t_bytes;
    }
    if (struct_temps) free(struct_temps);

    // 8. Write Chain Header
    uint32_t tot_res = res_count;
    memcpy(out, &tot_res, 4); out += 4;
    char cid[2] = { chain_id, '\0' };
    memcpy(out, cid, 2); out += 2;
    int16_t resSeqShifted = (first_seq_id << 1) | 1;
    memcpy(out, &resSeqShifted, 2); out += 2;
    memcpy(out, fmt_widths, 8); out += 8;

    for (int i = 0; i < res_count; i++) *out++ = residues[i].code;
    for (int i = 0; i < pad; i++) *out++ = 0x00;

    // 9. Serialize and return
    size_t stream_size = serialize_coordinates(residues, res_count, (int16_t*)out);
    
    free(residues);
    return (out - *out_buf) + stream_size;
}
