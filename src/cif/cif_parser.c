#include "cif/cif_parser.h"
#include "cif/cif_lexer.h"
#include "chem/dict.h"
#include "chem/struct_conf.h"
#include "util/str_util.h"
#include <stdlib.h>
#include <string.h>

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

static void parse_global_ids(const char* data, const char* end, char* accession, int* frag_num) {
    const char* hdr = find_token(data, end, "data_AF-");
    if (!hdr) return;

    hdr += 8;
    const char* dash = memchr(hdr, '-', (size_t)(end - hdr));
    if (!dash) return;

    size_t acc_len = (size_t)(dash - hdr);
    if (acc_len > 12) acc_len = 12;
    memcpy(accession, hdr, acc_len);
    accession[acc_len] = '\0';

    if (dash + 2 < end && dash[1] == 'F') {
        *frag_num = atoi(dash + 2);
    }
}

static float parse_global_plddt(const char* data, const char* end) {
    const char* ptr = find_token(data, end, "_ma_qa_metric_global.metric_value");
    if (!ptr) return -1.0f;
    ptr = cif_skip_spaces(ptr + strlen("_ma_qa_metric_global.metric_value"), end);
    return (ptr < end) ? strtof(ptr, NULL) : -1.0f;
}

static int parse_entity_seq(const char* data, const char* end, ParsedResidue** out_residues) {
    const char* ptr = find_token(data, end, "_entity_poly.pdbx_seq_one_letter_code");
    if (!ptr) return 0;

    ptr = cif_skip_spaces(ptr + strlen("_entity_poly.pdbx_seq_one_letter_code"), end);

    bool is_semi = false;
    if (ptr < end && *ptr == '\n') {
        ptr++;
        if (ptr < end && *ptr == ';') { is_semi = true; ptr++; }
    } else if (ptr < end && *ptr == ';') {
        is_semi = true;
        ptr++;
    }

    int capacity = 1024;
    int count = 0;
    ParsedResidue* residues = calloc((size_t)capacity, sizeof(ParsedResidue));
    if (!residues) return 0;

    while (ptr < end) {
        if (is_semi) {
            if (*ptr == ';' && *(ptr - 1) == '\n') break;
        } else {
            if (*ptr == ' ' || *ptr == '\n' || *ptr == '\r') break;
        }

        char c = *ptr++;
        if (c >= 'A' && c <= 'Z') {
            if (c == 'U') c = 'C';
            if (c == 'O') c = 'X';

            if (count >= capacity) {
                capacity *= 2;
                residues = realloc(residues, (size_t)capacity * sizeof(ParsedResidue));
                memset(residues + count, 0, (size_t)(capacity - count) * sizeof(ParsedResidue));
            }
            residues[count++].code = (uint8_t)c;
        }
    }

    *out_residues = residues;
    return count;
}

static bool detect_atom_columns(const char* line, int len, int* off, int* wid, uint8_t* fmt_widths) {
    int tok_starts[30] = {0};
    int tok_count = 0;
    int i = 0;

    while (i < len && tok_count < 30) {
        while (i < len && line[i] == ' ') i++;
        if (i == len) break;
        tok_starts[tok_count++] = i;
        while (i < len && line[i] != ' ') i++;
    }
    if (tok_count <= IDX_DB_NUM) return false;

    static const int cols[] = {
        IDX_ATOM_ID, IDX_ATOM_NAME, IDX_CHAIN, IDX_RES_SEQ,
        IDX_X, IDX_Y, IDX_Z, IDX_B, IDX_DB_NUM
    };
    for (size_t k = 0; k < sizeof(cols) / sizeof(cols[0]); k++) {
        int c = cols[k];
        off[c] = tok_starts[c];
        wid[c] = tok_starts[c + 1] - tok_starts[c];
    }

    fmt_widths[0] = (uint8_t)wid[IDX_ATOM_ID];
    fmt_widths[1] = (uint8_t)wid[IDX_ATOM_NAME];
    fmt_widths[2] = (uint8_t)wid[IDX_RES_SEQ];
    fmt_widths[3] = (uint8_t)wid[IDX_X];
    fmt_widths[4] = (uint8_t)wid[IDX_Y];
    fmt_widths[5] = (uint8_t)wid[IDX_Z];
    fmt_widths[6] = (uint8_t)wid[IDX_B];
    fmt_widths[7] = (uint8_t)wid[IDX_DB_NUM];
    return true;
}

static bool is_neg_zero_token(const char* s, int width, int32_t parsed_val) {
    if (parsed_val != 0) return false;
    for (int i = 0; i < width; i++) {
        if (s[i] == '-') return true;
        if (s[i] != ' ' && s[i] != '\t') break;
    }
    return false;
}

static void parse_atoms(const char* data, const char* end, ParsedChain* chain) {
    bool layout_ready = false;
    int off[30] = {0}, wid[30] = {0};

    const char* line = find_token(data, end, "ATOM ");
    while (line && line < end) {
        const char* eol = memchr(line, '\n', (size_t)(end - line));
        if (!eol) eol = end;
        int len = (int)(eol - line);

        if (len >= 5 && strncmp(line, "ATOM ", 5) == 0) {
            if (!layout_ready) {
                layout_ready = detect_atom_columns(line, len, off, wid, chain->fmt_widths);
            }
            if (layout_ready) {
                int seq_id = parse_number(line + off[IDX_RES_SEQ], wid[IDX_RES_SEQ]);
                if (chain->first_seq_id == -1) {
                    chain->first_seq_id = seq_id;
                    char chain_str[3] = {0};
                    strip_whitespace_fixed(line + off[IDX_CHAIN], wid[IDX_CHAIN], chain_str, 2);
                    chain->chain_id = chain_str[0];
                }

                int res_idx = seq_id - 1;
                if (res_idx >= 0 && res_idx < chain->res_count) {
                    char atom_name[5] = {0};
                    strip_whitespace_fixed(line + off[IDX_ATOM_NAME], wid[IDX_ATOM_NAME], atom_name, 4);

                    ParsedResidue* res = &chain->residues[res_idx];
                    int atom_idx = dict_get_atom_idx(res->code, atom_name, res_idx == chain->res_count - 1);
                    if (atom_idx != -1) {
                        int32_t vb = parse_number(line + off[IDX_B], wid[IDX_B]);
                        int32_t vx = parse_number(line + off[IDX_X], wid[IDX_X]);
                        int32_t vy = parse_number(line + off[IDX_Y], wid[IDX_Y]);
                        int32_t vz = parse_number(line + off[IDX_Z], wid[IDX_Z]);

                        res->coords[atom_idx][0] = vb;
                        res->coords[atom_idx][1] = vx;
                        res->coords[atom_idx][2] = vy;
                        res->coords[atom_idx][3] = vz;
                        res->present[atom_idx]   = 1;

                        uint8_t lsb = 0;
                        if (is_neg_zero_token(line + off[IDX_X], wid[IDX_X], vx)) lsb |= (1U << 0);
                        if (is_neg_zero_token(line + off[IDX_Y], wid[IDX_Y], vy)) lsb |= (1U << 1);
                        if (is_neg_zero_token(line + off[IDX_Z], wid[IDX_Z], vz)) lsb |= (1U << 2);
                        res->xyz_lsb[atom_idx] = lsb;
                    }
                }
            }
        }
        if (eol == end) break;
        line = eol + 1;
    }
}

typedef struct {
    int beg_comp, beg_seq, beg_chain, beg_icode;
    int end_comp, end_seq, end_chain, end_icode;
    int conf_type;
} StructConfCols;

static StructConfCols map_struct_conf_cols(const char** headers, const int* lens, int num_headers) {
    StructConfCols c = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
    for (int i = 0; i < num_headers; i++) {
        const char* h = headers[i] + 13;
        int hl = lens[i] - 13;
        if      (hl == 16 && strncmp(h, "beg_auth_comp_id", 16) == 0)      c.beg_comp  = i;
        else if (hl == 15 && strncmp(h, "beg_auth_seq_id", 15) == 0)       c.beg_seq   = i;
        else if (hl == 16 && strncmp(h, "beg_auth_asym_id", 16) == 0)      c.beg_chain = i;
        else if (hl == 21 && strncmp(h, "pdbx_beg_PDB_ins_code", 21) == 0) c.beg_icode = i;
        else if (hl == 16 && strncmp(h, "end_auth_comp_id", 16) == 0)      c.end_comp  = i;
        else if (hl == 15 && strncmp(h, "end_auth_seq_id", 15) == 0)       c.end_seq   = i;
        else if (hl == 16 && strncmp(h, "end_auth_asym_id", 16) == 0)      c.end_chain = i;
        else if (hl == 21 && strncmp(h, "pdbx_end_PDB_ins_code", 21) == 0) c.end_icode = i;
        else if (hl == 12 && strncmp(h, "conf_type_id", 12) == 0)          c.conf_type = i;
    }
    return c;
}

static void detect_loop_col_layout(const char* ptr, const char* end, int num_headers, int* col_off, int* col_wid) {
    const char* line_end = ptr;
    while (line_end < end && *line_end != '\n') line_end++;
    int line_len = (int)(line_end - ptr);

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
}

static size_t parse_struct_conf_single(const char* data, const char* end, MetaStructConf** out_confs) {
    const char* key = "_struct_conf.beg_auth_asym_id";
    const char* start = find_token(data, end, key);
    if (!start || !cif_has_inline_value(start + strlen(key), end)) return 0;

    const char* blk_end = cif_find_block_end(start, end);
    MetaStructConf* conf = calloc(1, sizeof(MetaStructConf));
    if (!conf) return 0;

    char buf[32];
    buf[0] = '\0'; cif_extract_kv_string(start, blk_end, "_struct_conf.beg_auth_comp_id", buf, sizeof(buf));
    if (buf[0]) conf->start.resType = res3_to_type_id(buf, (int)strlen(buf));

    buf[0] = '\0'; cif_extract_kv_string(start, blk_end, "_struct_conf.beg_auth_asym_id", buf, sizeof(buf));
    strip_whitespace_fixed(buf, (int)strlen(buf), conf->start.chain, 2);

    buf[0] = '\0'; cif_extract_kv_string(start, blk_end, "_struct_conf.beg_auth_seq_id", buf, sizeof(buf));
    conf->start.resSeq = (uint32_t)atoi(buf);

    buf[0] = '\0'; cif_extract_kv_string(start, blk_end, "_struct_conf.pdbx_beg_PDB_ins_code", buf, sizeof(buf));
    conf->start.iCode = (buf[0] == '?' || buf[0] == '\0') ? '\0' : buf[0];

    buf[0] = '\0'; cif_extract_kv_string(start, blk_end, "_struct_conf.end_auth_comp_id", buf, sizeof(buf));
    if (buf[0]) conf->end.resType = res3_to_type_id(buf, (int)strlen(buf));

    buf[0] = '\0'; cif_extract_kv_string(start, blk_end, "_struct_conf.end_auth_asym_id", buf, sizeof(buf));
    strip_whitespace_fixed(buf, (int)strlen(buf), conf->end.chain, 2);

    buf[0] = '\0'; cif_extract_kv_string(start, blk_end, "_struct_conf.end_auth_seq_id", buf, sizeof(buf));
    conf->end.resSeq = (uint32_t)atoi(buf);

    buf[0] = '\0'; cif_extract_kv_string(start, blk_end, "_struct_conf.pdbx_end_PDB_ins_code", buf, sizeof(buf));
    conf->end.iCode = (buf[0] == '?' || buf[0] == '\0') ? '\0' : buf[0];

    conf->length = (uint16_t)(conf->end.resSeq - conf->start.resSeq + 1);

    buf[0] = '\0'; cif_extract_kv_string(start, blk_end, "_struct_conf.conf_type_id", buf, sizeof(buf));
    if (buf[0]) conf->type_id = encode_struct_conf_type(buf);

    *out_confs = conf;
    return 1;
}

static size_t parse_struct_conf(const char* data, const char* end, MetaStructConf** out_confs) {
    const char* ptr = data;
    MetaStructConf* confs = NULL;
    size_t count = 0, cap = 0;

    while (ptr < end) {
        const char* loop_start = find_token(ptr, end, "loop_");
        if (!loop_start) break;
        ptr = loop_start + 5;

        const char* headers[MAX_LOOP_HEADERS];
        int header_lens[MAX_LOOP_HEADERS];
        int num_headers = 0;

        while (ptr < end && num_headers < MAX_LOOP_HEADERS) {
            while (ptr < end && (*ptr == ' ' || *ptr == '\r' || *ptr == '\n')) ptr++;
            if (ptr >= end || *ptr != '_') break;
            headers[num_headers] = ptr;
            while (ptr < end && *ptr != ' ' && *ptr != '\n' && *ptr != '\r') ptr++;
            header_lens[num_headers] = (int)(ptr - headers[num_headers]);
            num_headers++;
        }

        if (num_headers == 0 || ptr >= end) continue;

        if (header_lens[0] > 13 && strncmp(headers[0], "_struct_conf.", 13) == 0) {
            StructConfCols cols = map_struct_conf_cols(headers, header_lens, num_headers);
            int col_off[MAX_LOOP_HEADERS] = {0}, col_wid[MAX_LOOP_HEADERS] = {0};
            detect_loop_col_layout(ptr, end, num_headers, col_off, col_wid);

            while (!cif_is_block_end(ptr, end)) {
                if (*ptr == '\n' || *ptr == '\r' || *ptr == ' ') { ptr++; continue; }

                if (count >= cap) {
                    cap = (cap == 0) ? 32 : cap * 2;
                    confs = realloc(confs, cap * sizeof(MetaStructConf));
                }

                MetaStructConf* conf = &confs[count++];
                memset(conf, 0, sizeof(*conf));

                if (cols.beg_comp != -1)  conf->start.resType = res3_to_type_id(ptr + col_off[cols.beg_comp], col_wid[cols.beg_comp]);
                if (cols.beg_chain != -1) strip_whitespace_fixed(ptr + col_off[cols.beg_chain], col_wid[cols.beg_chain], conf->start.chain, 2);
                if (cols.beg_seq != -1)   conf->start.resSeq = (uint32_t)parse_number(ptr + col_off[cols.beg_seq], col_wid[cols.beg_seq]);
                if (cols.beg_icode != -1) { char ic = ptr[col_off[cols.beg_icode]]; conf->start.iCode = (ic == '?') ? '\0' : ic; }

                if (cols.end_comp != -1)  conf->end.resType = res3_to_type_id(ptr + col_off[cols.end_comp], col_wid[cols.end_comp]);
                if (cols.end_chain != -1) strip_whitespace_fixed(ptr + col_off[cols.end_chain], col_wid[cols.end_chain], conf->end.chain, 2);
                if (cols.end_seq != -1)   conf->end.resSeq = (uint32_t)parse_number(ptr + col_off[cols.end_seq], col_wid[cols.end_seq]);
                if (cols.end_icode != -1) { char ic = ptr[col_off[cols.end_icode]]; conf->end.iCode = (ic == '?') ? '\0' : ic; }

                conf->length = (uint16_t)(conf->end.resSeq - conf->start.resSeq + 1);

                if (cols.conf_type != -1) {
                    char type_buf[16] = {0};
                    strip_whitespace_fixed(ptr + col_off[cols.conf_type], col_wid[cols.conf_type], type_buf, 15);
                    conf->type_id = encode_struct_conf_type(type_buf);
                }

                ptr = cif_skip_line(ptr, end);
            }
        } else {
            while (ptr < end && !cif_is_block_end(ptr, end)) ptr++;
        }
    }

    if (count == 0) {
        return parse_struct_conf_single(data, end, out_confs);
    }

    *out_confs = confs;
    return count;
}

static size_t parse_templates(const char* data, const char* end, TemplateEntry** out_temps) {
    const char* key = "_ma_template_ref_db_details.db_accession_code";
    const char* ref_start = find_token(data, end, key);
    if (!ref_start) return 0;

    // --- Case 1: Single-row KV block (count == 1) ---
    if (cif_has_inline_value(ref_start + strlen(key), end)) {
        const char* ref_end = cif_find_block_end(ref_start, end);
        TemplateEntry* t = calloc(1, sizeof(TemplateEntry));
        if (!t) return 0;

        char id_buf[16] = {0};
        t->accession = cif_extract_kv_alloc(ref_start, ref_end, "_ma_template_ref_db_details.db_accession_code");
        t->db_name   = cif_extract_kv_alloc(ref_start, ref_end, "_ma_template_ref_db_details.db_name");
        cif_extract_kv_string(ref_start, ref_end, "_ma_template_ref_db_details.template_id", id_buf, sizeof(id_buf));
        t->template_id = (uint16_t)atoi(id_buf);

        const char* det_key = "_ma_template_details.template_auth_asym_id";
        const char* det_start = find_token(data, end, det_key);
        if (det_start && cif_has_inline_value(det_start + strlen(det_key), end)) {
            const char* det_end = cif_find_block_end(det_start, end);
            t->asym_id = cif_extract_kv_alloc(det_start, det_end, det_key);
        }
        if (!t->asym_id) t->asym_id = str_dup_n("A", 1);

        *out_temps = t;
        return 1;
    }

    // --- Case 2: Multi-row loop_ block ---
    const char* ptr = find_token(ref_start, end, "_ma_template_ref_db_details.template_id");
    if (!ptr) return 0;
    ptr = cif_skip_line(ptr, end);

    size_t cap = 4, count = 0;
    TemplateEntry* temps = malloc(cap * sizeof(TemplateEntry));

    while (!cif_is_block_end(ptr, end)) {
        if (*ptr == ' ' || *ptr == '\n' || *ptr == '\r') { ptr++; continue; }
        if (count >= cap) { cap *= 2; temps = realloc(temps, cap * sizeof(TemplateEntry)); }
        memset(&temps[count], 0, sizeof(TemplateEntry));

        char acc_buf[64] = {0}, db_buf[64] = {0}, id_buf[16] = {0};
        ptr = cif_next_token(ptr, end, acc_buf, sizeof(acc_buf), NULL);
        ptr = cif_next_token(ptr, end, db_buf, sizeof(db_buf), NULL);
        ptr = cif_next_token(ptr, end, id_buf, sizeof(id_buf), NULL);

        temps[count].accession   = str_dup_n(acc_buf, strlen(acc_buf));
        temps[count].db_name     = str_dup_n(db_buf, strlen(db_buf));
        temps[count].template_id = (uint16_t)parse_number(id_buf, (int)strlen(id_buf));
        temps[count].asym_id     = str_dup_n("A", 1);
        count++;
        ptr = cif_skip_line(ptr, end);
    }

    if (count == 0) {
        free(temps);
        return 0;
    }

    const char* det_ptr = find_token(data, end, "_ma_template_details.template_trans_matrix_id");
    if (det_ptr) {
        det_ptr = cif_skip_line(det_ptr, end);
        size_t d_count = 0;
        while (!cif_is_block_end(det_ptr, end)) {
            if (*det_ptr == ' ' || *det_ptr == '\n' || *det_ptr == '\r') { det_ptr++; continue; }
            char auth_buf[64] = {0};
            det_ptr = cif_next_token(det_ptr, end, NULL, 0, NULL); // ordinal_id
            det_ptr = cif_next_token(det_ptr, end, NULL, 0, NULL); // target_asym_id
            det_ptr = cif_next_token(det_ptr, end, auth_buf, sizeof(auth_buf), NULL);
            if (d_count < count && auth_buf[0] != '\0') {
                free(temps[d_count].asym_id);
                temps[d_count].asym_id = str_dup_n(auth_buf, strlen(auth_buf));
            }
            d_count++;
            det_ptr = cif_skip_line(det_ptr, end);
        }
    }

    *out_temps = temps;
    return count;
}

static bool parse_dbref(const char* data, const char* end, char chain_id, AfdbMetadata* meta) {
    const char* blk_start = find_token(data, end, "_ma_target_ref_db_details.");
    if (!blk_start) return false;
    const char* blk_end = cif_find_block_end(blk_start, end);

    char* acc = cif_extract_kv_alloc(blk_start, blk_end, "_ma_target_ref_db_details.db_accession");
    if (!acc) return false;

    DbRefEntry* ref = &meta->dbref;
    memset(ref, 0, sizeof(*ref));

    ref->chain          = chain_id;
    ref->iCode_begin    = '?';
    ref->iCode_end      = '?';
    ref->db_iCode_begin = '?';
    ref->db_iCode_end   = '?';

    ref->db_accession = acc;
    ref->db_id_code   = cif_extract_kv_alloc(blk_start, blk_end, "_ma_target_ref_db_details.db_code");
    ref->db_name      = cif_extract_kv_alloc(blk_start, blk_end, "_ma_target_ref_db_details.db_name");
    ref->version_date = cif_extract_kv_alloc(blk_start, blk_end, "_ma_target_ref_db_details.seq_db_sequence_version_date");
    ref->gene         = cif_extract_kv_alloc(blk_start, blk_end, "_ma_target_ref_db_details.gene_name");
    ref->organism     = cif_extract_kv_alloc(blk_start, blk_end, "_ma_target_ref_db_details.organism_scientific");

    char* iso = cif_extract_kv_alloc(blk_start, blk_end, "_ma_target_ref_db_details.seq_db_isoform");
    if (iso && strcmp(iso, "?") == 0) {
        free(iso);
        iso = NULL;
    }
    ref->isoform = iso;

    char tmp_val[32] = {0};
    cif_extract_kv_string(blk_start, blk_end, "_ma_target_ref_db_details.seq_db_align_begin", tmp_val, sizeof(tmp_val));
    ref->db_seq_begin = atoi(tmp_val);

    tmp_val[0] = '\0';
    cif_extract_kv_string(blk_start, blk_end, "_ma_target_ref_db_details.seq_db_align_end", tmp_val, sizeof(tmp_val));
    ref->db_seq_end = atoi(tmp_val);

    tmp_val[0] = '\0';
    cif_extract_kv_string(blk_start, blk_end, "_ma_target_ref_db_details.seq_db_sequence_checksum", tmp_val, sizeof(tmp_val));
    ref->checksum = parse_hex_uint64(tmp_val);

    tmp_val[0] = '\0';
    cif_extract_kv_string(blk_start, blk_end, "_ma_target_ref_db_details.ncbi_taxonomy_id", tmp_val, sizeof(tmp_val));
    ref->ncbi_tax_id = (uint32_t)atoi(tmp_val);

    return true;
}

static uint32_t parse_audit_revisions(const char* data, const char* end, AuditRevisionEntry** out_revs) {
    const char* ptr = find_token(data, end, "_pdbx_audit_revision_details.type");
    if (!ptr) return 0;
    ptr = cif_skip_line(ptr, end);

    size_t cap = 4;
    uint32_t count = 0;
    AuditRevisionEntry* revs = malloc(cap * sizeof(AuditRevisionEntry));

    while (!cif_is_block_end(ptr, end)) {
        if (*ptr == ' ' || *ptr == '\n' || *ptr == '\r') { ptr++; continue; }
        if (count >= cap) { cap *= 2; revs = realloc(revs, cap * sizeof(AuditRevisionEntry)); }
        memset(&revs[count], 0, sizeof(AuditRevisionEntry));

        char tmp_desc[1024] = {0};
        char tmp_type[256]  = {0};
        bool desc_quoted = false;

        ptr = cif_next_token(ptr, end, NULL, 0, NULL);
        ptr = cif_next_token(ptr, end, tmp_desc, sizeof(tmp_desc), &desc_quoted);
        if (!desc_quoted && strcmp(tmp_desc, "?") == 0) tmp_desc[0] = '\0';

        ptr = cif_next_token(ptr, end, NULL, 0, NULL);
        ptr = cif_next_token(ptr, end, NULL, 0, NULL);
        ptr = cif_next_token(ptr, end, NULL, 0, NULL);
        ptr = cif_next_token(ptr, end, tmp_type, sizeof(tmp_type), NULL);

        size_t dlen = strlen(tmp_desc);
        size_t tlen = strlen(tmp_type);
        revs[count].meta.desc_len = (uint16_t)dlen;
        revs[count].meta.type_len = (uint8_t)tlen;
        if (dlen > 0) revs[count].description = str_dup_n(tmp_desc, dlen);
        if (tlen > 0) revs[count].type        = str_dup_n(tmp_type, tlen);

        count++;
        ptr = cif_skip_line(ptr, end);
    }

    if (count == 0) {
        free(revs);
        return 0;
    }

    ptr = find_token(data, end, "_pdbx_audit_revision_history.revision_date");
    if (ptr) {
        ptr = cif_skip_line(ptr, end);
        uint32_t h_idx = 0;
        while (!cif_is_block_end(ptr, end)) {
            if (*ptr == ' ' || *ptr == '\n' || *ptr == '\r') { ptr++; continue; }

            char maj_s[16] = {0}, min_s[16] = {0}, date_s[16] = {0};
            ptr = cif_next_token(ptr, end, NULL, 0, NULL);
            ptr = cif_next_token(ptr, end, maj_s, sizeof(maj_s), NULL);
            ptr = cif_next_token(ptr, end, min_s, sizeof(min_s), NULL);
            ptr = cif_next_token(ptr, end, NULL, 0, NULL);
            ptr = cif_next_token(ptr, end, date_s, sizeof(date_s), NULL);

            if (h_idx < count) {
                revs[h_idx].meta.major = (uint8_t)atoi(maj_s);
                revs[h_idx].meta.minor = (uint8_t)atoi(min_s);
                revs[h_idx].meta.year  = (uint8_t)(atoi(date_s) - 1973);
                if (strlen(date_s) >= 7)  revs[h_idx].meta.month = (uint8_t)atoi(date_s + 5);
                if (strlen(date_s) >= 10) revs[h_idx].meta.day   = (uint8_t)atoi(date_s + 8);
            }
            h_idx++;
            ptr = cif_skip_line(ptr, end);
        }
    }

    *out_revs = revs;
    return count;
}

bool cif_parse_afdb(const char* cif_data, size_t cif_len,
                    AfdbMetadata* out_meta, ParsedChain* out_chain) {
    memset(out_meta, 0, sizeof(*out_meta));
    memset(out_chain, 0, sizeof(*out_chain));

    const char* end = cif_data + cif_len;
    out_meta->frag_num = 1;
    out_meta->flags = FLAG_AF_MODE | FLAG_AFDB | FLAG_CIF | FLAG_AF2 | FLAG_ENTITY;

    parse_global_ids(cif_data, end, out_meta->accession, &out_meta->frag_num);

    out_meta->global_plddt = parse_global_plddt(cif_data, end);
    if (out_meta->global_plddt >= 0.0f) {
        out_meta->flags |= FLAG_RESOLUTION;
    }

    out_meta->entity_desc = cif_extract_kv_alloc(cif_data, end, "_entity.pdbx_description");

    out_chain->res_count = parse_entity_seq(cif_data, end, &out_chain->residues);
    if (out_chain->res_count <= 0) {
        return false;
    }

    out_meta->conf_count = parse_struct_conf(cif_data, end, &out_meta->confs);
    if (out_meta->conf_count > 0) {
        out_meta->flags |= FLAG_STRUCT;
    }

    out_meta->template_count = parse_templates(cif_data, end, &out_meta->templates);
    if (out_meta->template_count > 0) {
        out_meta->flags |= FLAG_TEMPLATES;
    }

    out_chain->first_seq_id = -1;
    out_chain->chain_id     = 'A';
    parse_atoms(cif_data, end, out_chain);

    if (parse_dbref(cif_data, end, out_chain->chain_id, out_meta)) {
        out_meta->flags |= FLAG_DBREF;
    }

    out_meta->revision_count = parse_audit_revisions(cif_data, end, &out_meta->revisions);
    if (out_meta->revision_count > 0) {
        out_meta->flags |= FLAG_AUDIT;
    }

    return true;
}
