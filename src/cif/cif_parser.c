#include "cif/cif_parser.h"
#include "cif/cif_lexer.h"
#include "cif/afdb_strings.h"
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

typedef struct {
    const char* start;
    const char* end;
} CifSpan;

typedef struct {
    CifSpan header;
    CifSpan entity;
    CifSpan entity_poly;
    CifSpan qa_global;
    CifSpan dbref;
    CifSpan tpl_det;
    CifSpan tpl_ref;
    CifSpan audit_det;
    CifSpan audit_hist;
    CifSpan struct_conf;
    CifSpan atom_site;
} AfdbBlockIndex;

static void index_afdb_blocks(const char* data, const char* end, AfdbBlockIndex* idx) {
    memset(idx, 0, sizeof(*idx));
    const char* p = data;

    while (p < end) {
        while (p < end && (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' || *p == '#')) p++;
        if (p >= end) break;

        const char* blk_start = p;
        const char* blk_end   = cif_find_block_end(blk_start, end);
        const char* tag       = blk_start;

        if ((size_t)(blk_end - tag) >= 5 && memcmp(tag, "loop_", 5) == 0) {
            tag += 5;
            while (tag < blk_end && (*tag == ' ' || *tag == '\t' || *tag == '\r' || *tag == '\n')) tag++;
        }

        size_t rem = (size_t)(blk_end - tag);
        if (rem >= 8 && memcmp(tag, "data_AF-", 8) == 0) {
            idx->header = (CifSpan){ blk_start, blk_end };
        } else if (rem > 1 && tag[0] == '_') {
            switch (tag[1]) {
                case 'a':
                    if (rem >= 11 && memcmp(tag, "_atom_site.", 11) == 0)
                        idx->atom_site = (CifSpan){ blk_start, blk_end };
                    break;
                case 'e':
                    if (rem >= 8 && memcmp(tag, "_entity.", 8) == 0)
                        idx->entity = (CifSpan){ blk_start, blk_end };
                    else if (rem >= 13 && memcmp(tag, "_entity_poly.", 13) == 0)
                        idx->entity_poly = (CifSpan){ blk_start, blk_end };
                    break;
                case 'm':
                    if (rem >= 21 && memcmp(tag, "_ma_qa_metric_global.", 21) == 0)
                        idx->qa_global = (CifSpan){ blk_start, blk_end };
                    else if (rem >= 26 && memcmp(tag, "_ma_target_ref_db_details.", 26) == 0)
                        idx->dbref = (CifSpan){ blk_start, blk_end };
                    else if (rem >= 21 && memcmp(tag, "_ma_template_details.", 21) == 0)
                        idx->tpl_det = (CifSpan){ blk_start, blk_end };
                    else if (rem >= 28 && memcmp(tag, "_ma_template_ref_db_details.", 28) == 0)
                        idx->tpl_ref = (CifSpan){ blk_start, blk_end };
                    break;
                case 'p':
                    if (rem >= 29 && memcmp(tag, "_pdbx_audit_revision_details.", 29) == 0)
                        idx->audit_det = (CifSpan){ blk_start, blk_end };
                    else if (rem >= 29 && memcmp(tag, "_pdbx_audit_revision_history.", 29) == 0)
                        idx->audit_hist = (CifSpan){ blk_start, blk_end };
                    break;
                case 's':
                    if (rem >= 13 && memcmp(tag, "_struct_conf.", 13) == 0)
                        idx->struct_conf = (CifSpan){ blk_start, blk_end };
                    break;
                default:
                    break;
            }
        }
        p = (blk_end < end) ? blk_end + 1 : end;
    }
}

static void parse_global_ids(CifSpan span, char* accession, int* frag_num) {
    if (!span.start) return;
    const char* hdr = find_token(span.start, span.end, "data_AF-");
    if (!hdr) return;

    hdr += 8;
    const char* dash = (const char*)memchr(hdr, '-', (size_t)(span.end - hdr));
    if (!dash) return;

    size_t acc_len = (size_t)(dash - hdr);
    if (acc_len > 12) acc_len = 12;
    memcpy(accession, hdr, acc_len);
    accession[acc_len] = '\0';

    if (dash + 2 < span.end && dash[1] == 'F') {
        *frag_num = atoi(dash + 2);
    }
}

static bool parse_global_plddt(CifSpan span, uint16_t* out_val) {
    if (!span.start) return false;
    const char* key = "_ma_qa_metric_global.metric_value";
    const char* ptr = find_token(span.start, span.end, key);
    if (!ptr) return false;
    ptr = cif_skip_spaces(ptr + strlen(key), span.end);
    if (ptr >= span.end || *ptr == '\n' || *ptr == '\r') return false;

    const char* eol = (const char*)memchr(ptr, '\n', (size_t)(span.end - ptr));
    if (!eol) eol = span.end;
    int32_t v = parse_fixed_decimal2(ptr, (int)(eol - ptr));
    if (v < 0) return false;
    if (v > 65535) v = 65535;
    *out_val = (uint16_t)v;
    return true;
}

static int parse_entity_seq(CifSpan span, ParsedChain* chain) {
    if (!span.start) return 0;
    const char* key = "_entity_poly.pdbx_seq_one_letter_code";
    const char* ptr = find_token(span.start, span.end, key);
    if (!ptr) return 0;

    ptr = cif_skip_spaces(ptr + strlen(key), span.end);

    bool is_semi = false;
    if (ptr < span.end && *ptr == '\r') ptr++;
    if (ptr < span.end && *ptr == '\n') {
        ptr++;
        if (ptr < span.end && *ptr == ';') { is_semi = true; ptr++; }
    } else if (ptr < span.end && *ptr == ';') {
        is_semi = true;
        ptr++;
    }

    int count = 0;
    const char* scan = ptr;
    while (scan < span.end) {
        if (is_semi) {
            if (*scan == ';' && *(scan - 1) == '\n') break;
        } else {
            if (*scan == ' ' || *scan == '\n' || *scan == '\r') break;
        }
        char c = *scan++;
        if (c >= 'A' && c <= 'Z') count++;
    }
    if (count <= 0) return 0;

    if (chain->res_cap < count) {
        ParsedResidue* nr = realloc(chain->residues, (size_t)count * sizeof(ParsedResidue));
        if (!nr) return 0;
        chain->residues = nr;
        chain->res_cap  = count;
    }

    int idx = 0;
    while (ptr < span.end && idx < count) {
        if (is_semi) {
            if (*ptr == ';' && *(ptr - 1) == '\n') break;
        } else {
            if (*ptr == ' ' || *ptr == '\n' || *ptr == '\r') break;
        }
        char c = *ptr++;
        if (c >= 'A' && c <= 'Z') {
            if (c == 'U') c = 'C';
            if (c == 'O') c = 'X';
            ParsedResidue* r = &chain->residues[idx++];
            r->code         = (uint8_t)c;
            r->plddt        = 0;
            r->present_mask = 0;
        }
    }
    return idx;
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

static inline int32_t parse_int_fast(const char* s, int width) {
    const char* p = s;
    const char* e = s + width;
    while (p < e && *p == ' ') p++;
    int neg = 0;
    if (p < e && *p == '-') { neg = 1; p++; }
    int32_t val = 0;
    while (p < e) {
        unsigned d = (unsigned)(*p - '0');
        if (d >= 10U) break;
        val = val * 10 + (int32_t)d;
        p++;
    }
    return neg ? -val : val;
}

static inline int32_t parse_coord3_fast(const char* s, int width, uint8_t* lsb_out, uint8_t bit) {
    const char* p = s;
    const char* e = s + width;
    while (p < e && *p == ' ') p++;
    int neg = 0;
    if (p < e && *p == '-') { neg = 1; p++; }

    int32_t val = 0;
    while (p < e) {
        unsigned d = (unsigned)(*p - '0');
        if (d < 10U) {
            val = val * 10 + (int32_t)d;
            p++;
        } else if (*p == '.') {
            p++;
            if (p + 3 <= e &&
                (unsigned)(p[0] - '0') < 10U &&
                (unsigned)(p[1] - '0') < 10U &&
                (unsigned)(p[2] - '0') < 10U &&
                (p + 3 == e || p[3] == ' ')) {
                val = val * 1000
                    + (int32_t)(p[0] - '0') * 100
                    + (int32_t)(p[1] - '0') * 10
                    + (int32_t)(p[2] - '0');
                break;
            }
            while (p < e) {
                unsigned f = (unsigned)(*p - '0');
                if (f < 10U) val = val * 10 + (int32_t)f;
                else if (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\0') break;
                p++;
            }
            break;
        } else {
            break;
        }
    }
    if (neg) {
        if (val == 0) *lsb_out |= bit;
        return -val;
    }
    return val;
}

static void parse_atoms(CifSpan span, ParsedChain* chain) {
    if (!span.start) return;
    bool layout_ready = false;
    int off[30] = {0}, wid[30] = {0};

    const char* line = span.start;
    const char* end  = span.end;

    while (line < end && !(line + 5 <= end && memcmp(line, "ATOM ", 5) == 0)) {
        line = cif_skip_line(line, end);
    }

    int last_res_idx = -1;
    int cached_aa    = -1;
    int last_res     = chain->res_count - 1;

    while (line < end) {
        const char* eol = (const char*)memchr(line, '\n', (size_t)(end - line));
        if (!eol) eol = end;
        int len = (int)(eol - line);

        if (__builtin_expect(len >= 5 && memcmp(line, "ATOM ", 5) == 0, 1)) {
            if (__builtin_expect(!layout_ready, 0)) {
                layout_ready = detect_atom_columns(line, len, off, wid, chain->fmt_widths);
            }
            if (__builtin_expect(layout_ready, 1)) {
                int seq_id = parse_int_fast(line + off[IDX_RES_SEQ], wid[IDX_RES_SEQ]);
                if (__builtin_expect(chain->first_seq_id == -1, 0)) {
                    chain->first_seq_id = seq_id;
                    char chain_str[3] = {0};
                    strip_whitespace_fixed(line + off[IDX_CHAIN], wid[IDX_CHAIN], chain_str, 2);
                    chain->chain_id = chain_str[0];
                }

                int res_idx = seq_id - 1;
                if (__builtin_expect((unsigned)res_idx < (unsigned)chain->res_count, 1)) {
                    ParsedResidue* res = &chain->residues[res_idx];
                    if (res_idx != last_res_idx) {
                        last_res_idx = res_idx;
                        cached_aa    = dict_aa_index(res->code);
                    }

                    int atom_idx = dict_get_atom_idx_fast(
                        cached_aa, line + off[IDX_ATOM_NAME], wid[IDX_ATOM_NAME], res_idx == last_res
                    );
                    if (__builtin_expect(atom_idx >= 0, 1)) {
                        if (atom_idx == 0 || res->present_mask == 0) {
                            int32_t vb = parse_fixed_decimal2(line + off[IDX_B], wid[IDX_B]);
                            if (vb < 0) vb = 0;
                            if (vb > 65535) vb = 65535;
                            res->plddt = (uint16_t)vb;
                        }

                        uint8_t lsb = 0;
                        res->xyz[atom_idx][0] = parse_coord3_fast(line + off[IDX_X], wid[IDX_X], &lsb, 1U << 0);
                        res->xyz[atom_idx][1] = parse_coord3_fast(line + off[IDX_Y], wid[IDX_Y], &lsb, 1U << 1);
                        res->xyz[atom_idx][2] = parse_coord3_fast(line + off[IDX_Z], wid[IDX_Z], &lsb, 1U << 2);
                        res->xyz_lsb[atom_idx] = lsb;
                        res->present_mask |= (1U << atom_idx);
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
        if      (hl == 16 && memcmp(h, "beg_auth_comp_id", 16) == 0)      c.beg_comp  = i;
        else if (hl == 15 && memcmp(h, "beg_auth_seq_id", 15) == 0)       c.beg_seq   = i;
        else if (hl == 16 && memcmp(h, "beg_auth_asym_id", 16) == 0)      c.beg_chain = i;
        else if (hl == 21 && memcmp(h, "pdbx_beg_PDB_ins_code", 21) == 0) c.beg_icode = i;
        else if (hl == 16 && memcmp(h, "end_auth_comp_id", 16) == 0)      c.end_comp  = i;
        else if (hl == 15 && memcmp(h, "end_auth_seq_id", 15) == 0)       c.end_seq   = i;
        else if (hl == 16 && memcmp(h, "end_auth_asym_id", 16) == 0)      c.end_chain = i;
        else if (hl == 21 && memcmp(h, "pdbx_end_PDB_ins_code", 21) == 0) c.end_icode = i;
        else if (hl == 12 && memcmp(h, "conf_type_id", 12) == 0)          c.conf_type = i;
    }
    return c;
}

static void detect_loop_col_layout(const char* ptr, const char* end, int num_headers, int* col_off, int* col_wid) {
    const char* line_end = (const char*)memchr(ptr, '\n', (size_t)(end - ptr));
    if (!line_end) line_end = end;
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

static bool ensure_conf_cap(AfdbMetadata* meta, size_t need) {
    if (meta->conf_cap >= need) return true;
    size_t ncap = (meta->conf_cap == 0) ? 32U : meta->conf_cap;
    while (ncap < need) ncap *= 2U;
    MetaStructConf* nc = realloc(meta->confs, ncap * sizeof(MetaStructConf));
    if (!nc) return false;
    meta->confs    = nc;
    meta->conf_cap = ncap;
    return true;
}

static size_t parse_struct_conf_single(CifSpan span, AfdbMetadata* meta) {
    const char* key = "_struct_conf.beg_auth_asym_id";
    const char* start = find_token(span.start, span.end, key);
    if (!start || !cif_has_inline_value(start + strlen(key), span.end)) return 0;

    if (!ensure_conf_cap(meta, 1)) return 0;
    MetaStructConf* conf = &meta->confs[0];
    memset(conf, 0, sizeof(*conf));

    char buf[32];
    buf[0] = '\0'; cif_extract_kv_string(span.start, span.end, "_struct_conf.beg_auth_comp_id", buf, sizeof(buf));
    if (buf[0]) conf->start.resType = res3_to_type_id(buf, (int)strlen(buf));

    buf[0] = '\0'; cif_extract_kv_string(span.start, span.end, "_struct_conf.beg_auth_asym_id", buf, sizeof(buf));
    strip_whitespace_fixed(buf, (int)strlen(buf), conf->start.chain, 2);

    buf[0] = '\0'; cif_extract_kv_string(span.start, span.end, "_struct_conf.beg_auth_seq_id", buf, sizeof(buf));
    conf->start.resSeq = (uint32_t)atoi(buf);

    buf[0] = '\0'; cif_extract_kv_string(span.start, span.end, "_struct_conf.pdbx_beg_PDB_ins_code", buf, sizeof(buf));
    conf->start.iCode = (buf[0] == '?' || buf[0] == '\0') ? '\0' : buf[0];

    buf[0] = '\0'; cif_extract_kv_string(span.start, span.end, "_struct_conf.end_auth_comp_id", buf, sizeof(buf));
    if (buf[0]) conf->end.resType = res3_to_type_id(buf, (int)strlen(buf));

    buf[0] = '\0'; cif_extract_kv_string(span.start, span.end, "_struct_conf.end_auth_asym_id", buf, sizeof(buf));
    strip_whitespace_fixed(buf, (int)strlen(buf), conf->end.chain, 2);

    buf[0] = '\0'; cif_extract_kv_string(span.start, span.end, "_struct_conf.end_auth_seq_id", buf, sizeof(buf));
    conf->end.resSeq = (uint32_t)atoi(buf);

    buf[0] = '\0'; cif_extract_kv_string(span.start, span.end, "_struct_conf.pdbx_end_PDB_ins_code", buf, sizeof(buf));
    conf->end.iCode = (buf[0] == '?' || buf[0] == '\0') ? '\0' : buf[0];

    conf->length = (uint16_t)(conf->end.resSeq - conf->start.resSeq + 1);

    buf[0] = '\0'; cif_extract_kv_string(span.start, span.end, "_struct_conf.conf_type_id", buf, sizeof(buf));
    if (buf[0]) conf->type_id = encode_struct_conf_type(buf);

    return 1;
}

static size_t parse_struct_conf(CifSpan span, AfdbMetadata* meta) {
    if (!span.start) return 0;

    const char* ptr = span.start;
    const char* end = span.end;

    if ((size_t)(end - ptr) >= 5 && memcmp(ptr, "loop_", 5) == 0) {
        ptr += 5;
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

        if (num_headers == 0 || ptr >= end) return 0;

        StructConfCols cols = map_struct_conf_cols(headers, header_lens, num_headers);
        int col_off[MAX_LOOP_HEADERS] = {0}, col_wid[MAX_LOOP_HEADERS] = {0};
        detect_loop_col_layout(ptr, end, num_headers, col_off, col_wid);

        size_t count = 0;
        while (!cif_is_block_end(ptr, end)) {
            if (*ptr == '\n' || *ptr == '\r' || *ptr == ' ') { ptr++; continue; }
            if (!ensure_conf_cap(meta, count + 1)) return count;

            MetaStructConf* conf = &meta->confs[count++];
            memset(conf, 0, sizeof(*conf));

            if (cols.beg_comp != -1)  conf->start.resType = res3_to_type_id(ptr + col_off[cols.beg_comp], col_wid[cols.beg_comp]);
            if (cols.beg_chain != -1) strip_whitespace_fixed(ptr + col_off[cols.beg_chain], col_wid[cols.beg_chain], conf->start.chain, 2);
            if (cols.beg_seq != -1)   conf->start.resSeq = (uint32_t)parse_int_fast(ptr + col_off[cols.beg_seq], col_wid[cols.beg_seq]);
            if (cols.beg_icode != -1) { char ic = ptr[col_off[cols.beg_icode]]; conf->start.iCode = (ic == '?') ? '\0' : ic; }

            if (cols.end_comp != -1)  conf->end.resType = res3_to_type_id(ptr + col_off[cols.end_comp], col_wid[cols.end_comp]);
            if (cols.end_chain != -1) strip_whitespace_fixed(ptr + col_off[cols.end_chain], col_wid[cols.end_chain], conf->end.chain, 2);
            if (cols.end_seq != -1)   conf->end.resSeq = (uint32_t)parse_int_fast(ptr + col_off[cols.end_seq], col_wid[cols.end_seq]);
            if (cols.end_icode != -1) { char ic = ptr[col_off[cols.end_icode]]; conf->end.iCode = (ic == '?') ? '\0' : ic; }

            conf->length = (uint16_t)(conf->end.resSeq - conf->start.resSeq + 1);

            if (cols.conf_type != -1) {
                char type_buf[16] = {0};
                strip_whitespace_fixed(ptr + col_off[cols.conf_type], col_wid[cols.conf_type], type_buf, 15);
                conf->type_id = encode_struct_conf_type(type_buf);
            }

            ptr = cif_skip_line(ptr, end);
        }
        return count;
    }

    return parse_struct_conf_single(span, meta);
}

static size_t parse_templates(CifSpan ref_span, CifSpan det_span, TemplateEntry** out_temps) {
    if (!ref_span.start) return 0;

    const char* key = "_ma_template_ref_db_details.db_accession_code";
    const char* ref_start = find_token(ref_span.start, ref_span.end, key);
    if (!ref_start) return 0;

    // --- Case 1: Single-row KV block (count == 1) ---
    if (cif_has_inline_value(ref_start + strlen(key), ref_span.end)) {
        TemplateEntry* t = calloc(1, sizeof(TemplateEntry));
        if (!t) return 0;

        char id_buf[16] = {0};
        t->accession = cif_extract_kv_alloc(ref_span.start, ref_span.end, "_ma_template_ref_db_details.db_accession_code");
        t->db_name   = cif_extract_kv_alloc(ref_span.start, ref_span.end, "_ma_template_ref_db_details.db_name");
        cif_extract_kv_string(ref_span.start, ref_span.end, "_ma_template_ref_db_details.template_id", id_buf, sizeof(id_buf));
        t->template_id = (uint16_t)atoi(id_buf);

        if (det_span.start) {
            const char* det_key = "_ma_template_details.template_auth_asym_id";
            const char* det_start = find_token(det_span.start, det_span.end, det_key);
            if (det_start && cif_has_inline_value(det_start + strlen(det_key), det_span.end)) {
                t->asym_id = cif_extract_kv_alloc(det_span.start, det_span.end, det_key);
            }
        }
        if (!t->asym_id) t->asym_id = str_dup_n("A", 1);

        *out_temps = t;
        return 1;
    }

    // --- Case 2: Multi-row loop_ block ---
    const char* ptr = find_token(ref_start, ref_span.end, "_ma_template_ref_db_details.template_id");
    if (!ptr) return 0;
    ptr = cif_skip_line(ptr, ref_span.end);

    size_t cap = 4, count = 0;
    TemplateEntry* temps = malloc(cap * sizeof(TemplateEntry));
    if (!temps) return 0;

    while (!cif_is_block_end(ptr, ref_span.end)) {
        if (*ptr == ' ' || *ptr == '\n' || *ptr == '\r') { ptr++; continue; }
        if (count >= cap) { cap *= 2; temps = realloc(temps, cap * sizeof(TemplateEntry)); }
        memset(&temps[count], 0, sizeof(TemplateEntry));

        char acc_buf[64] = {0}, db_buf[64] = {0}, id_buf[16] = {0};
        ptr = cif_next_token(ptr, ref_span.end, acc_buf, sizeof(acc_buf), NULL);
        ptr = cif_next_token(ptr, ref_span.end, db_buf, sizeof(db_buf), NULL);
        ptr = cif_next_token(ptr, ref_span.end, id_buf, sizeof(id_buf), NULL);

        temps[count].accession   = str_dup_n(acc_buf, strlen(acc_buf));
        temps[count].db_name     = str_dup_n(db_buf, strlen(db_buf));
        temps[count].template_id = (uint16_t)parse_int_fast(id_buf, (int)strlen(id_buf));
        temps[count].asym_id     = str_dup_n("A", 1);
        count++;
        ptr = cif_skip_line(ptr, ref_span.end);
    }

    if (count == 0) {
        free(temps);
        return 0;
    }

    if (det_span.start) {
        const char* det_ptr = find_token(det_span.start, det_span.end, "_ma_template_details.template_trans_matrix_id");
        if (det_ptr) {
            det_ptr = cif_skip_line(det_ptr, det_span.end);
            size_t d_count = 0;
            while (!cif_is_block_end(det_ptr, det_span.end)) {
                if (*det_ptr == ' ' || *det_ptr == '\n' || *det_ptr == '\r') { det_ptr++; continue; }
                char auth_buf[64] = {0};
                det_ptr = cif_next_token(det_ptr, det_span.end, NULL, 0, NULL);
                det_ptr = cif_next_token(det_ptr, det_span.end, NULL, 0, NULL);
                det_ptr = cif_next_token(det_ptr, det_span.end, auth_buf, sizeof(auth_buf), NULL);
                if (d_count < count && auth_buf[0] != '\0') {
                    free(temps[d_count].asym_id);
                    temps[d_count].asym_id = str_dup_n(auth_buf, strlen(auth_buf));
                }
                d_count++;
                det_ptr = cif_skip_line(det_ptr, det_span.end);
            }
        }
    }

    *out_temps = temps;
    return count;
}

static bool parse_dbref(CifSpan span, char chain_id, AfdbMetadata* meta) {
    if (!span.start) return false;

    char* acc = cif_extract_kv_alloc(span.start, span.end, "_ma_target_ref_db_details.db_accession");
    if (!acc) return false;

    DbRefEntry* ref = &meta->dbref;
    memset(ref, 0, sizeof(*ref));

    ref->chain          = chain_id;
    ref->iCode_begin    = '?';
    ref->iCode_end      = '?';
    ref->db_iCode_begin = '?';
    ref->db_iCode_end   = '?';

    ref->db_accession = acc;
    ref->db_id_code   = cif_extract_kv_alloc(span.start, span.end, "_ma_target_ref_db_details.db_code");
    ref->db_name      = cif_extract_kv_alloc(span.start, span.end, "_ma_target_ref_db_details.db_name");
    ref->version_date = cif_extract_kv_alloc(span.start, span.end, "_ma_target_ref_db_details.seq_db_sequence_version_date");
    ref->gene         = cif_extract_kv_alloc(span.start, span.end, "_ma_target_ref_db_details.gene_name");
    ref->organism     = cif_extract_kv_alloc(span.start, span.end, "_ma_target_ref_db_details.organism_scientific");

    char* iso = cif_extract_kv_alloc(span.start, span.end, "_ma_target_ref_db_details.seq_db_isoform");
    if (iso && strcmp(iso, "?") == 0) {
        free(iso);
        iso = NULL;
    }
    ref->isoform = iso;

    char tmp_val[32] = {0};
    cif_extract_kv_string(span.start, span.end, "_ma_target_ref_db_details.seq_db_align_begin", tmp_val, sizeof(tmp_val));
    ref->db_seq_begin = atoi(tmp_val);

    tmp_val[0] = '\0';
    cif_extract_kv_string(span.start, span.end, "_ma_target_ref_db_details.seq_db_align_end", tmp_val, sizeof(tmp_val));
    ref->db_seq_end = atoi(tmp_val);

    tmp_val[0] = '\0';
    cif_extract_kv_string(span.start, span.end, "_ma_target_ref_db_details.seq_db_sequence_checksum", tmp_val, sizeof(tmp_val));
    ref->checksum = parse_hex_uint64(tmp_val);

    tmp_val[0] = '\0';
    cif_extract_kv_string(span.start, span.end, "_ma_target_ref_db_details.ncbi_taxonomy_id", tmp_val, sizeof(tmp_val));
    ref->ncbi_tax_id = (uint32_t)atoi(tmp_val);

    return true;
}

static int match_audit_preset(const MetaAuditRevision* m, const char* type, const char* desc) {
    for (size_t i = 0; i < NUM_AFDB_AUDIT_PRESETS; i++) {
        const AfdbAuditPreset* p = &AFDB_AUDIT_PRESETS[i];
        if (m->major == p->major && m->minor == p->minor &&
            m->year  == p->year  && m->month == p->month && m->day == p->day &&
            strcmp(type, p->type) == 0 && strcmp(desc, p->description) == 0) {
            return (int)i;
        }
    }
    return -1;
}

static bool parse_audit_revisions(CifSpan det_span, CifSpan hist_span, AfdbMetadata* meta) {
    if (!det_span.start) return false;

    const char* ptr = find_token(det_span.start, det_span.end, "_pdbx_audit_revision_details.type");
    if (!ptr) return false;
    ptr = cif_skip_line(ptr, det_span.end);

    // Fast path: stack-parse up to 6 rows and match AFDB v1..v6 presets (zero heap allocs)
    char st_desc[6][128];
    char st_type[6][32];
    MetaAuditRevision st_meta[6];
    uint32_t st_count = 0;
    bool overflow_stack = false;

    const char* save_det = ptr;
    while (!cif_is_block_end(ptr, det_span.end)) {
        if (*ptr == ' ' || *ptr == '\n' || *ptr == '\r') { ptr++; continue; }
        if (st_count >= 6) { overflow_stack = true; break; }

        char tmp_desc[1024] = {0};
        char tmp_type[256]  = {0};
        bool desc_quoted = false;

        ptr = cif_next_token(ptr, det_span.end, NULL, 0, NULL);
        ptr = cif_next_token(ptr, det_span.end, tmp_desc, sizeof(tmp_desc), &desc_quoted);
        if (!desc_quoted && strcmp(tmp_desc, "?") == 0) tmp_desc[0] = '\0';

        ptr = cif_next_token(ptr, det_span.end, NULL, 0, NULL);
        ptr = cif_next_token(ptr, det_span.end, NULL, 0, NULL);
        ptr = cif_next_token(ptr, det_span.end, NULL, 0, NULL);
        ptr = cif_next_token(ptr, det_span.end, tmp_type, sizeof(tmp_type), NULL);

        size_t dlen = strlen(tmp_desc);
        size_t tlen = strlen(tmp_type);
        if (dlen >= sizeof(st_desc[0]) || tlen >= sizeof(st_type[0])) {
            overflow_stack = true;
            break;
        }
        memcpy(st_desc[st_count], tmp_desc, dlen + 1);
        memcpy(st_type[st_count], tmp_type, tlen + 1);
        memset(&st_meta[st_count], 0, sizeof(MetaAuditRevision));
        st_meta[st_count].desc_len = (uint16_t)dlen;
        st_meta[st_count].type_len = (uint8_t)tlen;
        st_count++;
        ptr = cif_skip_line(ptr, det_span.end);
    }

    if (!overflow_stack && st_count > 0 && hist_span.start) {
        const char* hptr = find_token(hist_span.start, hist_span.end, "_pdbx_audit_revision_history.revision_date");
        if (hptr) {
            hptr = cif_skip_line(hptr, hist_span.end);
            uint32_t h_idx = 0;
            while (!cif_is_block_end(hptr, hist_span.end)) {
                if (*hptr == ' ' || *hptr == '\n' || *hptr == '\r') { hptr++; continue; }
                char maj_s[16] = {0}, min_s[16] = {0}, date_s[16] = {0};
                hptr = cif_next_token(hptr, hist_span.end, NULL, 0, NULL);
                hptr = cif_next_token(hptr, hist_span.end, maj_s, sizeof(maj_s), NULL);
                hptr = cif_next_token(hptr, hist_span.end, min_s, sizeof(min_s), NULL);
                hptr = cif_next_token(hptr, hist_span.end, NULL, 0, NULL);
                hptr = cif_next_token(hptr, hist_span.end, date_s, sizeof(date_s), NULL);

                if (h_idx < st_count) {
                    st_meta[h_idx].major = (uint8_t)atoi(maj_s);
                    st_meta[h_idx].minor = (uint8_t)atoi(min_s);
                    st_meta[h_idx].year  = (uint8_t)(atoi(date_s) - 1973);
                    if (strlen(date_s) >= 7)  st_meta[h_idx].month = (uint8_t)atoi(date_s + 5);
                    if (strlen(date_s) >= 10) st_meta[h_idx].day   = (uint8_t)atoi(date_s + 8);
                }
                h_idx++;
                hptr = cif_skip_line(hptr, hist_span.end);
            }

            if (h_idx == st_count) {
                uint8_t mask = 0;
                int last_p = -1;
                bool all_match = true;
                for (uint32_t i = 0; i < st_count; i++) {
                    int p = match_audit_preset(&st_meta[i], st_type[i], st_desc[i]);
                    if (p < 0 || p <= last_p) {
                        all_match = false;
                        break;
                    }
                    mask |= (uint8_t)(1U << p);
                    last_p = p;
                }
                if (all_match && mask != 0) {
                    meta->audit_mask     = mask;
                    meta->revisions      = NULL;
                    meta->revision_count = 0;
                    return true;
                }
            }
        }
    }

    // Fallback: parse arbitrary revisions into heap array (AUDIT_MASK_RAW = 0xFF)
    ptr = save_det;
    size_t cap = 4;
    uint32_t count = 0;
    AuditRevisionEntry* revs = malloc(cap * sizeof(AuditRevisionEntry));
    if (!revs) return false;

    while (!cif_is_block_end(ptr, det_span.end)) {
        if (*ptr == ' ' || *ptr == '\n' || *ptr == '\r') { ptr++; continue; }
        if (count >= cap) { cap *= 2; revs = realloc(revs, cap * sizeof(AuditRevisionEntry)); }
        memset(&revs[count], 0, sizeof(AuditRevisionEntry));

        char tmp_desc[1024] = {0};
        char tmp_type[256]  = {0};
        bool desc_quoted = false;

        ptr = cif_next_token(ptr, det_span.end, NULL, 0, NULL);
        ptr = cif_next_token(ptr, det_span.end, tmp_desc, sizeof(tmp_desc), &desc_quoted);
        if (!desc_quoted && strcmp(tmp_desc, "?") == 0) tmp_desc[0] = '\0';

        ptr = cif_next_token(ptr, det_span.end, NULL, 0, NULL);
        ptr = cif_next_token(ptr, det_span.end, NULL, 0, NULL);
        ptr = cif_next_token(ptr, det_span.end, NULL, 0, NULL);
        ptr = cif_next_token(ptr, det_span.end, tmp_type, sizeof(tmp_type), NULL);

        size_t dlen = strlen(tmp_desc);
        size_t tlen = strlen(tmp_type);
        revs[count].meta.desc_len = (uint16_t)dlen;
        revs[count].meta.type_len = (uint8_t)tlen;
        if (dlen > 0) revs[count].description = str_dup_n(tmp_desc, dlen);
        if (tlen > 0) revs[count].type        = str_dup_n(tmp_type, tlen);

        count++;
        ptr = cif_skip_line(ptr, det_span.end);
    }

    if (count == 0) {
        free(revs);
        return false;
    }

    if (hist_span.start) {
        const char* hptr = find_token(hist_span.start, hist_span.end, "_pdbx_audit_revision_history.revision_date");
        if (hptr) {
            hptr = cif_skip_line(hptr, hist_span.end);
            uint32_t h_idx = 0;
            while (!cif_is_block_end(hptr, hist_span.end)) {
                if (*hptr == ' ' || *hptr == '\n' || *hptr == '\r') { hptr++; continue; }

                char maj_s[16] = {0}, min_s[16] = {0}, date_s[16] = {0};
                hptr = cif_next_token(hptr, hist_span.end, NULL, 0, NULL);
                hptr = cif_next_token(hptr, hist_span.end, maj_s, sizeof(maj_s), NULL);
                hptr = cif_next_token(hptr, hist_span.end, min_s, sizeof(min_s), NULL);
                hptr = cif_next_token(hptr, hist_span.end, NULL, 0, NULL);
                hptr = cif_next_token(hptr, hist_span.end, date_s, sizeof(date_s), NULL);

                if (h_idx < count) {
                    revs[h_idx].meta.major = (uint8_t)atoi(maj_s);
                    revs[h_idx].meta.minor = (uint8_t)atoi(min_s);
                    revs[h_idx].meta.year  = (uint8_t)(atoi(date_s) - 1973);
                    if (strlen(date_s) >= 7)  revs[h_idx].meta.month = (uint8_t)atoi(date_s + 5);
                    if (strlen(date_s) >= 10) revs[h_idx].meta.day   = (uint8_t)atoi(date_s + 8);
                }
                h_idx++;
                hptr = cif_skip_line(hptr, hist_span.end);
            }
        }
    }

    meta->audit_mask     = AUDIT_MASK_RAW;
    meta->revisions      = revs;
    meta->revision_count = count;
    return true;
}

bool cif_parse_afdb(const char* cif_data, size_t cif_len,
                    AfdbMetadata* out_meta, ParsedChain* out_chain) {
    MetaStructConf* keep_confs    = out_meta->confs;
    size_t          keep_conf_cap = out_meta->conf_cap;
    uint32_t*       keep_nz       = out_meta->neg_zeros;
    uint32_t        keep_nz_cap   = out_meta->neg_zero_cap;
    memset(out_meta, 0, sizeof(*out_meta));
    out_meta->confs        = keep_confs;
    out_meta->conf_cap     = keep_conf_cap;
    out_meta->neg_zeros    = keep_nz;
    out_meta->neg_zero_cap = keep_nz_cap;

    ParsedResidue* keep_res     = out_chain->residues;
    int            keep_res_cap = out_chain->res_cap;
    memset(out_chain, 0, sizeof(*out_chain));
    out_chain->residues = keep_res;
    out_chain->res_cap  = keep_res_cap;

    const char* end = cif_data + cif_len;
    AfdbBlockIndex idx;
    index_afdb_blocks(cif_data, end, &idx);

    out_meta->frag_num = 1;
    out_meta->flags = FLAG_AF_MODE | FLAG_AFDB | FLAG_CIF | FLAG_AF2;

    parse_global_ids(idx.header, out_meta->accession, &out_meta->frag_num);
    if (out_meta->frag_num > 1) {
        out_meta->flags |= FLAG_FRAGMENTED;
    }

    if (parse_global_plddt(idx.qa_global, &out_meta->global_plddt)) {
        out_meta->flags |= FLAG_RESOLUTION;
    }

    if (idx.entity.start) {
        out_meta->entity_desc = cif_extract_kv_alloc(idx.entity.start, idx.entity.end, "_entity.pdbx_description");
        if (out_meta->entity_desc && out_meta->entity_desc[0] != '\0') {
            out_meta->flags |= FLAG_ENTITY;
        }
    }

    out_chain->res_count = parse_entity_seq(idx.entity_poly, out_chain);
    if (out_chain->res_count <= 0) {
        return false;
    }

    out_meta->conf_count = parse_struct_conf(idx.struct_conf, out_meta);
    if (out_meta->conf_count > 0) {
        out_meta->flags |= FLAG_STRUCT;
    }

    out_meta->template_count = parse_templates(idx.tpl_ref, idx.tpl_det, &out_meta->templates);
    if (out_meta->template_count > 0) {
        out_meta->flags |= FLAG_TEMPLATES;
    }

    out_chain->first_seq_id = -1;
    out_chain->chain_id     = 'A';
    parse_atoms(idx.atom_site, out_chain);

    if (parse_dbref(idx.dbref, out_chain->chain_id, out_meta)) {
        out_meta->flags |= FLAG_DBREF;
    }

    if (parse_audit_revisions(idx.audit_det, idx.audit_hist, out_meta)) {
        out_meta->flags |= FLAG_AUDIT;
    }

    return true;
}
