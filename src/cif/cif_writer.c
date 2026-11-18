#include "cif/cif_writer.h"
#include "cif/afdb_strings.h"
#include "binary/coord_codec.h"
#include "chem/dict.h"
#include "chem/struct_conf.h"
#include "util/str_util.h"
#include <stdbool.h>
#include <string.h>

static void write_wrapped_seq(FILE* out_f, const uint8_t* seq, uint32_t len) {
    if (len <= 80) {
        fprintf(out_f, "%.*s\n", (int)len, seq);
        return;
    }
    fprintf(out_f, "\n;");
    for (uint32_t i = 0; i < len; i++) {
        fputc(seq[i], out_f);
        if ((i + 1) % 80 == 0 && i != len - 1) fputc('\n', out_f);
    }
    fprintf(out_f, "\n;\n");
}

static void get_res3_clean(uint8_t one_letter, char out_res3[4]) {
    const ResidueDict* rd = dict_find_aa(one_letter);
    if (rd) {
        strip_whitespace_fixed(rd->name, 4, out_res3, 3);
    } else {
        memcpy(out_res3, "UNK", 4);
    }
}

static void write_chem_and_citation_blocks(FILE* out_f, const AfdbMetadata* meta, const ChainView* chain) {
    bool present[26] = {false};
    for (uint32_t i = 0; i < chain->total_res; i++) {
        uint8_t c = chain->seq[i];
        if (c >= 'A' && c <= 'Z') present[c - 'A'] = true;
    }

    fprintf(out_f, AFDB_ENTRY_ID_FORMAT, meta->accession, meta->frag_num, meta->accession, meta->frag_num);

    fprintf(out_f, "%s", AFDB_ATOM_TYPE_HEADER);
    if (present['C' - 'A'] || present['M' - 'A']) {
        fprintf(out_f, "S \n");
    }
    fprintf(out_f, "#\n");

    fprintf(out_f, "%s", AFDB_AUDIT_BLOCK);

    fprintf(out_f, "%s", AFDB_CHEM_COMP_HEADER);
    for (size_t k = 0; k < NUM_CHEM_COMP_ROWS; k++) {
        if (present[CHEM_COMP_ROWS[k].code - 'A']) {
            fprintf(out_f, "%s\n", CHEM_COMP_ROWS[k].str);
        }
    }
    fprintf(out_f, "#\n");

    fprintf(out_f, "%s", AFDB_CITATION_BLOCK);
    fprintf(out_f, AFDB_DATABASE_2_FORMAT, meta->accession, meta->frag_num);
}

static void write_entity_blocks(FILE* out_f, const AfdbMetadata* meta, const ChainView* chain, int num_width) {
    fprintf(out_f, "%s", AFDB_ENTITY_BLOCK_1);
    fprintf(out_f, "_entity.pdbx_description         \"%s\"\n", meta->entity_desc);
    fprintf(out_f, "%s", AFDB_ENTITY_BLOCK_2);

    fprintf(out_f, "_entity_poly.pdbx_seq_one_letter_code     ");
    write_wrapped_seq(out_f, chain->seq, chain->total_res);

    fprintf(out_f, "_entity_poly.pdbx_seq_one_letter_code_can ");
    write_wrapped_seq(out_f, chain->seq, chain->total_res);

    fprintf(out_f, "_entity_poly.pdbx_strand_id               %c\n", chain->chain_id);
    fprintf(out_f, "_entity_poly.type                         polypeptide(L)\n#\n");

    fprintf(out_f, "%s", AFDB_ENTITY_POLY_SEQ_HEADER);
    for (uint32_t i = 0; i < chain->total_res; i++) {
        char res3[4];
        get_res3_clean(chain->seq[i], res3);
        fprintf(out_f, "1 n %-3s %-*u \n", res3, num_width, i + 1);
    }
    fprintf(out_f, "#\n");
}

static void write_qa_metric_blocks(FILE* out_f, const AfdbMetadata* meta, const ChainView* chain, int num_width) {
    fprintf(out_f, "%s", (meta->template_count > 0) ? AFDB_MA_DATA_MULTI : AFDB_MA_DATA_SINGLE);
    fprintf(out_f, "%s", AFDB_MA_MODEL_LIST);

    if (meta->flags & FLAG_RESOLUTION) {
        fprintf(out_f, AFDB_QA_METRIC_GLOBAL_FORMAT, meta->global_plddt);
    }

    fprintf(out_f, "%s", AFDB_QA_METRIC_LOCAL_HEADER);

    CoordStream cs;
    coord_stream_init(&cs, chain);

    for (uint32_t i = 0; i < chain->total_res; i++) {
        const ResidueDict* rd = dict_find_aa(chain->seq[i]);
        if (!rd) continue;

        char res3[4];
        strip_whitespace_fixed(rd->name, 4, res3, 3);

        int num_atoms = (i == chain->total_res - 1) ? rd->numHeavyAtomsTerm : rd->numHeavyAtoms;
        int32_t b = 0;

        for (int a = 0; a < num_atoms; a++) {
            int32_t tb = 0, tx = 0, ty = 0, tz = 0;
            coord_stream_step(&cs, &tb, &tx, &ty, &tz);
            if (a == 0) b = tb;
        }

        char sb[16];
        format_fixed_decimal(sb, b, 2);
        fprintf(out_f, "%c %-3s %-*d 2 %-5s 1 %-*u \n",
                chain->chain_id, res3, num_width, chain->first_seq_id + (int)i, sb, num_width, i + 1);
    }
    fprintf(out_f, "#\n");
    fprintf(out_f, "%s", AFDB_MA_POST_BLOCK);
}

static void write_dbref_block(FILE* out_f, const AfdbMetadata* meta) {
    if (!(meta->flags & FLAG_DBREF)) return;

    const MetaDBRef*     ref = &meta->dbref;
    const AFDBTargetExt* ext = &meta->target_ext;

    fprintf(out_f, "_ma_target_ref_db_details.db_accession                 %.12s\n", ref->db_accession);
    fprintf(out_f, "_ma_target_ref_db_details.db_code                      %.12s\n", ref->db_id_code);
    fprintf(out_f, "_ma_target_ref_db_details.db_name                      %.8s\n", ref->db_name);
    if (meta->gene)     fprintf(out_f, "_ma_target_ref_db_details.gene_name                    %s\n", meta->gene);
    fprintf(out_f, "_ma_target_ref_db_details.ncbi_taxonomy_id             %u\n", ext->ncbi_tax_id);
    if (meta->organism) fprintf(out_f, "_ma_target_ref_db_details.organism_scientific          \"%s\"\n", meta->organism);
    fprintf(out_f, "_ma_target_ref_db_details.seq_db_align_begin           %d\n", ref->db_seq_begin);
    fprintf(out_f, "_ma_target_ref_db_details.seq_db_align_end             %d\n", ref->db_seq_end);
    fprintf(out_f, "_ma_target_ref_db_details.seq_db_isoform               %s\n", meta->isoform ? meta->isoform : "?");
    fprintf(out_f, "_ma_target_ref_db_details.seq_db_sequence_checksum     %016llX\n", (unsigned long long)ext->checksum);
    fprintf(out_f, "_ma_target_ref_db_details.seq_db_sequence_version_date %.*s\n", 10, ext->version_date);
    fprintf(out_f, "_ma_target_ref_db_details.target_entity_id             1\n#\n");
}

static void write_template_blocks(FILE* out_f, const AfdbMetadata* meta, char chain_id) {
    if (meta->template_count == 0) return;

    fprintf(out_f, "%s", AFDB_MA_TEMPLATE_DETAILS_HEADER);
    for (size_t i = 0; i < meta->template_count; i++) {
        char auth = meta->templates[i].template_asym_id ? meta->templates[i].template_asym_id : 'A';
        fprintf(out_f, "%u %c %c 2 polymer %u 1 \"reference database\" 1 \n",
                (uint32_t)(i + 1), chain_id, auth, meta->templates[i].template_id);
    }
    fprintf(out_f, "#\n");

    fprintf(out_f, "%s", AFDB_MA_TEMPLATE_REF_DB_DETAILS_HEADER);
    for (size_t i = 0; i < meta->template_count; i++) {
        fprintf(out_f, "%s %s %u \n",
                meta->templates[i].accession, meta->templates[i].db_name, meta->templates[i].template_id);
    }
    fprintf(out_f, "#\n");
    fprintf(out_f, "%s", AFDB_MA_TEMPLATE_TRANS_MATRIX);
}

static void format_cif_quoted(const char* src, size_t len, char* dst, size_t dst_size) {
    if (len == 0) {
        snprintf(dst, dst_size, "?");
    } else if (strchr(src, ' ')) {
        snprintf(dst, dst_size, "\"%s\"", src);
    } else {
        snprintf(dst, dst_size, "%s", src);
    }
}

static void write_audit_revision_blocks(FILE* out_f, const AfdbMetadata* meta) {
    if (meta->revision_count == 0) return;

    int w_desc = 1, w_type = 1, w_maj = 1, w_min = 1;
    for (uint32_t i = 0; i < meta->revision_count; i++) {
        const AuditRevisionEntry* r = &meta->revisions[i];
        int l_desc = (r->meta.desc_len == 0) ? 1 : (int)r->meta.desc_len + (strchr(r->description, ' ') ? 2 : 0);
        int l_type = (r->meta.type_len == 0) ? 1 : (int)r->meta.type_len + (strchr(r->type, ' ') ? 2 : 0);
        int l_maj  = snprintf(NULL, 0, "%d", r->meta.major);
        int l_min  = snprintf(NULL, 0, "%d", r->meta.minor);

        if (l_desc > w_desc) w_desc = l_desc;
        if (l_type > w_type) w_type = l_type;
        if (l_maj  > w_maj)  w_maj  = l_maj;
        if (l_min  > w_min)  w_min  = l_min;
    }

    fprintf(out_f, "%s", AFDB_AUDIT_REVISION_DETAILS_HEADER);
    for (uint32_t i = 0; i < meta->revision_count; i++) {
        const AuditRevisionEntry* r = &meta->revisions[i];
        char q_desc[1028], q_type[260];
        format_cif_quoted(r->description, r->meta.desc_len, q_desc, sizeof(q_desc));
        format_cif_quoted(r->type, r->meta.type_len, q_type, sizeof(q_type));

        fprintf(out_f, "\"Structure model\" %-*s %-*d repository %-*d %-*s \n",
                w_desc, q_desc, w_maj, r->meta.major, w_maj, r->meta.major, w_type, q_type);
    }
    fprintf(out_f, "#\n");

    fprintf(out_f, "%s", AFDB_AUDIT_REVISION_HISTORY_HEADER);
    for (uint32_t i = 0; i < meta->revision_count; i++) {
        const MetaAuditRevision* m = &meta->revisions[i].meta;
        fprintf(out_f, "\"Structure model\" %-*d %-*d %-*d %04d-%02d-%02d \n",
                w_maj, m->major, w_min, m->minor, w_maj, m->major,
                m->year + 1973, m->month, m->day);
    }
    fprintf(out_f, "#\n");
}

static void write_poly_seq_scheme_block(FILE* out_f, const ChainView* chain, int num_width) {
    fprintf(out_f, "%s", AFDB_POLY_SEQ_SCHEME_HEADER);
    for (uint32_t i = 0; i < chain->total_res; i++) {
        char res3[4];
        get_res3_clean(chain->seq[i], res3);
        int seq_id = chain->first_seq_id + (int)i;
        fprintf(out_f, "%c %-*d 1 n %-3s . %-3s %-*d %c %-*u \n",
                chain->chain_id, num_width, seq_id, res3, res3,
                num_width, seq_id, chain->chain_id, num_width, i + 1);
    }
    fprintf(out_f, "#\n");
}

static void write_struct_conf_single(FILE* out_f, const MetaStructConf* c) {
    const char* type_str = decode_struct_conf_type(c->type_id);
    const char* r3_s = type_id_to_res3(c->start.resType);
    const char* r3_e = type_id_to_res3(c->end.resType);

    fprintf(out_f, "_struct_conf.beg_auth_asym_id      %s\n", c->start.chain);
    fprintf(out_f, "_struct_conf.beg_auth_comp_id      %s\n", r3_s);
    fprintf(out_f, "_struct_conf.beg_auth_seq_id       %u\n", c->start.resSeq);
    fprintf(out_f, "_struct_conf.beg_label_asym_id     %s\n", c->start.chain);
    fprintf(out_f, "_struct_conf.beg_label_comp_id     %s\n", r3_s);
    fprintf(out_f, "_struct_conf.beg_label_seq_id      %u\n", c->start.resSeq);
    fprintf(out_f, "_struct_conf.conf_type_id          %s\n", type_str);
    fprintf(out_f, "_struct_conf.end_auth_asym_id      %s\n", c->end.chain);
    fprintf(out_f, "_struct_conf.end_auth_comp_id      %s\n", r3_e);
    fprintf(out_f, "_struct_conf.end_auth_seq_id       %u\n", c->end.resSeq);
    fprintf(out_f, "_struct_conf.end_label_asym_id     %s\n", c->end.chain);
    fprintf(out_f, "_struct_conf.end_label_comp_id     %s\n", r3_e);
    fprintf(out_f, "_struct_conf.end_label_seq_id      %u\n", c->end.resSeq);
    fprintf(out_f, "_struct_conf.id                    %s1\n", type_str);
    fprintf(out_f, "_struct_conf.pdbx_beg_PDB_ins_code %c\n", c->start.iCode ? c->start.iCode : '?');
    fprintf(out_f, "_struct_conf.pdbx_end_PDB_ins_code %c\n", c->end.iCode ? c->end.iCode : '?');
}

static void write_struct_conf_loop(FILE* out_f, const MetaStructConf* confs, size_t count) {
    fprintf(out_f, "%s", AFDB_STRUCT_CONF_HEADER);

    int counts[256] = {0};
    int w_bc = 1, w_br = 1, w_bs = 1, w_ec = 1, w_er = 1, w_es = 1, w_type = 1, w_id = 1;

    for (size_t i = 0; i < count; i++) {
        const char* type_str = decode_struct_conf_type(confs[i].type_id);
        counts[confs[i].type_id]++;

        int l;
        l = (int)strlen(confs[i].start.chain);                  if (l > w_bc)   w_bc   = l;
        l = (int)strlen(type_id_to_res3(confs[i].start.resType)); if (l > w_br) w_br   = l;
        l = snprintf(NULL, 0, "%u", confs[i].start.resSeq);     if (l > w_bs)   w_bs   = l;
        l = (int)strlen(confs[i].end.chain);                    if (l > w_ec)   w_ec   = l;
        l = (int)strlen(type_id_to_res3(confs[i].end.resType)); if (l > w_er)   w_er   = l;
        l = snprintf(NULL, 0, "%u", confs[i].end.resSeq);       if (l > w_es)   w_es   = l;
        l = (int)strlen(type_str);                              if (l > w_type) w_type = l;
        l = snprintf(NULL, 0, "%s%d", type_str, counts[confs[i].type_id]); if (l > w_id) w_id = l;
    }

    memset(counts, 0, sizeof(counts));

    for (size_t i = 0; i < count; i++) {
        const char* type_str = decode_struct_conf_type(confs[i].type_id);
        char id_str[32];
        snprintf(id_str, sizeof(id_str), "%s%d", type_str, ++counts[confs[i].type_id]);

        char i_s[2] = { confs[i].start.iCode ? confs[i].start.iCode : '?', 0 };
        char i_e[2] = { confs[i].end.iCode   ? confs[i].end.iCode   : '?', 0 };
        const char* r3_s = type_id_to_res3(confs[i].start.resType);
        const char* r3_e = type_id_to_res3(confs[i].end.resType);

        fprintf(out_f, "%-*s %-*s %-*u %-*s %-*s %-*u %-*s %-*s %-*s %-*u %-*s %-*s %-*u %-*s %s %s \n",
                w_bc, confs[i].start.chain, w_br, r3_s, w_bs, confs[i].start.resSeq,
                w_bc, confs[i].start.chain, w_br, r3_s, w_bs, confs[i].start.resSeq,
                w_type, type_str,
                w_ec, confs[i].end.chain, w_er, r3_e, w_es, confs[i].end.resSeq,
                w_ec, confs[i].end.chain, w_er, r3_e, w_es, confs[i].end.resSeq,
                w_id, id_str, i_s, i_e);
    }
}

static void write_struct_conf_blocks(FILE* out_f, const MetaStructConf* confs, size_t count) {
    if (count == 0 || !confs) return;

    if (count == 1) write_struct_conf_single(out_f, &confs[0]);
    else            write_struct_conf_loop(out_f, confs, count);
    fprintf(out_f, "#\n");

    uint8_t unique_types[256];
    int seen[256] = {0};
    size_t unique_count = 0;
    int w_ct = 1;

    for (size_t i = 0; i < count; i++) {
        uint8_t tid = confs[i].type_id;
        if (!seen[tid]) {
            seen[tid] = 1;
            unique_types[unique_count++] = tid;
            int l = (int)strlen(decode_struct_conf_type(tid));
            if (l > w_ct) w_ct = l;
        }
    }

    if (unique_count == 1) {
        fprintf(out_f, "_struct_conf_type.criteria DSSP\n");
        fprintf(out_f, "_struct_conf_type.id       %s\n", decode_struct_conf_type(unique_types[0]));
    } else {
        fprintf(out_f, "%s", AFDB_STRUCT_CONF_TYPE_HEADER);
        for (size_t i = 0; i < unique_count; i++) {
            fprintf(out_f, "DSSP %-*s \n", w_ct, decode_struct_conf_type(unique_types[i]));
        }
    }
    fprintf(out_f, "#\n");
}

static void write_struct_ref_blocks(FILE* out_f, const AfdbMetadata* meta, const ChainView* chain) {
    if (!(meta->flags & FLAG_DBREF)) return;

    const MetaDBRef* ref = &meta->dbref;
    fprintf(out_f, AFDB_STRUCT_REF_FORMAT,
            ref->db_id_code, ref->db_name,
            ref->db_seq_begin, ref->db_seq_end,
            ref->db_accession, meta->isoform ? meta->isoform : "?");

    write_wrapped_seq(out_f, chain->seq, chain->total_res);
    fprintf(out_f, "#\n");

    int last_seq_id = chain->first_seq_id + (int)chain->total_res - 1;
    fprintf(out_f, AFDB_STRUCT_REF_SEQ_FORMAT,
            ref->db_seq_begin, ref->db_seq_end,
            meta->accession, meta->frag_num,
            chain->first_seq_id, last_seq_id,
            ref->db_accession, chain->chain_id,
            chain->first_seq_id, last_seq_id);
}

void cif_write_afdb_metadata(FILE* out_f, const AfdbMetadata* meta, const ChainView* chain) {
    int num_width = snprintf(NULL, 0, "%u", chain->total_res);

    write_chem_and_citation_blocks(out_f, meta, chain);
    write_entity_blocks(out_f, meta, chain, num_width);
    write_qa_metric_blocks(out_f, meta, chain, num_width);
    write_dbref_block(out_f, meta);
    write_template_blocks(out_f, meta, chain->chain_id);
    write_audit_revision_blocks(out_f, meta);

    fprintf(out_f, "%s", AFDB_PDBX_DATA_USAGE);

    char last_date[16] = "2021-07-01";
    if (meta->revision_count > 0) {
        const MetaAuditRevision* m = &meta->revisions[meta->revision_count - 1].meta;
        snprintf(last_date, sizeof(last_date), "%04d-%02d-%02d", m->year + 1973, m->month, m->day);
    }
    fprintf(out_f, AFDB_DATABASE_STATUS_FORMAT, meta->accession, meta->frag_num, last_date);

    write_poly_seq_scheme_block(out_f, chain, num_width);
    fprintf(out_f, "%s", AFDB_SOFTWARE_ASYM);
    fprintf(out_f, AFDB_STRUCT_ASYM_FORMAT, chain->chain_id);

    write_struct_conf_blocks(out_f, meta->confs, meta->conf_count);
    write_struct_ref_blocks(out_f, meta, chain);
}

void cif_write_atom_site(FILE* out_f, const ChainView* chain, const char* accession) {
    const uint8_t* w = chain->fmt_widths;
    fprintf(out_f, "%s", AFDB_ATOM_SITE_HEADER);

    int atom_id = 1;
    CoordStream cs;
    coord_stream_init(&cs, chain);

    for (uint32_t i = 0; i < chain->total_res; i++) {
        uint8_t code = chain->seq[i];
        const ResidueDict* rd = dict_find_aa(code);
        if (!rd) continue;

        int num_atoms = (i == chain->total_res - 1) ? rd->numHeavyAtomsTerm : rd->numHeavyAtoms;
        int seq_id = chain->first_seq_id + (int)i;

        int32_t dec_b[30] = {0}, dec_x[30] = {0}, dec_y[30] = {0}, dec_z[30] = {0};
        for (int a = 0; a < num_atoms; a++) {
            coord_stream_step(&cs, &dec_b[a], &dec_x[a], &dec_y[a], &dec_z[a]);
        }

        for (int j = 0; j < 37; j++) {
            for (int a = 0; a < num_atoms; a++) {
                char clean_atom[5];
                strip_whitespace_fixed(rd->atomNames[a], 4, clean_atom, 4);
                if (strcmp(clean_atom, ATOM37_ORDER[j]) != 0) continue;

                char elem[3];
                strip_whitespace_fixed(rd->elements[a], 2, elem, 2);

                char sx[16], sy[16], sz[16], sb[16];
                format_fixed_decimal(sx, dec_x[a], 3);
                format_fixed_decimal(sy, dec_y[a], 3);
                format_fixed_decimal(sz, dec_z[a], 3);
                format_fixed_decimal(sb, dec_b[a], 2);

                fprintf(out_f,
                        "ATOM %-*d%-2s%-*s. %-4s%c 1 %-*d? %-*s%-*s%-*s1.00 %-*s? %-*d%-4s%c %-*s1 %s UNP %-*d%c \n",
                        w[0], atom_id, elem, w[1], clean_atom, rd->name, chain->chain_id,
                        w[2], seq_id, w[3], sx, w[4], sy, w[5], sz, w[6], sb,
                        w[2], seq_id, rd->name, chain->chain_id, w[1], clean_atom,
                        accession, w[7], seq_id, code);
                atom_id++;
                break;
            }
        }
    }
}
