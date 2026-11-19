#include "binary/meta_codec.h"
#include "util/str_util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char UNCHAR_PROTEIN[] = "Uncharacterized protein";

static void afdb_meta_free_dynamic_fields(AfdbMetadata* meta) {
    free(meta->entity_desc);
    if (meta->templates) {
        for (size_t i = 0; i < meta->template_count; i++) {
            free(meta->templates[i].accession);
            free(meta->templates[i].db_name);
            free(meta->templates[i].asym_id);
        }
        free(meta->templates);
    }
    free(meta->dbref.db_name);
    free(meta->dbref.db_accession);
    free(meta->dbref.db_id_code);
    free(meta->dbref.version_date);
    free(meta->dbref.gene);
    free(meta->dbref.organism);
    free(meta->dbref.isoform);
    if (meta->revisions) {
        for (uint32_t i = 0; i < meta->revision_count; i++) {
            free(meta->revisions[i].description);
            free(meta->revisions[i].type);
        }
        free(meta->revisions);
    }
}

void afdb_meta_free(AfdbMetadata* meta) {
    if (!meta) return;
    afdb_meta_free_dynamic_fields(meta);
    free(meta->confs);
    free(meta->neg_zeros);
    memset(meta, 0, sizeof(*meta));
}

void afdb_meta_reset_reuse(AfdbMetadata* meta) {
    if (!meta) return;
    afdb_meta_free_dynamic_fields(meta);
    MetaStructConf* keep_confs    = meta->confs;
    size_t          keep_conf_cap = meta->conf_cap;
    uint32_t*       keep_nz       = meta->neg_zeros;
    uint32_t        keep_nz_cap   = meta->neg_zero_cap;
    memset(meta, 0, sizeof(*meta));
    meta->confs        = keep_confs;
    meta->conf_cap     = keep_conf_cap;
    meta->neg_zeros    = keep_nz;
    meta->neg_zero_cap = keep_nz_cap;
}

static size_t tlv_varstr_size(const char* s) {
    if (!s || s[0] == '\0') return 0;
    size_t len = strlen(s);
    if (len > 65535) len = 65535;
    return 1 + bw_varlen_size(len) + len;
}

static void write_tlv_varstr(ByteWriter* w, uint8_t tag, const char* s) {
    if (!s || s[0] == '\0') return;
    size_t len = strlen(s);
    if (len > 65535) len = 65535;
    bw_write_u8(w, tag);
    bw_write_varlen(w, len);
    bw_write(w, s, len);
}

static bool pack_iso_date(const char* s, uint8_t out[3]) {
    if (!s || strlen(s) != 10 || s[4] != '-' || s[7] != '-') return false;
    for (int i = 0; i < 10; i++) {
        if (i == 4 || i == 7) continue;
        if (s[i] < '0' || s[i] > '9') return false;
    }
    int yr = (s[0] - '0') * 1000 + (s[1] - '0') * 100 + (s[2] - '0') * 10 + (s[3] - '0');
    int mo = (s[5] - '0') * 10 + (s[6] - '0');
    int dy = (s[8] - '0') * 10 + (s[9] - '0');
    if (yr < 1973 || yr > 2228 || mo < 1 || mo > 12 || dy < 1 || dy > 31) return false;
    out[0] = (uint8_t)(yr - 1973);
    out[1] = (uint8_t)mo;
    out[2] = (uint8_t)dy;
    return true;
}

static bool has_nondefault_icodes(const DbRefEntry* ref) {
    return (ref->iCode_begin && ref->iCode_begin != '?') ||
           (ref->iCode_end && ref->iCode_end != '?') ||
           (ref->db_iCode_begin && ref->db_iCode_begin != '?') ||
           (ref->db_iCode_end && ref->db_iCode_end != '?');
}

static const char* dbref_id_suffix(const char* accession, const char* id_code) {
    if (!accession || !accession[0] || !id_code) return NULL;
    size_t alen = strlen(accession);
    if (strncmp(id_code, accession, alen) == 0 && id_code[alen] == '_' && id_code[alen + 1] != '\0') {
        return id_code + alen + 1;
    }
    return NULL;
}

static bool dbref_align_matches_chain(const DbRefEntry* ref, const ParsedChain* chain) {
    return ref->db_seq_begin == chain->first_seq_id &&
           ref->db_seq_end   == chain->first_seq_id + chain->res_count - 1;
}

static size_t dbref_tlv_payload_bytes(const AfdbMetadata* meta, const ParsedChain* chain) {
    const DbRefEntry* ref = &meta->dbref;
    size_t bytes = 1; // DBREF_TAG_END

    if (ref->db_name && ref->db_name[0] && strcmp(ref->db_name, "UNP") != 0) {
        bytes += tlv_varstr_size(ref->db_name);
    }
    if (ref->db_accession && ref->db_accession[0] && strcmp(ref->db_accession, meta->accession) != 0) {
        bytes += tlv_varstr_size(ref->db_accession);
    }
    if (ref->db_id_code && ref->db_id_code[0]) {
        const char* suf = dbref_id_suffix(meta->accession, ref->db_id_code);
        bytes += tlv_varstr_size(suf ? suf : ref->db_id_code);
    }
    if (!dbref_align_matches_chain(ref, chain)) {
        bytes += 2 + 8; // DBREF_TAG_DB_ALIGN
    }
    if (ref->seq_begin != 0 || ref->seq_end != 0) bytes += 2 + 8;
    if (ref->chain && ref->chain != 'A')          bytes += 2 + 1;
    if (has_nondefault_icodes(ref))               bytes += 2 + 4;
    if (ref->checksum != 0)                       bytes += 2 + 8;
    if (ref->ncbi_tax_id != 0)                    bytes += 2 + 4;

    if (ref->version_date && ref->version_date[0]) {
        uint8_t dt[3];
        if (pack_iso_date(ref->version_date, dt)) bytes += 2 + 3;
        else                                      bytes += tlv_varstr_size(ref->version_date);
    }
    bytes += tlv_varstr_size(ref->gene);
    bytes += tlv_varstr_size(ref->organism);
    bytes += tlv_varstr_size(ref->isoform);
    return bytes;
}

static void write_dbref_tlv(ByteWriter* w, const AfdbMetadata* meta, const ParsedChain* chain) {
    const DbRefEntry* ref = &meta->dbref;
    size_t payload_bytes = dbref_tlv_payload_bytes(meta, chain);
    bw_write_varlen(w, payload_bytes);

    if (ref->db_name && ref->db_name[0] && strcmp(ref->db_name, "UNP") != 0) {
        write_tlv_varstr(w, DBREF_TAG_DB_NAME, ref->db_name);
    }
    if (ref->db_accession && ref->db_accession[0] && strcmp(ref->db_accession, meta->accession) != 0) {
        write_tlv_varstr(w, DBREF_TAG_ACCESSION, ref->db_accession);
    }
    if (ref->db_id_code && ref->db_id_code[0]) {
        const char* suf = dbref_id_suffix(meta->accession, ref->db_id_code);
        if (suf) write_tlv_varstr(w, DBREF_TAG_ID_SUFFIX, suf);
        else     write_tlv_varstr(w, DBREF_TAG_ID_CODE,   ref->db_id_code);
    }

    if (!dbref_align_matches_chain(ref, chain)) {
        bw_write_u8(w, DBREF_TAG_DB_ALIGN);
        bw_write_u8(w, 8);
        bw_write(w, &ref->db_seq_begin, 4);
        bw_write(w, &ref->db_seq_end, 4);
    }
    if (ref->seq_begin != 0 || ref->seq_end != 0) {
        bw_write_u8(w, DBREF_TAG_SEQ_ALIGN);
        bw_write_u8(w, 8);
        bw_write(w, &ref->seq_begin, 4);
        bw_write(w, &ref->seq_end, 4);
    }
    if (ref->chain && ref->chain != 'A') {
        bw_write_u8(w, DBREF_TAG_CHAIN);
        bw_write_u8(w, 1);
        bw_write_u8(w, (uint8_t)ref->chain);
    }
    if (has_nondefault_icodes(ref)) {
        char ic[4] = { ref->iCode_begin, ref->iCode_end, ref->db_iCode_begin, ref->db_iCode_end };
        bw_write_u8(w, DBREF_TAG_ICODES);
        bw_write_u8(w, 4);
        bw_write(w, ic, 4);
    }
    if (ref->checksum != 0) {
        bw_write_u8(w, DBREF_TAG_CHECKSUM);
        bw_write_u8(w, 8);
        bw_write(w, &ref->checksum, 8);
    }
    if (ref->ncbi_tax_id != 0) {
        bw_write_u8(w, DBREF_TAG_TAX_ID);
        bw_write_u8(w, 4);
        bw_write(w, &ref->ncbi_tax_id, 4);
    }
    if (ref->version_date && ref->version_date[0]) {
        uint8_t dt[3];
        if (pack_iso_date(ref->version_date, dt)) {
            bw_write_u8(w, DBREF_TAG_VERSION_DATE);
            bw_write_u8(w, 3);
            bw_write(w, dt, 3);
        } else {
            write_tlv_varstr(w, DBREF_TAG_VERSION_DATE, ref->version_date);
        }
    }
    write_tlv_varstr(w, DBREF_TAG_GENE,     ref->gene);
    write_tlv_varstr(w, DBREF_TAG_ORGANISM, ref->organism);
    write_tlv_varstr(w, DBREF_TAG_ISOFORM,  ref->isoform);
    bw_write_u8(w, DBREF_TAG_END);
}

static bool read_dbref_tlv(ByteReader* r, const ChainView* chain, AfdbMetadata* out_meta) {
    size_t payload_bytes = 0;
    if (!br_read_varlen(r, &payload_bytes)) return false;
    if (payload_bytes == 0 || payload_bytes > br_remaining(r)) {
        r->error = true;
        return false;
    }

    ByteReader sub = br_init(r->ptr, payload_bytes);
    r->ptr += payload_bytes;

    DbRefEntry* ref = &out_meta->dbref;
    ref->chain          = 'A';
    ref->iCode_begin    = '?';
    ref->iCode_end      = '?';
    ref->db_iCode_begin = '?';
    ref->db_iCode_end   = '?';
    ref->db_seq_begin   = chain->first_seq_id;
    ref->db_seq_end     = chain->first_seq_id + (int32_t)chain->total_res - 1;

    while (true) {
        uint8_t tag = 0;
        if (!br_read(&sub, &tag, 1)) { r->error = true; return false; }
        if (tag == DBREF_TAG_END) break;

        size_t len = 0;
        if (!br_read_varlen(&sub, &len)) { r->error = true; return false; }

        switch (tag) {
            case DBREF_TAG_DB_NAME:
                free(ref->db_name);
                ref->db_name = br_read_str(&sub, len, false);
                break;
            case DBREF_TAG_ACCESSION:
                free(ref->db_accession);
                ref->db_accession = br_read_str(&sub, len, false);
                break;
            case DBREF_TAG_ID_CODE:
                free(ref->db_id_code);
                ref->db_id_code = br_read_str(&sub, len, false);
                break;
            case DBREF_TAG_ID_SUFFIX: {
                if (len > br_remaining(&sub)) { r->error = true; return false; }
                size_t alen = strlen(out_meta->accession);
                char* full = malloc(alen + 1 + len + 1);
                if (!full) { r->error = true; return false; }
                memcpy(full, out_meta->accession, alen);
                full[alen] = '_';
                if (len > 0) memcpy(full + alen + 1, sub.ptr, len);
                full[alen + 1 + len] = '\0';
                sub.ptr += len;
                free(ref->db_id_code);
                ref->db_id_code = full;
                break;
            }
            case DBREF_TAG_DB_ALIGN:
                if (len != 8 || !br_read(&sub, &ref->db_seq_begin, 4) || !br_read(&sub, &ref->db_seq_end, 4)) {
                    r->error = true; return false;
                }
                break;
            case DBREF_TAG_SEQ_ALIGN:
                if (len != 8 || !br_read(&sub, &ref->seq_begin, 4) || !br_read(&sub, &ref->seq_end, 4)) {
                    r->error = true; return false;
                }
                break;
            case DBREF_TAG_CHAIN:
                if (len != 1 || !br_read(&sub, &ref->chain, 1)) { r->error = true; return false; }
                break;
            case DBREF_TAG_ICODES: {
                char ic[4];
                if (len != 4 || !br_read(&sub, ic, 4)) { r->error = true; return false; }
                ref->iCode_begin    = ic[0];
                ref->iCode_end      = ic[1];
                ref->db_iCode_begin = ic[2];
                ref->db_iCode_end   = ic[3];
                break;
            }
            case DBREF_TAG_CHECKSUM:
                if (len != 8 || !br_read(&sub, &ref->checksum, 8)) { r->error = true; return false; }
                break;
            case DBREF_TAG_TAX_ID:
                if (len != 4 || !br_read(&sub, &ref->ncbi_tax_id, 4)) { r->error = true; return false; }
                break;
            case DBREF_TAG_VERSION_DATE:
                free(ref->version_date);
                if (len == 3) {
                    uint8_t dt[3];
                    if (!br_read(&sub, dt, 3)) { r->error = true; return false; }
                    ref->version_date = malloc(16);
                    if (!ref->version_date) { r->error = true; return false; }
                    snprintf(ref->version_date, 16, "%04u-%02u-%02u",
                             (unsigned)dt[0] + 1973U, (unsigned)(dt[1] % 100U), (unsigned)(dt[2] % 100U));
                } else {
                    ref->version_date = br_read_str(&sub, len, false);
                }
                break;
            case DBREF_TAG_GENE:
                free(ref->gene);
                ref->gene = br_read_str(&sub, len, false);
                break;
            case DBREF_TAG_ORGANISM:
                free(ref->organism);
                ref->organism = br_read_str(&sub, len, false);
                break;
            case DBREF_TAG_ISOFORM:
                free(ref->isoform);
                ref->isoform = br_read_str(&sub, len, false);
                break;
            default:
                if (len > br_remaining(&sub)) { r->error = true; return false; }
                sub.ptr += len;
                break;
        }
        if (sub.error) { r->error = true; return false; }
    }

    if (!ref->db_name) {
        ref->db_name = str_dup_n("UNP", 3);
        if (!ref->db_name) { r->error = true; return false; }
    }
    if (!ref->db_accession && out_meta->accession[0] != '\0') {
        ref->db_accession = str_dup_n(out_meta->accession, strlen(out_meta->accession));
        if (!ref->db_accession) { r->error = true; return false; }
    }
    return true;
}

static uint8_t clamp_tlv_len(const char* s) {
    if (!s) return 0;
    size_t len = strlen(s);
    return (len > 255) ? 255 : (uint8_t)len;
}

static size_t templates_tlv_payload_bytes(const AfdbMetadata* meta) {
    size_t bytes = bw_varlen_size(meta->template_count);
    for (size_t i = 0; i < meta->template_count; i++) {
        const TemplateEntry* t = &meta->templates[i];
        if (t->template_id != (uint16_t)(i + 1)) bytes += 4;
        uint8_t acc_len  = clamp_tlv_len(t->accession);
        uint8_t db_len   = (t->db_name && strcmp(t->db_name, "PDB") != 0) ? clamp_tlv_len(t->db_name) : 0;
        uint8_t asym_len = (t->asym_id && strcmp(t->asym_id, "A") != 0)   ? clamp_tlv_len(t->asym_id) : 0;
        if (acc_len > 0)  bytes += 2 + acc_len;
        if (db_len > 0)   bytes += 2 + db_len;
        if (asym_len > 0) bytes += 2 + asym_len;
        bytes += 1; // TPL_TAG_END
    }
    return bytes;
}

static size_t struct_afdb_encoded_size(const AfdbMetadata* meta, const ParsedChain* chain) {
    size_t bytes = bw_varlen_size(meta->conf_count);
    uint32_t prev_end = 0;
    for (size_t i = 0; i < meta->conf_count; i++) {
        const MetaStructConf* c = &meta->confs[i];
        uint32_t start_idx = c->start.resSeq - (uint32_t)chain->first_seq_id;
        uint32_t end_idx   = (c->end.resSeq >= c->start.resSeq)
                           ? (c->end.resSeq - (uint32_t)chain->first_seq_id) : start_idx;
        uint32_t gap  = (start_idx >= prev_end) ? (start_idx - prev_end) : 0;
        uint32_t diff = end_idx - start_idx;
        bytes += bw_varlen_size(gap) + bw_varlen_size(diff) + 1U;
        prev_end = end_idx;
    }
    return bytes;
}

size_t meta_cold_encoded_size(const AfdbMetadata* meta, const ParsedChain* chain) {
    size_t size = 8; // uint8_t fmt_widths[8]

    if (meta->flags & FLAG_RESOLUTION) size += 2; // uint16_t global_plddt x100

    if (meta->flags & FLAG_ENTITY) {
        if (!meta->entity_desc || strcmp(meta->entity_desc, UNCHAR_PROTEIN) == 0) {
            size += 1; // 0x00 sentinel
        } else {
            size_t len = strlen(meta->entity_desc);
            if (len > 65535) len = 65535;
            size += bw_varlen_size(len) + len;
        }
    }

    if (meta->flags & FLAG_STRUCT) {
        if (meta->flags & FLAG_AFDB) size += struct_afdb_encoded_size(meta, chain);
        else                         size += 8 + (meta->conf_count * sizeof(MetaStructConf));
    }
    if (meta->flags & FLAG_TEMPLATES) {
        size_t plen = templates_tlv_payload_bytes(meta);
        size += bw_varlen_size(plen) + plen;
    }
    if (meta->flags & FLAG_DBREF) {
        size_t plen = dbref_tlv_payload_bytes(meta, chain);
        size += bw_varlen_size(plen) + plen;
    }
    if (meta->flags & FLAG_AUDIT) {
        size += 1; // uint8_t audit_mask
        if (meta->audit_mask == AUDIT_MASK_RAW) {
            size += 4;
            for (uint32_t i = 0; i < meta->revision_count; i++) {
                size += sizeof(MetaAuditRevision)
                      + meta->revisions[i].meta.desc_len
                      + meta->revisions[i].meta.type_len;
            }
        }
    }
    if (meta->flags & FLAG_NEG_ZERO) {
        size += 4 + ((size_t)meta->neg_zero_count * 4U);
    }
    return size;
}

static void write_tlv_str(ByteWriter* w, uint8_t tag, const char* s) {
    uint8_t len = clamp_tlv_len(s);
    if (len == 0) return;
    bw_write_u8(w, tag);
    bw_write_u8(w, len);
    bw_write(w, s, len);
}

static void write_templates_tlv(ByteWriter* w, const AfdbMetadata* meta) {
    size_t payload_bytes = templates_tlv_payload_bytes(meta);
    bw_write_varlen(w, payload_bytes);
    bw_write_varlen(w, meta->template_count);

    for (size_t i = 0; i < meta->template_count; i++) {
        const TemplateEntry* t = &meta->templates[i];
        if (t->template_id != (uint16_t)(i + 1)) {
            bw_write_u8(w, TPL_TAG_ID);
            bw_write_u8(w, 2);
            bw_write(w, &t->template_id, 2);
        }
        write_tlv_str(w, TPL_TAG_ACCESSION, t->accession);
        if (t->db_name && strcmp(t->db_name, "PDB") != 0) {
            write_tlv_str(w, TPL_TAG_DB_NAME, t->db_name);
        }
        if (t->asym_id && strcmp(t->asym_id, "A") != 0) {
            write_tlv_str(w, TPL_TAG_ASYM_ID, t->asym_id);
        }
        bw_write_u8(w, TPL_TAG_END);
    }
}

void meta_write_cold(ByteWriter* w, const AfdbMetadata* meta, const ParsedChain* chain) {
    bw_write(w, chain->fmt_widths, 8);

    if (meta->flags & FLAG_RESOLUTION) {
        bw_write(w, &meta->global_plddt, 2);
    }
    if (meta->flags & FLAG_ENTITY) {
        if (!meta->entity_desc || strcmp(meta->entity_desc, UNCHAR_PROTEIN) == 0) {
            bw_write_u8(w, 0x00);
        } else {
            size_t dlen = strlen(meta->entity_desc);
            if (dlen > 65535) dlen = 65535;
            bw_write_varlen(w, dlen);
            bw_write(w, meta->entity_desc, dlen);
        }
    }
    if (meta->flags & FLAG_STRUCT) {
        if (meta->flags & FLAG_AFDB) {
            bw_write_varlen(w, meta->conf_count);
            uint32_t prev_end = 0;
            for (size_t i = 0; i < meta->conf_count; i++) {
                const MetaStructConf* c = &meta->confs[i];
                uint32_t start_idx = c->start.resSeq - (uint32_t)chain->first_seq_id;
                uint32_t end_idx   = (c->end.resSeq >= c->start.resSeq)
                                   ? (c->end.resSeq - (uint32_t)chain->first_seq_id) : start_idx;
                uint32_t gap  = (start_idx >= prev_end) ? (start_idx - prev_end) : 0;
                uint32_t diff = end_idx - start_idx;
                bw_write_varlen(w, gap);
                bw_write_varlen(w, diff);
                bw_write_u8(w, c->type_id);
                prev_end = end_idx;
            }
        } else {
            uint64_t count = (uint64_t)meta->conf_count;
            bw_write(w, &count, 8);
            bw_write(w, meta->confs, meta->conf_count * sizeof(MetaStructConf));
        }
    }
    if (meta->flags & FLAG_TEMPLATES) {
        write_templates_tlv(w, meta);
    }
    if (meta->flags & FLAG_DBREF) {
        write_dbref_tlv(w, meta, chain);
    }
    if (meta->flags & FLAG_AUDIT) {
        bw_write_u8(w, meta->audit_mask);
        if (meta->audit_mask == AUDIT_MASK_RAW) {
            bw_write(w, &meta->revision_count, 4);
            for (uint32_t i = 0; i < meta->revision_count; i++) {
                const AuditRevisionEntry* rev = &meta->revisions[i];
                bw_write(w, &rev->meta, sizeof(MetaAuditRevision));
                if (rev->meta.desc_len > 0) bw_write(w, rev->description, rev->meta.desc_len);
                if (rev->meta.type_len > 0) bw_write(w, rev->type, rev->meta.type_len);
            }
        }
    }
    if (meta->flags & FLAG_NEG_ZERO) {
        bw_write(w, &meta->neg_zero_count, 4);
        if (meta->neg_zero_count > 0) {
            bw_write(w, meta->neg_zeros, (size_t)meta->neg_zero_count * 4U);
        }
    }
}

static bool read_templates_tlv(ByteReader* r, AfdbMetadata* out_meta) {
    size_t payload_bytes = 0;
    if (!br_read_varlen(r, &payload_bytes)) return false;
    if (payload_bytes == 0 || payload_bytes > br_remaining(r)) {
        r->error = true;
        return false;
    }

    ByteReader sub = br_init(r->ptr, payload_bytes);
    r->ptr += payload_bytes;

    size_t count = 0;
    if (!br_read_varlen(&sub, &count)) { r->error = true; return false; }
    if (count > br_remaining(&sub)) { r->error = true; return false; }

    out_meta->templates = calloc(count, sizeof(TemplateEntry));
    if (!out_meta->templates && count > 0) { r->error = true; return false; }
    out_meta->template_count = count;

    for (size_t i = 0; i < count; i++) {
        TemplateEntry* t = &out_meta->templates[i];
        t->template_id = (uint16_t)(i + 1);
        while (true) {
            uint8_t tag = 0;
            if (!br_read(&sub, &tag, 1)) { r->error = true; return false; }
            if (tag == TPL_TAG_END) break;

            uint8_t len = 0;
            if (!br_read(&sub, &len, 1)) { r->error = true; return false; }

            switch (tag) {
                case TPL_TAG_ID:
                    if (len != 2 || !br_read(&sub, &t->template_id, 2)) { r->error = true; return false; }
                    break;
                case TPL_TAG_ACCESSION:
                    free(t->accession);
                    t->accession = br_read_str(&sub, len, false);
                    break;
                case TPL_TAG_DB_NAME:
                    free(t->db_name);
                    t->db_name = br_read_str(&sub, len, false);
                    break;
                case TPL_TAG_ASYM_ID:
                    free(t->asym_id);
                    t->asym_id = br_read_str(&sub, len, false);
                    break;
                default:
                    if (len > br_remaining(&sub)) { r->error = true; return false; }
                    sub.ptr += len;
                    break;
            }
            if (sub.error) { r->error = true; return false; }
        }
        if (!t->db_name) {
            t->db_name = str_dup_n("PDB", 3);
            if (!t->db_name) { r->error = true; return false; }
        }
        if (!t->asym_id) {
            t->asym_id = str_dup_n("A", 1);
            if (!t->asym_id) { r->error = true; return false; }
        }
    }
    return true;
}

static bool read_audit_block(ByteReader* r, AfdbMetadata* out_meta) {
    uint8_t mask = 0;
    if (!br_read(r, &mask, 1)) return false;
    out_meta->audit_mask = mask;
    if (mask != AUDIT_MASK_RAW) return true;

    uint32_t count = 0;
    if (!br_read(r, &count, 4)) return false;
    if ((size_t)count * sizeof(MetaAuditRevision) > br_remaining(r)) {
        r->error = true;
        return false;
    }

    out_meta->revisions = calloc(count, sizeof(AuditRevisionEntry));
    if (!out_meta->revisions && count > 0) return false;
    out_meta->revision_count = count;

    for (uint32_t i = 0; i < count; i++) {
        AuditRevisionEntry* rev = &out_meta->revisions[i];
        if (!br_read(r, &rev->meta, sizeof(MetaAuditRevision))) return false;
        rev->description = br_read_str(r, rev->meta.desc_len, false);
        rev->type        = br_read_str(r, rev->meta.type_len, false);
        if (r->error) return false;
    }
    return true;
}

bool meta_read_cold(ByteReader* r, ChainView* chain, AfdbMetadata* out_meta) {
    if (!br_read(r, chain->fmt_widths, 8)) return false;

    if (out_meta->flags & FLAG_RESOLUTION) {
        if (!br_read(r, &out_meta->global_plddt, 2)) return false;
    }
    if (out_meta->flags & FLAG_ENTITY) {
        size_t dlen = 0;
        if (!br_read_varlen(r, &dlen)) return false;
        if (dlen == 0) {
            out_meta->entity_desc = str_dup_n(UNCHAR_PROTEIN, sizeof(UNCHAR_PROTEIN) - 1);
            if (!out_meta->entity_desc) { r->error = true; return false; }
        } else {
            out_meta->entity_desc = br_read_str(r, dlen, false);
        }
        if (r->error) return false;
    }
    if (out_meta->flags & FLAG_STRUCT) {
        if (out_meta->flags & FLAG_AFDB) {
            size_t count = 0;
            if (!br_read_varlen(r, &count)) return false;
            if (count * 3U > br_remaining(r)) { r->error = true; return false; }
            out_meta->confs = calloc(count, sizeof(MetaStructConf));
            if (!out_meta->confs && count > 0) { r->error = true; return false; }
            out_meta->conf_count = count;

            uint32_t prev_end = 0;
            for (size_t i = 0; i < count; i++) {
                size_t  gap = 0, diff = 0;
                uint8_t type_id = 0;
                if (!br_read_varlen(r, &gap) || !br_read_varlen(r, &diff) || !br_read(r, &type_id, 1)) {
                    return false;
                }
                uint32_t start_idx = prev_end + (uint32_t)gap;
                uint32_t end_idx   = start_idx + (uint32_t)diff;
                if (start_idx >= chain->total_res || end_idx >= chain->total_res) {
                    r->error = true;
                    return false;
                }
                MetaStructConf* c = &out_meta->confs[i];
                c->start.resType  = (uint8_t)(chain->seq[start_idx] - 0x40);
                c->start.chain[0] = chain->chain_id;
                c->start.resSeq   = (uint32_t)chain->first_seq_id + start_idx;
                c->end.resType    = (uint8_t)(chain->seq[end_idx] - 0x40);
                c->end.chain[0]   = chain->chain_id;
                c->end.resSeq     = (uint32_t)chain->first_seq_id + end_idx;
                c->length         = (uint16_t)(diff + 1U);
                c->type_id        = type_id;
                prev_end = end_idx;
            }
        } else {
            uint64_t count = 0;
            if (!br_read(r, &count, 8)) return false;
            size_t bytes = (size_t)count * sizeof(MetaStructConf);
            if (bytes > br_remaining(r)) return false;
            out_meta->confs = malloc(bytes);
            if (!out_meta->confs && bytes > 0) return false;
            if (!br_read(r, out_meta->confs, bytes)) return false;
            out_meta->conf_count = (size_t)count;
        }
    }
    if (out_meta->flags & FLAG_TEMPLATES) {
        if (!read_templates_tlv(r, out_meta)) return false;
    }
    if (out_meta->flags & FLAG_DBREF) {
        if (!read_dbref_tlv(r, chain, out_meta)) return false;
    }
    if (out_meta->flags & FLAG_AUDIT) {
        if (!read_audit_block(r, out_meta)) return false;
    }
    if (out_meta->flags & FLAG_NEG_ZERO) {
        uint32_t count = 0;
        if (!br_read(r, &count, 4)) return false;
        size_t bytes = (size_t)count * 4U;
        if (bytes > br_remaining(r)) { r->error = true; return false; }
        out_meta->neg_zeros = count ? malloc(bytes) : NULL;
        if (count && !out_meta->neg_zeros) { r->error = true; return false; }
        if (count && !br_read(r, out_meta->neg_zeros, bytes)) return false;
        out_meta->neg_zero_count = count;
    }
    return !r->error;
}
