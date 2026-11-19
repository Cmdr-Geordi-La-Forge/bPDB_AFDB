#include "binary/meta_codec.h"
#include "util/str_util.h"
#include <stdlib.h>
#include <string.h>

void afdb_meta_free(AfdbMetadata* meta) {
    if (!meta) return;
    free(meta->entity_desc);
    free(meta->confs);
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
    free(meta->neg_zeros);
    memset(meta, 0, sizeof(*meta));
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

static bool has_nondefault_icodes(const DbRefEntry* ref) {
    return (ref->iCode_begin && ref->iCode_begin != '?') ||
           (ref->iCode_end && ref->iCode_end != '?') ||
           (ref->db_iCode_begin && ref->db_iCode_begin != '?') ||
           (ref->db_iCode_end && ref->db_iCode_end != '?');
}

static size_t dbref_tlv_payload_bytes(const DbRefEntry* ref) {
    size_t bytes = 1; // DBREF_TAG_END
    bytes += tlv_varstr_size(ref->db_name);
    bytes += tlv_varstr_size(ref->db_accession);
    bytes += tlv_varstr_size(ref->db_id_code);
    bytes += 2 + 8; // DBREF_TAG_DB_ALIGN
    if (ref->seq_begin != 0 || ref->seq_end != 0) bytes += 2 + 8;
    if (ref->chain && ref->chain != 'A')          bytes += 2 + 1;
    if (has_nondefault_icodes(ref))               bytes += 2 + 4;
    if (ref->checksum != 0)                       bytes += 2 + 8;
    if (ref->ncbi_tax_id != 0)                    bytes += 2 + 4;
    bytes += tlv_varstr_size(ref->version_date);
    bytes += tlv_varstr_size(ref->gene);
    bytes += tlv_varstr_size(ref->organism);
    bytes += tlv_varstr_size(ref->isoform);
    return bytes;
}

static void write_dbref_tlv(ByteWriter* w, const DbRefEntry* ref) {
    uint32_t payload_bytes = (uint32_t)dbref_tlv_payload_bytes(ref);
    bw_write(w, &payload_bytes, 4);

    write_tlv_varstr(w, DBREF_TAG_DB_NAME,   ref->db_name);
    write_tlv_varstr(w, DBREF_TAG_ACCESSION, ref->db_accession);
    write_tlv_varstr(w, DBREF_TAG_ID_CODE,   ref->db_id_code);

    bw_write_u8(w, DBREF_TAG_DB_ALIGN);
    bw_write_u8(w, 8);
    bw_write(w, &ref->db_seq_begin, 4);
    bw_write(w, &ref->db_seq_end, 4);

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

    write_tlv_varstr(w, DBREF_TAG_VERSION_DATE, ref->version_date);
    write_tlv_varstr(w, DBREF_TAG_GENE,         ref->gene);
    write_tlv_varstr(w, DBREF_TAG_ORGANISM,     ref->organism);
    write_tlv_varstr(w, DBREF_TAG_ISOFORM,      ref->isoform);
    bw_write_u8(w, DBREF_TAG_END);
}

static bool read_dbref_tlv(ByteReader* r, AfdbMetadata* out_meta) {
    uint32_t payload_bytes = 0;
    if (!br_read(r, &payload_bytes, 4)) return false;
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
                ref->version_date = br_read_str(&sub, len, false);
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
    return true;
}

static uint8_t clamp_tlv_len(const char* s) {
    if (!s) return 0;
    size_t len = strlen(s);
    return (len > 255) ? 255 : (uint8_t)len;
}

static size_t templates_tlv_payload_bytes(const AfdbMetadata* meta) {
    size_t bytes = 4; // uint32_t template_count
    for (size_t i = 0; i < meta->template_count; i++) {
        const TemplateEntry* t = &meta->templates[i];
        bytes += 4; // TPL_TAG_ID (1) + len (1) + uint16_t (2)
        uint8_t acc_len  = clamp_tlv_len(t->accession);
        uint8_t db_len   = clamp_tlv_len(t->db_name);
        uint8_t asym_len = clamp_tlv_len(t->asym_id);
        if (acc_len > 0)  bytes += 2 + acc_len;
        if (db_len > 0)   bytes += 2 + db_len;
        if (asym_len > 0) bytes += 2 + asym_len;
        bytes += 1; // TPL_TAG_END
    }
    return bytes;
}

size_t meta_cold_encoded_size(const AfdbMetadata* meta) {
    size_t size = 8; // uint8_t fmt_widths[8]

    if (meta->flags & FLAG_RESOLUTION) size += 4;
    if (meta->flags & FLAG_ENTITY)     size += 4 + (meta->entity_desc ? strlen(meta->entity_desc) : 0);
    if (meta->flags & FLAG_STRUCT) {
        if (meta->flags & FLAG_AFDB) size += 4 + (meta->conf_count * 5U);
        else                         size += 8 + (meta->conf_count * sizeof(MetaStructConf));
    }
    if (meta->flags & FLAG_TEMPLATES)  size += 4 + templates_tlv_payload_bytes(meta);
    if (meta->flags & FLAG_DBREF)      size += 4 + dbref_tlv_payload_bytes(&meta->dbref);

    if (meta->flags & FLAG_AUDIT) {
        size += 4;
        for (uint32_t i = 0; i < meta->revision_count; i++) {
            size += sizeof(MetaAuditRevision)
                  + meta->revisions[i].meta.desc_len
                  + meta->revisions[i].meta.type_len;
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
    uint32_t payload_bytes = (uint32_t)templates_tlv_payload_bytes(meta);
    uint32_t count         = (uint32_t)meta->template_count;
    bw_write(w, &payload_bytes, 4);
    bw_write(w, &count, 4);

    for (size_t i = 0; i < meta->template_count; i++) {
        const TemplateEntry* t = &meta->templates[i];
        bw_write_u8(w, TPL_TAG_ID);
        bw_write_u8(w, 2);
        bw_write(w, &t->template_id, 2);

        write_tlv_str(w, TPL_TAG_ACCESSION, t->accession);
        write_tlv_str(w, TPL_TAG_DB_NAME,   t->db_name);
        write_tlv_str(w, TPL_TAG_ASYM_ID,   t->asym_id);
        bw_write_u8(w, TPL_TAG_END);
    }
}

void meta_write_cold(ByteWriter* w, const AfdbMetadata* meta, const ParsedChain* chain) {
    bw_write(w, chain->fmt_widths, 8);

    if (meta->flags & FLAG_RESOLUTION) {
        bw_write(w, &meta->global_plddt, 4);
    }
    if (meta->flags & FLAG_ENTITY) {
        uint32_t dlen = meta->entity_desc ? (uint32_t)strlen(meta->entity_desc) : 0;
        bw_write(w, &dlen, 4);
        if (dlen > 0) bw_write(w, meta->entity_desc, dlen);
    }
    if (meta->flags & FLAG_STRUCT) {
        if (meta->flags & FLAG_AFDB) {
            uint32_t count = (uint32_t)meta->conf_count;
            bw_write(w, &count, 4);
            for (size_t i = 0; i < meta->conf_count; i++) {
                const MetaStructConf* c = &meta->confs[i];
                uint16_t start_idx = (uint16_t)(c->start.resSeq - (uint32_t)chain->first_seq_id);
                uint16_t length    = (c->end.resSeq >= c->start.resSeq) ? (uint16_t)(c->end.resSeq - c->start.resSeq + 1U) : 0;
                bw_write(w, &start_idx, 2);
                bw_write(w, &length, 2);
                bw_write_u8(w, c->type_id);
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
        write_dbref_tlv(w, &meta->dbref);
    }
    if (meta->flags & FLAG_AUDIT) {
        bw_write(w, &meta->revision_count, 4);
        for (uint32_t i = 0; i < meta->revision_count; i++) {
            const AuditRevisionEntry* rev = &meta->revisions[i];
            bw_write(w, &rev->meta, sizeof(MetaAuditRevision));
            if (rev->meta.desc_len > 0) bw_write(w, rev->description, rev->meta.desc_len);
            if (rev->meta.type_len > 0) bw_write(w, rev->type, rev->meta.type_len);
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
    uint32_t payload_bytes = 0;
    if (!br_read(r, &payload_bytes, 4)) return false;
    if (payload_bytes < 4 || payload_bytes > br_remaining(r)) {
        r->error = true;
        return false;
    }

    ByteReader sub = br_init(r->ptr, payload_bytes);
    r->ptr += payload_bytes;

    uint32_t count = 0;
    if (!br_read(&sub, &count, 4)) { r->error = true; return false; }

    out_meta->templates = calloc(count, sizeof(TemplateEntry));
    if (!out_meta->templates && count > 0) { r->error = true; return false; }
    out_meta->template_count = count;

    for (uint32_t i = 0; i < count; i++) {
        TemplateEntry* t = &out_meta->templates[i];
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
                    // Forward-compatible skip of unknown TLV fields
                    if (len > br_remaining(&sub)) { r->error = true; return false; }
                    sub.ptr += len;
                    break;
            }
            if (sub.error) { r->error = true; return false; }
        }
    }
    return true;
}

static bool read_audit_block(ByteReader* r, AfdbMetadata* out_meta) {
    uint32_t count = 0;
    if (!br_read(r, &count, 4)) return false;

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
        if (!br_read(r, &out_meta->global_plddt, 4)) return false;
    }
    if (out_meta->flags & FLAG_ENTITY) {
        uint32_t dlen = 0;
        if (!br_read(r, &dlen, 4)) return false;
        out_meta->entity_desc = br_read_str(r, dlen, false);
        if (r->error) return false;
    }
    if (out_meta->flags & FLAG_STRUCT) {
        if (out_meta->flags & FLAG_AFDB) {
            uint32_t count = 0;
            if (!br_read(r, &count, 4)) return false;
            if ((size_t)count * 5U > br_remaining(r)) { r->error = true; return false; }
            out_meta->confs = calloc(count, sizeof(MetaStructConf));
            if (!out_meta->confs && count > 0) { r->error = true; return false; }
            out_meta->conf_count = count;

            for (uint32_t i = 0; i < count; i++) {
                uint16_t start_idx = 0, length = 0;
                uint8_t  type_id = 0;
                if (!br_read(r, &start_idx, 2) || !br_read(r, &length, 2) || !br_read(r, &type_id, 1)) {
                    return false;
                }
                uint32_t end_idx = (length > 0) ? ((uint32_t)start_idx + length - 1U) : start_idx;
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
                c->length         = length;
                c->type_id        = type_id;
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
        if (!read_dbref_tlv(r, out_meta)) return false;
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
