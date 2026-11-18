#include "binary/meta_codec.h"
#include "util/str_util.h"
#include <stdlib.h>
#include <string.h>

void afdb_meta_free(AfdbMetadata* meta) {
    if (!meta) return;
    free(meta->entity_desc);
    free(meta->confs);
    free(meta->templates);
    free(meta->gene);
    free(meta->organism);
    free(meta->isoform);
    if (meta->revisions) {
        for (uint32_t i = 0; i < meta->revision_count; i++) {
            free(meta->revisions[i].description);
            free(meta->revisions[i].type);
        }
        free(meta->revisions);
    }
    memset(meta, 0, sizeof(*meta));
}

size_t meta_encoded_size(const AfdbMetadata* meta) {
    size_t size = 64 + 8; // 64-byte global header + 8-byte bitmask

    if (meta->flags & FLAG_RESOLUTION) {
        size += 4;
    }
    if (meta->flags & FLAG_ENTITY) {
        size += 4 + (meta->entity_desc ? strlen(meta->entity_desc) : 0);
    }
    if (meta->flags & FLAG_STRUCT) {
        size += 8 + (meta->conf_count * sizeof(MetaStructConf));
    }
    if (meta->flags & FLAG_TEMPLATES) {
        size += 8 + (meta->template_count * sizeof(MetaTemplate));
    }
    if (meta->flags & FLAG_DBREF) {
        size += sizeof(MetaDBRef) + sizeof(AFDBTargetExt)
              + meta->target_ext.gene_len
              + meta->target_ext.org_len
              + meta->target_ext.isoform_len;
    }
    if (meta->flags & FLAG_AUDIT) {
        size += 4;
        for (uint32_t i = 0; i < meta->revision_count; i++) {
            size += sizeof(MetaAuditRevision)
                  + meta->revisions[i].meta.desc_len + 1
                  + meta->revisions[i].meta.type_len + 1;
        }
    }
    return size;
}

static void write_global_header(ByteWriter* w, const AfdbMetadata* meta) {
    uint8_t hdr[64];
    memset(hdr, 0, sizeof(hdr));
    memcpy(hdr, "bPDB", 4);
    hdr[4] = 0x00;
    hdr[5] = 1;
    hdr[6] = 1;
    hdr[7] = 0;

    memcpy(hdr + 8, meta->accession, 12);

    uint16_t tot_mod = 0;
    uint16_t cur_mod = (meta->frag_num > 0) ? (uint16_t)(meta->frag_num - 1) : 0;
    memcpy(hdr + 20, &tot_mod, 2);
    memcpy(hdr + 22, &cur_mod, 2);

    bw_write(w, hdr, sizeof(hdr));
}

void meta_write_binary(ByteWriter* w, const AfdbMetadata* meta) {
    write_global_header(w, meta);
    bw_write(w, &meta->flags, 8);

    if (meta->flags & FLAG_RESOLUTION) {
        bw_write(w, &meta->global_plddt, 4);
    }

    if (meta->flags & FLAG_ENTITY) {
        uint32_t dlen = meta->entity_desc ? (uint32_t)strlen(meta->entity_desc) : 0;
        bw_write(w, &dlen, 4);
        if (dlen > 0) bw_write(w, meta->entity_desc, dlen);
    }

    if (meta->flags & FLAG_STRUCT) {
        uint64_t count = (uint64_t)meta->conf_count;
        bw_write(w, &count, 8);
        bw_write(w, meta->confs, meta->conf_count * sizeof(MetaStructConf));
    }

    if (meta->flags & FLAG_TEMPLATES) {
        uint64_t count = (uint64_t)meta->template_count;
        bw_write(w, &count, 8);
        bw_write(w, meta->templates, meta->template_count * sizeof(MetaTemplate));
    }

    if (meta->flags & FLAG_DBREF) {
        bw_write(w, &meta->dbref, sizeof(MetaDBRef));
        bw_write(w, &meta->target_ext, sizeof(AFDBTargetExt));
        if (meta->target_ext.gene_len > 0)    bw_write(w, meta->gene, meta->target_ext.gene_len);
        if (meta->target_ext.org_len > 0)     bw_write(w, meta->organism, meta->target_ext.org_len);
        if (meta->target_ext.isoform_len > 0) bw_write(w, meta->isoform, meta->target_ext.isoform_len);
    }

    if (meta->flags & FLAG_AUDIT) {
        bw_write(w, &meta->revision_count, 4);
        for (uint32_t i = 0; i < meta->revision_count; i++) {
            const AuditRevisionEntry* rev = &meta->revisions[i];
            bw_write(w, &rev->meta, sizeof(MetaAuditRevision));
            if (rev->meta.desc_len > 0) bw_write(w, rev->description, rev->meta.desc_len);
            bw_write_u8(w, '\0');
            if (rev->meta.type_len > 0) bw_write(w, rev->type, rev->meta.type_len);
            bw_write_u8(w, '\0');
        }
    }
}

static bool read_global_header(ByteReader* r, AfdbMetadata* out_meta) {
    uint8_t hdr[64];
    if (!br_read(r, hdr, sizeof(hdr))) return false;

    memcpy(out_meta->accession, hdr + 8, 12);
    out_meta->accession[12] = '\0';

    uint16_t cur_mod = 0;
    memcpy(&cur_mod, hdr + 22, 2);
    out_meta->frag_num = (int)cur_mod + 1;
    return true;
}

static bool read_dbref_block(ByteReader* r, AfdbMetadata* out_meta) {
    if (!br_read(r, &out_meta->dbref, sizeof(MetaDBRef))) return false;
    if (!br_read(r, &out_meta->target_ext, sizeof(AFDBTargetExt))) return false;

    if (out_meta->target_ext.gene_len > 0) {
        out_meta->gene = br_read_str(r, out_meta->target_ext.gene_len, false);
    }
    if (out_meta->target_ext.org_len > 0) {
        out_meta->organism = br_read_str(r, out_meta->target_ext.org_len, false);
    }
    if (out_meta->target_ext.isoform_len > 0) {
        out_meta->isoform = br_read_str(r, out_meta->target_ext.isoform_len, false);
    }
    return !r->error;
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
        rev->description = br_read_str(r, rev->meta.desc_len, true);
        rev->type        = br_read_str(r, rev->meta.type_len, true);
        if (r->error) return false;
    }
    return true;
}

bool meta_read_binary(ByteReader* r, AfdbMetadata* out_meta) {
    memset(out_meta, 0, sizeof(*out_meta));

    if (br_remaining(r) < 80) return false;
    if (!read_global_header(r, out_meta)) return false;
    if (!br_read(r, &out_meta->flags, 8)) return false;

    if (out_meta->flags & FLAG_RESOLUTION) {
        if (!br_read(r, &out_meta->global_plddt, 4)) return false;
    }

    if (out_meta->flags & FLAG_ENTITY) {
        uint32_t dlen = 0;
        if (!br_read(r, &dlen, 4)) return false;
        if (dlen > 0) {
            out_meta->entity_desc = br_read_str(r, dlen, false);
            if (r->error) return false;
        }
    }
    if (!out_meta->entity_desc) {
        out_meta->entity_desc = str_dup_n("Unknown", 7);
    }

    if (out_meta->flags & FLAG_STRUCT) {
        uint64_t count = 0;
        if (!br_read(r, &count, 8)) return false;
        size_t bytes = (size_t)count * sizeof(MetaStructConf);
        if (bytes > br_remaining(r)) return false;
        out_meta->confs = malloc(bytes);
        if (!out_meta->confs && bytes > 0) return false;
        if (!br_read(r, out_meta->confs, bytes)) return false;
        out_meta->conf_count = (size_t)count;
    }

    if (out_meta->flags & FLAG_TEMPLATES) {
        uint64_t count = 0;
        if (!br_read(r, &count, 8)) return false;
        size_t bytes = (size_t)count * sizeof(MetaTemplate);
        if (bytes > br_remaining(r)) return false;
        out_meta->templates = malloc(bytes);
        if (!out_meta->templates && bytes > 0) return false;
        if (!br_read(r, out_meta->templates, bytes)) return false;
        out_meta->template_count = (size_t)count;
    }

    if (out_meta->flags & FLAG_DBREF) {
        if (!read_dbref_block(r, out_meta)) return false;
    }

    if (out_meta->flags & FLAG_AUDIT) {
        if (!read_audit_block(r, out_meta)) return false;
    }

    return !r->error;
}
