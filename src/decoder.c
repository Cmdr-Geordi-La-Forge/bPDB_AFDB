#include "bpdb.h"
#include "bpdb_internal.h"
#include "binary/byte_buf.h"
#include "binary/meta_codec.h"
#include "binary/coord_codec.h"
#include "cif/cif_writer.h"
#include <string.h>

void decode_binary_to_afdb(const uint8_t* buf, size_t buf_len, FILE* out_f) {
    BpdbV2Header hdr;
    ChainView    chain = {0};
    if (!coord_read_v2_hot(buf, buf_len, &hdr, &chain)) return;

    AfdbMetadata meta = {0};
    meta.flags    = hdr.flags;
    meta.frag_num = hdr.frag_num;
    memcpy(meta.accession, hdr.accession, 12);
    meta.accession[12] = '\0';

    ByteReader cold = br_init(buf + hdr.cold_offset, buf_len - hdr.cold_offset);
    if (meta_read_cold(&cold, &chain, &meta)) {
        if ((meta.flags & FLAG_AFDB) && (meta.flags & FLAG_CIF)) {
            cif_write_afdb_metadata(out_f, &meta, &chain);
        }
        cif_write_atom_site(out_f, &meta, &chain);
        fprintf(out_f, "#\n");
    }

    afdb_meta_free(&meta);
}
