#include "bpdb.h"
#include "bpdb_internal.h"
#include "binary/byte_buf.h"
#include "binary/meta_codec.h"
#include "binary/coord_codec.h"
#include "cif/cif_writer.h"

void decode_binary_to_afdb(const uint8_t* buf, size_t buf_len, FILE* out_f) {
    ByteReader   r    = br_init(buf, buf_len);
    AfdbMetadata meta = {0};

    if (!meta_read_binary(&r, &meta)) {
        afdb_meta_free(&meta);
        return;
    }

    ChainView chain;
    if (coord_read_chain_view(&r, &chain)) {
        if ((meta.flags & FLAG_AFDB) && (meta.flags & FLAG_CIF)) {
            cif_write_afdb_metadata(out_f, &meta, &chain);
        }
        cif_write_atom_site(out_f, &chain, meta.accession);
        fprintf(out_f, "#\n");
    }

    afdb_meta_free(&meta);
}
