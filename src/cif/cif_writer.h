#ifndef BPDB_CIF_WRITER_H
#define BPDB_CIF_WRITER_H

#include <stdio.h>
#include "bpdb_internal.h"

void cif_write_afdb_metadata(FILE* out_f, const AfdbMetadata* meta, const ChainView* chain);
void cif_write_atom_site(FILE* out_f, const AfdbMetadata* meta, const ChainView* chain);

#endif // BPDB_CIF_WRITER_H
