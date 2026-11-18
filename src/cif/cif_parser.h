#ifndef BPDB_CIF_PARSER_H
#define BPDB_CIF_PARSER_H

#include <stdbool.h>
#include <stddef.h>
#include "bpdb_internal.h"

bool cif_parse_afdb(const char* cif_data, size_t cif_len,
                    AfdbMetadata* out_meta, ParsedChain* out_chain);

#endif // BPDB_CIF_PARSER_H
