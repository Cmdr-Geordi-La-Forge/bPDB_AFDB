#ifndef BPDB_CHEM_DICT_H
#define BPDB_CHEM_DICT_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "bpdb_internal.h"

extern const ResidueDict DICT[];
extern const size_t      NUM_DICT_RES;
extern const char* const ATOM37_ORDER[37];

const ResidueDict* dict_find_aa(uint8_t one_letter);
int                dict_get_atom_idx(uint8_t one_letter, const char* atom_name, bool is_terminal);
uint8_t            res3_to_type_id(const char* src, int width);
const char*        type_id_to_res3(uint8_t res_type);

#endif // BPDB_CHEM_DICT_H
