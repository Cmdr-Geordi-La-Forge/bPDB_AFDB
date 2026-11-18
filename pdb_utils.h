#ifndef BPDB_PDB_UTILS_H
#define BPDB_PDB_UTILS_H

#include "bpdb.h"
#include <stdbool.h>
#include <stdio.h>

uint8_t res3_to_type_id(const char* src, int width);
const char* type_id_to_res3(uint8_t resType);

uint8_t encode_struct_conf_type(const char* type_str);
const char* decode_struct_conf_type(uint8_t type_id);

int get_atom_idx(uint8_t one_letter_code, const char* atom_name, bool is_terminal);

size_t parse_meta_struct_conf(const char* cif_data, const char* end, MetaStructConf** out_confs);
void write_meta_struct_conf(FILE* out_f, const MetaStructConf* confs, size_t count);

#endif // BPDB_PDB_UTILS_H
