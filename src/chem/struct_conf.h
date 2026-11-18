#ifndef BPDB_CHEM_STRUCT_CONF_H
#define BPDB_CHEM_STRUCT_CONF_H

#include <stdint.h>

uint8_t     encode_struct_conf_type(const char* type_str);
const char* decode_struct_conf_type(uint8_t type_id);

#endif // BPDB_CHEM_STRUCT_CONF_H
