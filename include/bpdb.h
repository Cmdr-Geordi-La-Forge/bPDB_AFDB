#ifndef BPDB_H
#define BPDB_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

size_t encode_afdb_to_binary(const char* cif_data, size_t cif_len, uint8_t** out_buf);
void   decode_binary_to_afdb(const uint8_t* buf, size_t buf_len, FILE* out_f);

#endif // BPDB_H
