#ifndef BPDB_UTILS_H
#define BPDB_UTILS_H

#include <stdint.h>

const char* find_token(const char* ptr, const char* end, const char* token);
void extract_string_fixed(const char* src, int width, char* dest, int max_len);
int32_t parse_number(const char* str, int width);

#endif
