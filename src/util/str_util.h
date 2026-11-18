#ifndef BPDB_STR_UTIL_H
#define BPDB_STR_UTIL_H

#include <stddef.h>
#include <stdint.h>

const char* find_token(const char* ptr, const char* end, const char* token);
void        strip_whitespace_fixed(const char* src, int width, char* dst, int max_len);
char*       str_dup_n(const char* src, size_t len);
int32_t     parse_number(const char* str, int width);
uint64_t    parse_hex_uint64(const char* str);
void        format_fixed_decimal(char* buf, int32_t val, int decimals);

#endif // BPDB_STR_UTIL_H
