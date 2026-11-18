#ifndef BPDB_CIF_LEXER_H
#define BPDB_CIF_LEXER_H

#include <stdbool.h>
#include <stddef.h>

const char* cif_skip_spaces(const char* ptr, const char* end);
const char* cif_skip_line(const char* ptr, const char* end);
bool        cif_is_block_end(const char* ptr, const char* end);
const char* cif_next_token(const char* ptr, const char* end, char* out, size_t out_size, bool* was_quoted);
void        cif_extract_kv_string(const char* data, const char* end, const char* key, char* out, int max_len);

#endif // BPDB_CIF_LEXER_H
