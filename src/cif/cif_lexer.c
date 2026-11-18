#include "cif/cif_lexer.h"
#include "util/str_util.h"
#include <string.h>

const char* cif_skip_spaces(const char* ptr, const char* end) {
    while (ptr < end && (*ptr == ' ' || *ptr == '\t')) ptr++;
    return ptr;
}

const char* cif_skip_line(const char* ptr, const char* end) {
    while (ptr < end && *ptr != '\n') ptr++;
    if (ptr < end) ptr++;
    return ptr;
}

bool cif_is_block_end(const char* ptr, const char* end) {
    if (ptr >= end) return true;
    if (*ptr == '#' || *ptr == '_') return true;
    if ((size_t)(end - ptr) >= 5 && strncmp(ptr, "loop_", 5) == 0) return true;
    return false;
}

const char* cif_next_token(const char* ptr, const char* end, char* out, size_t out_size, bool* was_quoted) {
    ptr = cif_skip_spaces(ptr, end);
    if (was_quoted) *was_quoted = false;
    if (out && out_size > 0) out[0] = '\0';
    if (ptr >= end || *ptr == '\n' || *ptr == '\r') return ptr;

    size_t d = 0;
    if (*ptr == '"' || *ptr == '\'') {
        char q = *ptr++;
        if (was_quoted) *was_quoted = true;
        while (ptr < end && *ptr != q) {
            if (out && d + 1 < out_size) out[d++] = *ptr;
            ptr++;
        }
        if (ptr < end) ptr++;
    } else {
        while (ptr < end && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r') {
            if (out && d + 1 < out_size) out[d++] = *ptr;
            ptr++;
        }
    }

    if (out && out_size > 0) out[d] = '\0';
    return ptr;
}

void cif_extract_kv_string(const char* data, const char* end, const char* key, char* out, int max_len) {
    const char* ptr = find_token(data, end, key);
    if (!ptr) return;
    ptr += strlen(key);
    cif_next_token(ptr, end, out, (size_t)max_len, NULL);
}
