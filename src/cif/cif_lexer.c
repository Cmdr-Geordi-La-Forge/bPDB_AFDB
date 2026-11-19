// src/cif/cif_lexer.c
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

bool cif_has_inline_value(const char* after_key, const char* end) {
    const char* p = cif_skip_spaces(after_key, end);
    if (p < end && *p != '\n' && *p != '\r') return true;
    // Also detect semicolon-delimited value starting on the next line
    if (p < end && *p == '\r') p++;
    if (p < end && *p == '\n') p++;
    return (p < end && *p == ';');
}

const char* cif_find_block_end(const char* start, const char* end) {
    const char* p = start;
    while (p < end) {
        if (*p == '#' && (p == start || *(p - 1) == '\n')) return p;
        p++;
    }
    return end;
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
        while (ptr < end) {
            if (*ptr == q && (ptr + 1 >= end || ptr[1] == ' ' || ptr[1] == '\t' || ptr[1] == '\n' || ptr[1] == '\r')) {
                ptr++;
                break;
            }
            if (out && d + 1 < out_size) out[d++] = *ptr;
            ptr++;
        }
    } else {
        while (ptr < end && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r') {
            if (out && d + 1 < out_size) out[d++] = *ptr;
            ptr++;
        }
    }

    if (out && out_size > 0) out[d] = '\0';
    return ptr;
}

static bool extract_semicolon_span(const char* ptr, const char* end, const char** out_start, size_t* out_len) {
    const char* p = ptr;
    if (p < end && *p == '\r') p++;
    if (p < end && *p == '\n') p++;
    if (p >= end || *p != ';') return false;
    p++; // skip opening ';'
    const char* start = p;
    while (p < end) {
        if (*p == ';' && *(p - 1) == '\n') {
            size_t len = (size_t)((p - 1) - start);
            if (len > 0 && start[len - 1] == '\r') len--;
            *out_start = start;
            *out_len   = len;
            return true;
        }
        p++;
    }
    return false;
}

void cif_extract_kv_string(const char* data, const char* end, const char* key, char* out, int max_len) {
    if (max_len <= 0) return;
    const char* ptr = find_token(data, end, key);
    if (!ptr) return;
    ptr = cif_skip_spaces(ptr + strlen(key), end);

    const char* semi_start = NULL;
    size_t semi_len = 0;
    if (extract_semicolon_span(ptr, end, &semi_start, &semi_len)) {
        size_t copy = (semi_len < (size_t)(max_len - 1)) ? semi_len : (size_t)(max_len - 1);
        memcpy(out, semi_start, copy);
        out[copy] = '\0';
        return;
    }
    cif_next_token(ptr, end, out, (size_t)max_len, NULL);
}

char* cif_extract_kv_alloc(const char* data, const char* end, const char* key) {
    const char* ptr = find_token(data, end, key);
    if (!ptr) return NULL;
    ptr = cif_skip_spaces(ptr + strlen(key), end);

    const char* semi_start = NULL;
    size_t semi_len = 0;
    if (extract_semicolon_span(ptr, end, &semi_start, &semi_len)) {
        return (semi_len > 0) ? str_dup_n(semi_start, semi_len) : NULL;
    }

    if (ptr >= end || *ptr == '\n' || *ptr == '\r') return NULL;
    char quote = 0;
    if (*ptr == '"' || *ptr == '\'') quote = *ptr++;

    const char* start = ptr;
    while (ptr < end) {
        if (quote) {
            if (*ptr == quote && (ptr + 1 >= end || ptr[1] == ' ' || ptr[1] == '\t' || ptr[1] == '\n' || ptr[1] == '\r')) break;
        } else {
            if (*ptr == ' ' || *ptr == '\t' || *ptr == '\n' || *ptr == '\r') break;
        }
        ptr++;
    }
    size_t len = (size_t)(ptr - start);
    return (len > 0) ? str_dup_n(start, len) : NULL;
}
