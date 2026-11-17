#include "utils.h"
#include <string.h>

const char* find_token(const char* ptr, const char* end, const char* token) {
    size_t len = strlen(token);
    for (const char* p = ptr; p <= end - len; p++) {
        if (memcmp(p, token, len) == 0) return p;
    }
    return NULL;
}

void extract_string_fixed(const char* src, int width, char* dest, int max_len) {
    int d = 0;
    for (int i = 0; i < width; i++) {
        if (src[i] != ' ' && src[i] != '\n' && src[i] != '\r' && src[i] != '\0') {
            if (d < max_len) dest[d++] = src[i];
        }
    }
    dest[d] = '\0';
}

int32_t parse_number(const char* str, int width) {
    int32_t val = 0, sign = 1;
    for (int i = 0; i < width; i++) {
        char c = str[i];
        if (c == '\0' || c == '\n' || c == '\r' || c == ' ') continue;
        if (c == '-') sign = -1;
        else if (c >= '0' && c <= '9') val = val * 10 + (c - '0');
    }
    return sign * val;
}
