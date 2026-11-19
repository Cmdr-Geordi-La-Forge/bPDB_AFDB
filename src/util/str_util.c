#include "util/str_util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char* find_token(const char* ptr, const char* end, const char* token) {
    size_t len = strlen(token);
    if (len == 0 || ptr >= end || (size_t)(end - ptr) < len) {
        return NULL;
    }
    unsigned char first = (unsigned char)token[0];
    char          last  = token[len - 1];
    const char*   limit = end - len;
    const char*   p     = ptr;

    while (p <= limit) {
        p = (const char*)memchr(p, first, (size_t)(limit - p + 1));
        if (!p) return NULL;
        if (p[len - 1] == last &&
            (len <= 2 || memcmp(p + 1, token + 1, len - 2) == 0)) {
            return p;
        }
        p++;
    }
    return NULL;
}

void strip_whitespace_fixed(const char* src, int width, char* dst, int max_len) {
    int d = 0;
    for (int i = 0; i < width && src[i] != '\0'; i++) {
        char c = src[i];
        if (c != ' ' && c != '\n' && c != '\r') {
            if (d < max_len) {
                dst[d++] = c;
            }
        }
    }
    dst[d] = '\0';
}

char* str_dup_n(const char* src, size_t len) {
    char* copy = malloc(len + 1);
    if (!copy) return NULL;
    if (len > 0) memcpy(copy, src, len);
    copy[len] = '\0';
    return copy;
}

int32_t parse_number(const char* str, int width) {
    int32_t val = 0;
    int32_t sign = 1;
    for (int i = 0; i < width; i++) {
        char c = str[i];
        if (c == '\0' || c == '\n' || c == '\r' || c == ' ') continue;
        if (c == '-') sign = -1;
        else if (c >= '0' && c <= '9') val = val * 10 + (c - '0');
    }
    return sign * val;
}

int32_t parse_fixed_decimal2(const char* str, int width) {
    const char* p = str;
    const char* e = str + width;
    while (p < e && (*p == ' ' || *p == '\t')) p++;
    int neg = 0;
    if (p < e && *p == '-') { neg = 1; p++; }

    int32_t val = 0;
    int frac_digits = 0;
    int seen_dot = 0;
    while (p < e) {
        char c = *p++;
        unsigned d = (unsigned)(c - '0');
        if (d < 10U) {
            if (!seen_dot) {
                val = val * 10 + (int32_t)d;
            } else if (frac_digits < 2) {
                val = val * 10 + (int32_t)d;
                frac_digits++;
            }
        } else if (c == '.') {
            seen_dot = 1;
        } else {
            break;
        }
    }
    while (frac_digits < 2) {
        val *= 10;
        frac_digits++;
    }
    return neg ? -val : val;
}

uint64_t parse_hex_uint64(const char* str) {
    uint64_t val = 0;
    while (*str) {
        char c = *str++;
        if (c >= '0' && c <= '9')      val = (val << 4) | (uint64_t)(c - '0');
        else if (c >= 'A' && c <= 'F') val = (val << 4) | (uint64_t)(c - 'A' + 10);
        else if (c >= 'a' && c <= 'f') val = (val << 4) | (uint64_t)(c - 'a' + 10);
        else break;
    }
    return val;
}

void format_fixed_decimal(char* buf, int32_t val, int decimals) {
    uint32_t uval = (val < 0) ? (uint32_t)(-(int64_t)val) : (uint32_t)val;
    if (decimals == 3) {
        sprintf(buf, "%s%u.%03u", (val < 0) ? "-" : "", uval / 1000U, uval % 1000U);
    } else {
        sprintf(buf, "%s%u.%02u", (val < 0) ? "-" : "", uval / 100U, uval % 100U);
    }
}
