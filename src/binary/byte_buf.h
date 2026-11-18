#ifndef BPDB_BYTE_BUF_H
#define BPDB_BYTE_BUF_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const uint8_t* ptr;
    const uint8_t* end;
    bool           error;
} ByteReader;

static inline ByteReader br_init(const uint8_t* buf, size_t len) {
    ByteReader r = { buf, buf + len, false };
    return r;
}

static inline size_t br_remaining(const ByteReader* r) {
    return r->error ? 0 : (size_t)(r->end - r->ptr);
}

static inline bool br_read(ByteReader* r, void* dst, size_t n) {
    if (r->error || (size_t)(r->end - r->ptr) < n) {
        r->error = true;
        return false;
    }
    if (n > 0) memcpy(dst, r->ptr, n);
    r->ptr += n;
    return true;
}

static inline char* br_read_str(ByteReader* r, size_t len, bool has_nul_pad) {
    size_t wire_len = len + (has_nul_pad ? 1 : 0);
    if (r->error || (size_t)(r->end - r->ptr) < wire_len) {
        r->error = true;
        return NULL;
    }
    char* s = malloc(len + 1);
    if (!s) {
        r->error = true;
        return NULL;
    }
    if (len > 0) memcpy(s, r->ptr, len);
    s[len] = '\0';
    r->ptr += wire_len;
    return s;
}

typedef struct {
    uint8_t* start;
    uint8_t* ptr;
    uint8_t* end;
} ByteWriter;

static inline ByteWriter bw_init(uint8_t* buf, size_t cap) {
    ByteWriter w = { buf, buf, buf + cap };
    return w;
}

static inline size_t bw_len(const ByteWriter* w) {
    return (size_t)(w->ptr - w->start);
}

static inline void bw_write(ByteWriter* w, const void* src, size_t n) {
    if (n > 0) {
        memcpy(w->ptr, src, n);
        w->ptr += n;
    }
}

static inline void bw_write_u8(ByteWriter* w, uint8_t v) {
    *w->ptr++ = v;
}

static inline void bw_write_zeros(ByteWriter* w, size_t n) {
    if (n > 0) {
        memset(w->ptr, 0, n);
        w->ptr += n;
    }
}

#endif // BPDB_BYTE_BUF_H
