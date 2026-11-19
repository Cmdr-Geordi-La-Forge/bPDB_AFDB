#define _POSIX_C_SOURCE 200809L
#define _FILE_OFFSET_BITS 64

#include "gcs_stream.h"
#include <archive.h>
#include <archive_entry.h>
#include <curl/curl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define TOKEN_REFRESH_SEC    2400
#define SMALL_TAR_MEM_LIMIT  (16ULL * 1024ULL * 1024ULL)
#define MAX_CIF_MEMBER_BYTES (256ULL * 1024ULL * 1024ULL)
#define MAX_HTTP_RETRIES     5

/* ===================== OAuth Token Manager ===================== */

static bool fetch_gcloud_token_locked(TokenManager* tm) {
    const char* env_tok = getenv("BPDB_GCS_TOKEN");
    if (env_tok && env_tok[0] != '\0') {
        snprintf(tm->token, sizeof(tm->token), "%s", env_tok);
        tm->fetched_at = time(NULL);
        tm->generation++;
        return true;
    }

    FILE* fp = popen("gcloud auth print-access-token 2>/dev/null", "r");
    if (!fp) return false;

    char buf[2048] = {0};
    size_t n = fread(buf, 1, sizeof(buf) - 1, fp);
    int rc = pclose(fp);
    if (rc != 0 || n == 0) return false;

    while (n > 0 && (buf[n - 1] == '\n' || buf[n - 1] == '\r' || buf[n - 1] == ' ' || buf[n - 1] == '\t')) {
        buf[--n] = '\0';
    }
    if (n < 16) return false;

    memcpy(tm->token, buf, n + 1);
    tm->fetched_at = time(NULL);
    tm->generation++;
    return true;
}

bool token_manager_init(TokenManager* tm) {
    memset(tm, 0, sizeof(*tm));
    pthread_mutex_init(&tm->mu, NULL);
    pthread_mutex_lock(&tm->mu);
    bool ok = fetch_gcloud_token_locked(tm);
    pthread_mutex_unlock(&tm->mu);
    return ok;
}

void token_manager_destroy(TokenManager* tm) {
    if (!tm) return;
    pthread_mutex_destroy(&tm->mu);
    memset(tm, 0, sizeof(*tm));
}

uint64_t token_manager_get(TokenManager* tm, char* dst, size_t dst_cap) {
    if (!tm || !dst || dst_cap == 0) return 0;
    pthread_mutex_lock(&tm->mu);
    time_t now = time(NULL);
    if (tm->token[0] == '\0' || (now - tm->fetched_at) >= TOKEN_REFRESH_SEC) {
        fetch_gcloud_token_locked(tm);
    }
    snprintf(dst, dst_cap, "%s", tm->token);
    uint64_t gen = tm->generation;
    pthread_mutex_unlock(&tm->mu);
    return gen;
}

bool token_manager_refresh(TokenManager* tm, uint64_t failed_gen) {
    if (!tm) return false;
    pthread_mutex_lock(&tm->mu);
    bool ok = true;
    if (tm->generation == failed_gen) {
        ok = fetch_gcloud_token_locked(tm);
    }
    pthread_mutex_unlock(&tm->mu);
    return ok;
}

/* ===================== Per-Worker Streamer ===================== */

struct GcsStreamer {
    CURL*    small_easy;
    CURLM*   multi;
    CURL*    stream_easy;
    uint8_t* mem_buf;
    size_t   mem_cap;
    uint8_t* entry_buf;
    size_t   entry_cap;
};

static bool ensure_buf(uint8_t** buf, size_t* cap, size_t need) {
    if (need <= *cap) return true;
    size_t ncap = (*cap == 0) ? (256U * 1024U) : *cap;
    while (ncap < need) ncap *= 2U;
    uint8_t* nb = realloc(*buf, ncap);
    if (!nb) return false;
    *buf = nb;
    *cap = ncap;
    return true;
}

GcsStreamer* gcs_streamer_create(void) {
    GcsStreamer* gs = calloc(1, sizeof(*gs));
    if (!gs) return NULL;
    gs->small_easy  = curl_easy_init();
    gs->multi       = curl_multi_init();
    gs->stream_easy = curl_easy_init();
    if (!gs->small_easy || !gs->multi || !gs->stream_easy) {
        gcs_streamer_destroy(gs);
        return NULL;
    }
    ensure_buf(&gs->mem_buf, &gs->mem_cap, 512U * 1024U);
    ensure_buf(&gs->entry_buf, &gs->entry_cap, 256U * 1024U);
    return gs;
}

void gcs_streamer_destroy(GcsStreamer* gs) {
    if (!gs) return;
    if (gs->stream_easy) curl_easy_cleanup(gs->stream_easy);
    if (gs->multi)       curl_multi_cleanup(gs->multi);
    if (gs->small_easy)  curl_easy_cleanup(gs->small_easy);
    free(gs->mem_buf);
    free(gs->entry_buf);
    free(gs);
}

/* ===================== libarchive Member Loop ===================== */

static bool has_suffix(const char* s, const char* suf) {
    size_t sl = strlen(s), sufl = strlen(suf);
    return (sl >= sufl) && (memcmp(s + sl - sufl, suf, sufl) == 0);
}

static bool is_target_cif_name(const char* name) {
    if (!name) return false;
    const char* base = strrchr(name, '/');
    base = base ? (base + 1) : name;
    if (base[0] == '\0' || (base[0] == '.' && base[1] == '_')) return false;
    return has_suffix(base, ".cif.gz") || has_suffix(base, ".cif");
}

static bool extract_archive_cifs(GcsStreamer* gs, struct archive* a,
                                 TarCifCallback cb, void* user_data,
                                 uint64_t* out_cif_count, uint64_t* out_skipped) {
    struct archive_entry* entry = NULL;
    uint64_t cifs = 0, skipped = 0;

    while (true) {
        int r = archive_read_next_header(a, &entry);
        if (r == ARCHIVE_EOF) break;
        if (r != ARCHIVE_OK && r != ARCHIVE_WARN) return false;

        const char* name = archive_entry_pathname(entry);
        bool is_reg = (archive_entry_filetype(entry) == AE_IFREG);
        int64_t sz  = archive_entry_size_is_set(entry) ? archive_entry_size(entry) : -1;

        if (!is_reg || !is_target_cif_name(name) || sz == 0 ||
            (sz > 0 && (uint64_t)sz > MAX_CIF_MEMBER_BYTES)) {
            skipped++;
            if (archive_read_data_skip(a) != ARCHIVE_OK) return false;
            continue;
        }

        size_t total = 0;
        if (sz > 0 && !ensure_buf(&gs->entry_buf, &gs->entry_cap, (size_t)sz + 1)) return false;

        while (true) {
            if (!ensure_buf(&gs->entry_buf, &gs->entry_cap, total + 65536U + 1U)) return false;
            la_ssize_t got = archive_read_data(a, gs->entry_buf + total, gs->entry_cap - 1U - total);
            if (got < 0) return false;
            if (got == 0) break;
            total += (size_t)got;
            if (total > MAX_CIF_MEMBER_BYTES) return false;
        }
        gs->entry_buf[total] = '\0';

        if (!cb(name, gs->entry_buf, total, user_data)) return false;
        cifs++;
    }

    if (out_cif_count) *out_cif_count += cifs;
    if (out_skipped)   *out_skipped   += skipped;
    return true;
}

/* ===================== Small-Tar Buffered Download ===================== */

typedef struct {
    GcsStreamer* gs;
    size_t       len;
    bool         oom;
} MemWriteCtx;

static size_t mem_write_cb(char* ptr, size_t size, size_t nmemb, void* userdata) {
    MemWriteCtx* ctx = (MemWriteCtx*)userdata;
    size_t n = size * nmemb;
    if (n == 0) return 0;
    if (!ensure_buf(&ctx->gs->mem_buf, &ctx->gs->mem_cap, ctx->len + n)) {
        ctx->oom = true;
        return 0;
    }
    memcpy(ctx->gs->mem_buf + ctx->len, ptr, n);
    ctx->len += n;
    return n;
}

static struct curl_slist* make_auth_header(TokenManager* tm, bool use_auth, uint64_t* out_gen) {
    *out_gen = 0;
    if (!use_auth || !tm) return NULL;
    char tok[2048];
    *out_gen = token_manager_get(tm, tok, sizeof(tok));
    if (tok[0] == '\0') return NULL;
    char hdr[2100];
    snprintf(hdr, sizeof(hdr), "Authorization: Bearer %s", tok);
    return curl_slist_append(NULL, hdr);
}

static void configure_common_curl(CURL* easy, const char* url, struct curl_slist* hdrs) {
    curl_easy_reset(easy);
    curl_easy_setopt(easy, CURLOPT_URL, url);
    if (hdrs) curl_easy_setopt(easy, CURLOPT_HTTPHEADER, hdrs);
    curl_easy_setopt(easy, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(easy, CURLOPT_HTTP_VERSION, (long)CURL_HTTP_VERSION_2TLS);
    curl_easy_setopt(easy, CURLOPT_NOSIGNAL, 1L);
    curl_easy_setopt(easy, CURLOPT_TCP_KEEPALIVE, 1L);
    curl_easy_setopt(easy, CURLOPT_CONNECTTIMEOUT, 30L);
    curl_easy_setopt(easy, CURLOPT_LOW_SPEED_LIMIT, 1024L);
    curl_easy_setopt(easy, CURLOPT_LOW_SPEED_TIME, 90L);
    curl_easy_setopt(easy, CURLOPT_BUFFERSIZE, 262144L);
}

static bool fetch_small_tar_and_extract(GcsStreamer* gs, TokenManager* tm,
                                        const char* url, bool use_auth,
                                        TarCifCallback cb, void* user_data,
                                        uint64_t* out_cif_count, uint64_t* out_skipped) {
    for (int attempt = 0; attempt < MAX_HTTP_RETRIES; attempt++) {
        uint64_t tok_gen = 0;
        struct curl_slist* hdrs = make_auth_header(tm, use_auth, &tok_gen);
        MemWriteCtx wctx = { .gs = gs, .len = 0, .oom = false };

        configure_common_curl(gs->small_easy, url, hdrs);
        curl_easy_setopt(gs->small_easy, CURLOPT_WRITEFUNCTION, mem_write_cb);
        curl_easy_setopt(gs->small_easy, CURLOPT_WRITEDATA, &wctx);

        CURLcode res = curl_easy_perform(gs->small_easy);
        long http_code = 0;
        curl_easy_getinfo(gs->small_easy, CURLINFO_RESPONSE_CODE, &http_code);
        curl_slist_free_all(hdrs);

        if (res == CURLE_OK && http_code == 200 && !wctx.oom) {
            struct archive* a = archive_read_new();
            if (!a) return false;
            archive_read_support_format_tar(a);
            archive_read_support_filter_all(a);
            bool ok = false;
            if (archive_read_open_memory(a, gs->mem_buf, wctx.len) == ARCHIVE_OK) {
                ok = extract_archive_cifs(gs, a, cb, user_data, out_cif_count, out_skipped);
            }
            archive_read_free(a);
            return ok;
        }

        if (http_code == 401 && use_auth && tm) {
            token_manager_refresh(tm, tok_gen);
            continue;
        }
        if (http_code == 404 || wctx.oom) return false;
        sleep((attempt < 4) ? (1U << attempt) : 16U);
    }
    return false;
}

/* ===================== Large-Tar Streaming Bridge (curl_multi -> libarchive) ===================== */

typedef struct {
    GcsStreamer*       gs;
    struct curl_slist* hdrs;
    size_t             buf_len;
    bool               served_current;
    bool               done;
    bool               oom;
    bool               attached;
    CURLcode           curl_err;
    long               http_code;
} StreamBridge;

static size_t stream_write_cb(char* ptr, size_t size, size_t nmemb, void* userdata) {
    StreamBridge* st = (StreamBridge*)userdata;
    size_t n = size * nmemb;
    if (n == 0) return 0;
    if (!ensure_buf(&st->gs->mem_buf, &st->gs->mem_cap, st->buf_len + n)) {
        st->oom = true;
        return 0;
    }
    memcpy(st->gs->mem_buf + st->buf_len, ptr, n);
    st->buf_len += n;
    return n;
}

static void stream_pump_once(StreamBridge* st) {
    int running = 0;
    CURLMcode mc = curl_multi_perform(st->gs->multi, &running);
    if (mc != CURLM_OK) {
        st->done     = true;
        st->curl_err = (CURLcode)1;
        return;
    }

    int msgs = 0;
    CURLMsg* msg;
    while ((msg = curl_multi_info_read(st->gs->multi, &msgs)) != NULL) {
        if (msg->msg == CURLMSG_DONE && msg->easy_handle == st->gs->stream_easy) {
            st->done     = true;
            st->curl_err = msg->data.result;
            curl_easy_getinfo(st->gs->stream_easy, CURLINFO_RESPONSE_CODE, &st->http_code);
        }
    }

    if (!st->done && running == 0) {
        st->done = true;
        curl_easy_getinfo(st->gs->stream_easy, CURLINFO_RESPONSE_CODE, &st->http_code);
    }
    if (!st->done && st->buf_len == 0) {
        curl_multi_poll(st->gs->multi, NULL, 0, 1000, NULL);
    }
}

static la_ssize_t stream_archive_read_cb(struct archive* a, void* client_data, const void** out_buf) {
    (void)a;
    StreamBridge* st = (StreamBridge*)client_data;

    if (st->served_current) {
        st->buf_len        = 0;
        st->served_current = false;
    }

    while (st->buf_len == 0 && !st->done) {
        stream_pump_once(st);
    }

    if (st->oom) return -1;
    if (st->buf_len > 0) {
        *out_buf           = st->gs->mem_buf;
        st->served_current = true;
        return (la_ssize_t)st->buf_len;
    }

    if (st->curl_err != CURLE_OK || (st->http_code != 0 && st->http_code != 200)) return -1;
    return 0;
}

static int stream_archive_close_cb(struct archive* a, void* client_data) {
    (void)a;
    StreamBridge* st = (StreamBridge*)client_data;
    if (st->attached) {
        curl_multi_remove_handle(st->gs->multi, st->gs->stream_easy);
        st->attached = false;
    }
    if (st->hdrs) {
        curl_slist_free_all(st->hdrs);
        st->hdrs = NULL;
    }
    return ARCHIVE_OK;
}

static bool stream_large_tar_and_extract(GcsStreamer* gs, TokenManager* tm,
                                         const char* url, bool use_auth,
                                         TarCifCallback cb, void* user_data,
                                         uint64_t* out_cif_count, uint64_t* out_skipped) {
    StreamBridge st;
    bool opened = false;

    for (int attempt = 0; attempt < MAX_HTTP_RETRIES; attempt++) {
        memset(&st, 0, sizeof(st));
        st.gs       = gs;
        st.curl_err = CURLE_OK;

        uint64_t tok_gen = 0;
        st.hdrs = make_auth_header(tm, use_auth, &tok_gen);
        configure_common_curl(gs->stream_easy, url, st.hdrs);
        curl_easy_setopt(gs->stream_easy, CURLOPT_WRITEFUNCTION, stream_write_cb);
        curl_easy_setopt(gs->stream_easy, CURLOPT_WRITEDATA, &st);

        if (curl_multi_add_handle(gs->multi, gs->stream_easy) != CURLM_OK) {
            curl_slist_free_all(st.hdrs);
            return false;
        }
        st.attached = true;

        while (st.buf_len == 0 && !st.done) {
            stream_pump_once(&st);
        }
        curl_easy_getinfo(gs->stream_easy, CURLINFO_RESPONSE_CODE, &st.http_code);

        if (!st.oom && st.http_code == 200 && (!st.done || st.curl_err == CURLE_OK)) {
            opened = true;
            break;
        }

        curl_multi_remove_handle(gs->multi, gs->stream_easy);
        st.attached = false;
        curl_slist_free_all(st.hdrs);
        st.hdrs = NULL;

        if (st.http_code == 401 && use_auth && tm) {
            token_manager_refresh(tm, tok_gen);
            continue;
        }
        if (st.http_code == 404 || st.oom) return false;
        sleep((attempt < 4) ? (1U << attempt) : 16U);
    }

    if (!opened) return false;

    struct archive* a = archive_read_new();
    if (!a) {
        stream_archive_close_cb(NULL, &st);
        return false;
    }
    archive_read_support_format_tar(a);
    archive_read_support_filter_all(a);

    bool ok = false;
    if (archive_read_open2(a, &st, NULL, stream_archive_read_cb, NULL, stream_archive_close_cb) == ARCHIVE_OK) {
        ok = extract_archive_cifs(gs, a, cb, user_data, out_cif_count, out_skipped);
    } else {
        stream_archive_close_cb(a, &st);
    }
    archive_read_free(a);

    if (st.curl_err != CURLE_OK || st.oom || (st.http_code != 0 && st.http_code != 200)) ok = false;
    return ok;
}

bool gcs_stream_tar_url(GcsStreamer* gs, TokenManager* tm,
                        const char* url, bool use_auth, uint64_t known_size,
                        TarCifCallback cb, void* user_data,
                        uint64_t* out_cif_count, uint64_t* out_skipped) {
    if (!gs || !url || !cb) return false;
    if (known_size > 0 && known_size <= SMALL_TAR_MEM_LIMIT) {
        return fetch_small_tar_and_extract(gs, tm, url, use_auth, cb, user_data, out_cif_count, out_skipped);
    }
    return stream_large_tar_and_extract(gs, tm, url, use_auth, cb, user_data, out_cif_count, out_skipped);
}

bool gcs_stream_tar_file(GcsStreamer* gs, const char* path,
                         TarCifCallback cb, void* user_data,
                         uint64_t* out_cif_count, uint64_t* out_skipped) {
    if (!gs || !path || !cb) return false;
    struct archive* a = archive_read_new();
    if (!a) return false;
    archive_read_support_format_tar(a);
    archive_read_support_filter_all(a);

    bool ok = false;
    if (archive_read_open_filename(a, path, 65536) == ARCHIVE_OK) {
        ok = extract_archive_cifs(gs, a, cb, user_data, out_cif_count, out_skipped);
    }
    archive_read_free(a);
    return ok;
}
