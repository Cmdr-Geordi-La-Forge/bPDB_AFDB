#define _POSIX_C_SOURCE 200809L
#define _FILE_OFFSET_BITS 64

#include "pack_writer.h"
#include "bpdb_internal.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

static const uint8_t ZERO_PAD16[16] = {0};

static void format_part_path(char* out, size_t cap, const char* base_path,
                             uint16_t planned_parts, uint16_t part_idx, bool is_tmp) {
    if (planned_parts <= 1) {
        snprintf(out, cap, "%s.bpdbpack%s", base_path, is_tmp ? ".tmp" : "");
    } else {
        snprintf(out, cap, "%s-p%03u.bpdbpack%s", base_path, (unsigned)part_idx, is_tmp ? ".tmp" : "");
    }
}

bool pack_mkdir_p(const char* dir_path) {
    if (!dir_path || dir_path[0] == '\0') return true;
    char tmp[512];
    size_t len = strlen(dir_path);
    if (len >= sizeof(tmp)) return false;
    memcpy(tmp, dir_path, len + 1);
    for (char* p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            if (mkdir(tmp, 0755) != 0 && errno != EEXIST) {
                fprintf(stderr, "Error: mkdir('%s') failed: %s\n", tmp, strerror(errno));
                return false;
            }
            *p = '/';
        }
    }
    if (mkdir(tmp, 0755) != 0 && errno != EEXIST) {
        fprintf(stderr, "Error: mkdir('%s') failed: %s\n", tmp, strerror(errno));
        return false;
    }
    return true;
}

static bool ensure_parent_dir(const char* file_path) {
    const char* slash = strrchr(file_path, '/');
    if (!slash || slash == file_path) return true;
    char dir[512];
    size_t dlen = (size_t)(slash - file_path);
    if (dlen >= sizeof(dir)) return false;
    memcpy(dir, file_path, dlen);
    dir[dlen] = '\0';
    return pack_mkdir_p(dir);
}

static bool validate_pack_file(const char* path, uint16_t expected_part_idx, BpdbPackHeader* out_hdr) {
    FILE* f = fopen(path, "rb");
    if (!f) return false;
    BpdbPackHeader hdr;
    bool ok = (fread(&hdr, 1, sizeof(hdr), f) == sizeof(hdr));
    fclose(f);
    if (!ok) return false;
    if (memcmp(hdr.magic, BPDB_PACK_MAGIC, 4) != 0 || hdr.version != BPDB_PACK_VERSION) return false;
    if (hdr.part_idx != expected_part_idx || hdr.part_count == 0) return false;
    if (out_hdr) *out_hdr = hdr;
    return true;
}

uint16_t pack_plan_parts(uint64_t est_bpdb_bytes, uint64_t target_pack_bytes) {
    if (target_pack_bytes == 0 || est_bpdb_bytes == 0) return 1;
    uint64_t single_limit = target_pack_bytes + (target_pack_bytes * 35ULL) / 100ULL;
    if (est_bpdb_bytes <= single_limit) return 1;

    uint64_t k = (est_bpdb_bytes + target_pack_bytes / 2ULL) / target_pack_bytes;
    if (k < 2) k = 2;
    if (k > BPDB_PACK_MAX_PARTS) k = BPDB_PACK_MAX_PARTS;
    return (uint16_t)k;
}

bool pack_already_complete(const char* base_path, uint16_t planned_parts) {
    char p0[512];
    format_part_path(p0, sizeof(p0), base_path, planned_parts, 0, false);

    BpdbPackHeader hdr0;
    if (!validate_pack_file(p0, 0, &hdr0)) return false;
    if (planned_parts <= 1) return hdr0.part_count == 1;

    uint16_t actual_k = hdr0.part_count;
    for (uint16_t p = 1; p < actual_k; p++) {
        char pn[512];
        format_part_path(pn, sizeof(pn), base_path, planned_parts, p, false);
        BpdbPackHeader hdrn;
        if (!validate_pack_file(pn, p, &hdrn) || hdrn.part_count != actual_k) {
            return false;
        }
    }
    return true;
}

static bool open_part_if_needed(BpdbPackWriter* pw, uint16_t p) {
    PackPartState* ps = &pw->parts[p];
    if (ps->fp) return true;

    format_part_path(ps->tmp_path, sizeof(ps->tmp_path), pw->base_path, pw->planned_parts, p, true);
    format_part_path(ps->final_path, sizeof(ps->final_path), pw->base_path, pw->planned_parts, p, false);

    if (!ensure_parent_dir(ps->tmp_path)) return false;
    ps->fp = fopen(ps->tmp_path, "wb");
    if (!ps->fp) {
        fprintf(stderr, "Error: fopen('%s') failed: %s\n", ps->tmp_path, strerror(errno));
        return false;
    }

    setvbuf(ps->fp, NULL, _IOFBF, 262144);

    BpdbPackHeader zhdr;
    memset(&zhdr, 0, sizeof(zhdr));
    if (fwrite(&zhdr, 1, sizeof(zhdr), ps->fp) != sizeof(zhdr)) {
        fclose(ps->fp);
        ps->fp = NULL;
        return false;
    }
    ps->cur_offset    = sizeof(BpdbPackHeader);
    ps->payload_bytes = 0;
    return true;
}

static bool write_record_to_part(BpdbPackWriter* pw, uint16_t p,
                                 const char accession[12], uint16_t frag_num, uint16_t frag_total,
                                 const uint8_t* data, size_t len) {
    if (!open_part_if_needed(pw, p)) return false;
    PackPartState* ps = &pw->parts[p];

    if ((ps->cur_offset >> 4) > 0xFFFFFFFFULL || len > 0xFFFFFFFFULL) return false;

    if (ps->toc_count >= ps->toc_cap) {
        size_t ncap = (ps->toc_cap == 0) ? 1024U : (ps->toc_cap * 2U);
        BpdbPackTocEntry* nt = realloc(ps->toc, ncap * sizeof(BpdbPackTocEntry));
        if (!nt) return false;
        ps->toc     = nt;
        ps->toc_cap = ncap;
    }

    size_t pad = (16U - (len & 15U)) & 15U;
    if (fwrite(data, 1, len, ps->fp) != len) return false;
    if (pad > 0 && fwrite(ZERO_PAD16, 1, pad, ps->fp) != pad) return false;

    BpdbPackTocEntry* e = &ps->toc[ps->toc_count++];
    memcpy(e->accession, accession, 12);
    e->frag_num     = frag_num;
    e->frag_total   = frag_total;
    e->offset_div16 = (uint32_t)(ps->cur_offset >> 4);
    e->byte_len     = (uint32_t)len;

    ps->cur_offset    += (uint64_t)(len + pad);
    ps->payload_bytes += (uint64_t)len;
    pw->total_records++;
    pw->total_payload_bytes += (uint64_t)len;
    return true;
}

bool pack_writer_init(BpdbPackWriter* pw, const char* base_path, uint32_t tax_id,
                      uint16_t pack_flags, uint64_t est_bpdb_bytes, uint64_t target_pack_bytes) {
    memset(pw, 0, sizeof(*pw));
    if (!base_path || strlen(base_path) >= sizeof(pw->base_path)) return false;
    snprintf(pw->base_path, sizeof(pw->base_path), "%s", base_path);

    pw->tax_id        = tax_id;
    pw->planned_parts = pack_plan_parts(est_bpdb_bytes, target_pack_bytes);
    pw->pack_flags    = pack_flags | ((pw->planned_parts > 1) ? BPDB_PACK_FLAG_SPLIT_PART : 0);
    pw->part_target   = (pw->planned_parts > 1) ? (est_bpdb_bytes / pw->planned_parts) : 0;
    pw->cur_part      = 0;

    pw->parts = calloc(pw->planned_parts, sizeof(PackPartState));
    if (!pw->parts) return false;

    if (!open_part_if_needed(pw, 0)) {
        free(pw->parts);
        pw->parts = NULL;
        return false;
    }
    pthread_mutex_init(&pw->mu, NULL);
    return true;
}

bool pack_writer_add_record(BpdbPackWriter* pw, const uint8_t* bpdb_data, size_t bpdb_len) {
    if (!pw || !bpdb_data || bpdb_len < sizeof(BpdbV2Header)) return false;

    BpdbV2Header hdr;
    memcpy(&hdr, bpdb_data, sizeof(hdr));
    if (memcmp(hdr.magic, "bPDB", 4) != 0 || hdr.version != 2) return false;

    bool is_frag_candidate = (hdr.frag_num > 1) || (hdr.total_res == BPDB_FRAG_WINDOW_RES);

    pthread_mutex_lock(&pw->mu);
    if (pw->failed) {
        pthread_mutex_unlock(&pw->mu);
        return false;
    }

    bool ok = true;
    if (is_frag_candidate) {
        if (pw->frag_count >= pw->frag_cap) {
            size_t ncap = (pw->frag_cap == 0) ? 256U : (pw->frag_cap * 2U);
            HeldFragRecord* nf = realloc(pw->frags, ncap * sizeof(HeldFragRecord));
            if (!nf) {
                pw->failed = true;
                pthread_mutex_unlock(&pw->mu);
                return false;
            }
            pw->frags    = nf;
            pw->frag_cap = ncap;
        }
        uint8_t* copy = malloc(bpdb_len);
        if (!copy) {
            pw->failed = true;
            pthread_mutex_unlock(&pw->mu);
            return false;
        }
        memcpy(copy, bpdb_data, bpdb_len);
        HeldFragRecord* fr = &pw->frags[pw->frag_count++];
        memcpy(fr->accession, hdr.accession, 12);
        fr->frag_num = hdr.frag_num;
        fr->len      = (uint32_t)bpdb_len;
        fr->data     = copy;
        pw->frag_bytes += (uint64_t)bpdb_len;
    } else {
        // Standalone -F1: safe to roll over strictly on this fresh -F1 boundary
        if (pw->planned_parts > 1 &&
            pw->cur_part + 1 < pw->planned_parts &&
            pw->parts[pw->cur_part].cur_offset >= pw->part_target) {
            pw->cur_part++;
        }
        ok = write_record_to_part(pw, pw->cur_part, hdr.accession, hdr.frag_num, 1, bpdb_data, bpdb_len);
        if (!ok) pw->failed = true;
    }

    pthread_mutex_unlock(&pw->mu);
    return ok;
}

static int cmp_held_frags(const void* a, const void* b) {
    const HeldFragRecord* fa = (const HeldFragRecord*)a;
    const HeldFragRecord* fb = (const HeldFragRecord*)b;
    return bpdb_acc_frag_cmp(fa->accession, fa->frag_num, fb->accession, fb->frag_num);
}

static int cmp_toc_entries(const void* a, const void* b) {
    const BpdbPackTocEntry* ea = (const BpdbPackTocEntry*)a;
    const BpdbPackTocEntry* eb = (const BpdbPackTocEntry*)b;
    return bpdb_acc_frag_cmp(ea->accession, ea->frag_num, eb->accession, eb->frag_num);
}

static void cleanup_writer_resources(BpdbPackWriter* pw, bool remove_tmps) {
    for (size_t i = 0; i < pw->frag_count; i++) free(pw->frags[i].data);
    free(pw->frags);
    pw->frags      = NULL;
    pw->frag_count = 0;

    if (pw->parts) {
        for (uint16_t p = 0; p < pw->planned_parts; p++) {
            if (pw->parts[p].fp) {
                fclose(pw->parts[p].fp);
                pw->parts[p].fp = NULL;
            }
            if (remove_tmps && pw->parts[p].tmp_path[0] != '\0') {
                unlink(pw->parts[p].tmp_path);
            }
            free(pw->parts[p].toc);
            pw->parts[p].toc = NULL;
        }
        free(pw->parts);
        pw->parts = NULL;
    }
    pthread_mutex_destroy(&pw->mu);
}

bool pack_writer_finish(BpdbPackWriter* pw) {
    if (!pw) return false;
    if (pw->failed) {
        cleanup_writer_resources(pw, true);
        return false;
    }

    // 1. Sort held fragments and emit every accession's F1..Fn contiguously in one part
    if (pw->frag_count > 0) {
        qsort(pw->frags, pw->frag_count, sizeof(HeldFragRecord), cmp_held_frags);

        size_t i = 0;
        while (i < pw->frag_count) {
            size_t j = i + 1;
            while (j < pw->frag_count &&
                   memcmp(pw->frags[j].accession, pw->frags[i].accession, 12) == 0) {
                j++;
            }

            uint16_t unique_frags = 0;
            for (size_t k = i; k < j; k++) {
                if (k == i || pw->frags[k].frag_num != pw->frags[k - 1].frag_num) {
                    unique_frags++;
                }
            }
            if (unique_frags > 1) pw->multi_frag_proteins++;

            // Roll over ONLY between complete accessions (fresh -F1 boundary)
            if (pw->planned_parts > 1 &&
                pw->cur_part + 1 < pw->planned_parts &&
                pw->parts[pw->cur_part].cur_offset >= pw->part_target) {
                pw->cur_part++;
            }

            for (size_t k = i; k < j; k++) {
                if (k > i && pw->frags[k].frag_num == pw->frags[k - 1].frag_num) continue;
                if (unique_frags > 1) {
                    uint64_t flags = 0;
                    memcpy(&flags, pw->frags[k].data + 8, sizeof(flags));
                    flags |= FLAG_FRAGMENTED;
                    memcpy(pw->frags[k].data + 8, &flags, sizeof(flags));
                }
                if (!write_record_to_part(pw, pw->cur_part,
                                          pw->frags[k].accession, pw->frags[k].frag_num, unique_frags,
                                          pw->frags[k].data, pw->frags[k].len)) {
                    cleanup_writer_resources(pw, true);
                    return false;
                }
            }
            i = j;
        }
    }

    // 2. Trim unused trailing parts (zero empty parts)
    uint16_t actual_parts = 0;
    for (uint16_t p = 0; p < pw->planned_parts; p++) {
        if (pw->parts[p].fp != NULL && pw->parts[p].toc_count > 0) {
            actual_parts = (uint16_t)(p + 1);
        }
    }
    if (actual_parts == 0) {
        if (!open_part_if_needed(pw, 0)) {
            cleanup_writer_resources(pw, true);
            return false;
        }
        actual_parts = 1;
    }

    // 3. Finalize each active part: sort TOC, append TOC, write 64B header at offset 0
    for (uint16_t p = 0; p < actual_parts; p++) {
        PackPartState* ps = &pw->parts[p];
        if (!ps->fp && !open_part_if_needed(pw, p)) {
            cleanup_writer_resources(pw, true);
            return false;
        }
        if (ps->toc_count > 1) {
            qsort(ps->toc, ps->toc_count, sizeof(BpdbPackTocEntry), cmp_toc_entries);
        }

        uint64_t toc_off = ps->cur_offset;
        if (ps->toc_count > 0 &&
            fwrite(ps->toc, sizeof(BpdbPackTocEntry), ps->toc_count, ps->fp) != ps->toc_count) {
            cleanup_writer_resources(pw, true);
            return false;
        }

        BpdbPackHeader hdr;
        memset(&hdr, 0, sizeof(hdr));
        memcpy(hdr.magic, BPDB_PACK_MAGIC, 4);
        hdr.version       = BPDB_PACK_VERSION;
        hdr.flags         = pw->pack_flags;
        hdr.tax_id        = pw->tax_id;
        hdr.part_idx      = p;
        hdr.part_count    = actual_parts;
        hdr.record_count  = (uint32_t)ps->toc_count;
        hdr.toc_offset    = toc_off;
        hdr.payload_bytes = ps->payload_bytes;

        if (fseeko(ps->fp, 0, SEEK_SET) != 0 ||
            fwrite(&hdr, 1, sizeof(hdr), ps->fp) != sizeof(hdr) ||
            fflush(ps->fp) != 0) {
            cleanup_writer_resources(pw, true);
            return false;
        }
        fclose(ps->fp);
        ps->fp = NULL;
    }

    // 4. Atomic commit: rename parts 1..actual_parts-1 first, then part 0 LAST
    for (uint16_t p = 1; p < actual_parts; p++) {
        if (rename(pw->parts[p].tmp_path, pw->parts[p].final_path) != 0) {
            cleanup_writer_resources(pw, true);
            return false;
        }
    }
    if (rename(pw->parts[0].tmp_path, pw->parts[0].final_path) != 0) {
        cleanup_writer_resources(pw, true);
        return false;
    }

    cleanup_writer_resources(pw, false);
    return true;
}

void pack_writer_abort(BpdbPackWriter* pw) {
    if (!pw) return;
    cleanup_writer_resources(pw, true);
}
