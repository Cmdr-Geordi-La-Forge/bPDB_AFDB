#ifndef BPDB_PACK_WRITER_H
#define BPDB_PACK_WRITER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <pthread.h>
#include "pack_format.h"

#define BPDB_FRAG_WINDOW_RES 1400U
#define BPDB_PACK_MAX_PARTS  1024U

typedef struct {
    FILE*             fp;
    char              tmp_path[512];
    char              final_path[512];
    uint64_t          cur_offset;
    uint64_t          payload_bytes;
    BpdbPackTocEntry* toc;
    size_t            toc_count;
    size_t            toc_cap;
} PackPartState;

typedef struct {
    char     accession[12];
    uint16_t frag_num;
    uint32_t len;
    uint8_t* data;
} HeldFragRecord;

typedef struct {
    pthread_mutex_t mu;
    char            base_path[480];
    uint32_t        tax_id;
    uint16_t        pack_flags;
    uint16_t        planned_parts;
    uint16_t        cur_part;
    uint64_t        part_target;
    PackPartState*  parts;

    HeldFragRecord* frags;
    size_t          frag_count;
    size_t          frag_cap;
    uint64_t        frag_bytes;

    uint64_t        total_records;
    uint64_t        total_payload_bytes;
    uint64_t        multi_frag_proteins;
    bool            failed;
} BpdbPackWriter;

uint16_t pack_plan_parts(uint64_t est_bpdb_bytes, uint64_t target_pack_bytes);
bool     pack_already_complete(const char* base_path, uint16_t planned_parts);
bool     pack_mkdir_p(const char* dir_path);

bool pack_writer_init(BpdbPackWriter* pw, const char* base_path, uint32_t tax_id,
                      uint16_t pack_flags, uint64_t est_bpdb_bytes, uint64_t target_pack_bytes);
bool pack_writer_add_record(BpdbPackWriter* pw, const uint8_t* bpdb_data, size_t bpdb_len);
bool pack_writer_finish(BpdbPackWriter* pw);
void pack_writer_abort(BpdbPackWriter* pw);

#endif // BPDB_PACK_WRITER_H
