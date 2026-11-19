#ifndef BPDB_PACK_FORMAT_H
#define BPDB_PACK_FORMAT_H

#include <stdint.h>
#include <string.h>

#define BPDB_PACK_MAGIC   "bPK2"
#define BPDB_PACK_VERSION 2U

#define BPDB_IDX_MAGIC    "bIX2"
#define BPDB_IDX_VERSION  2U

// BpdbPackHeader.flags
#define BPDB_PACK_FLAG_SMALL_CLUSTER (1U << 0) // clustered multi-TaxID pack
#define BPDB_PACK_FLAG_SPLIT_PART    (1U << 1) // part of a balanced multi-part proteome
#define BPDB_PACK_FLAG_AFDB_V6       (1U << 2) // v6 overlay pack (shadows v4 in index)

typedef struct {
    char     magic[4];      //  0..3 : "bPK2"
    uint16_t version;       //  4..5 : 2
    uint16_t flags;         //  6..7 : BPDB_PACK_FLAG_*
    uint32_t tax_id;        //  8..11: NCBI TaxID (0 for multi-TaxID small clusters)
    uint16_t part_idx;      // 12..13: 0-based part index
    uint16_t part_count;    // 14..15: total parts for this TaxID (1 for single pack)
    uint32_t record_count;  // 16..19: number of bPDB records & TOC entries
    uint32_t reserved0;     // 20..23: 0
    uint64_t toc_offset;    // 24..31: 16B-aligned byte offset of BpdbPackTocEntry[record_count]
    uint64_t payload_bytes; // 32..39: sum of unpadded bPDB record byte lengths
    uint8_t  reserved[24];  // 40..63: 0
} BpdbPackHeader;

typedef char bpdb_pack_hdr_size_check[(sizeof(BpdbPackHeader) == 64) ? 1 : -1];

typedef struct {
    char     accession[12]; //  0..11: NUL-padded UniProt accession (sorted ascending)
    uint16_t frag_num;      // 12..13: 1-based fragment number (sorted ascending)
    uint16_t frag_total;    // 14..15: total fragments for this accession (united in this pack)
    uint32_t offset_div16;  // 16..19: byte offset in pack >> 4 (addresses up to 64 GiB)
    uint32_t byte_len;      // 20..23: exact unpadded bPDB v2 byte length
} BpdbPackTocEntry;

typedef char bpdb_pack_toc_size_check[(sizeof(BpdbPackTocEntry) == 24) ? 1 : -1];

typedef struct {
    char     magic[4];           //  0..3 : "bIX2"
    uint16_t version;            //  4..5 : 2
    uint16_t flags;              //  6..7 : 0
    uint64_t entry_count;        //  8..15: total BpdbGlobalIdxEntry count
    uint32_t shard_count;        // 16..19: number of pack paths in shard table
    uint32_t reserved0;          // 20..23: 0
    uint64_t shard_table_offset; // 24..31: byte offset of NUL-terminated relative pack paths
    uint8_t  reserved[32];       // 32..63: 0
} BpdbGlobalIdxHeader;

typedef char bpdb_idx_hdr_size_check[(sizeof(BpdbGlobalIdxHeader) == 64) ? 1 : -1];

typedef struct {
    char     accession[12]; //  0..11: NUL-padded UniProt accession (sorted ascending)
    uint16_t frag_num;      // 12..13: 1-based fragment number (sorted ascending)
    uint16_t frag_total;    // 14..15: total fragments for this accession
    uint32_t shard_id;      // 16..19: index into shard path table
    uint32_t offset_div16;  // 20..23: byte offset in shard >> 4
    uint32_t byte_len;      // 24..27: exact unpadded bPDB v2 byte length
} BpdbGlobalIdxEntry;

typedef char bpdb_idx_entry_size_check[(sizeof(BpdbGlobalIdxEntry) == 28) ? 1 : -1];

static inline int bpdb_acc_frag_cmp(const char a_acc[12], uint16_t a_frag,
                                    const char b_acc[12], uint16_t b_frag) {
    int c = memcmp(a_acc, b_acc, 12);
    if (c != 0) return c;
    if (a_frag < b_frag) return -1;
    if (a_frag > b_frag) return 1;
    return 0;
}

#endif // BPDB_PACK_FORMAT_H
