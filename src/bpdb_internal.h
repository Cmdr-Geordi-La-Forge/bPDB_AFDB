#ifndef BPDB_INTERNAL_H
#define BPDB_INTERNAL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

// --- Structural Conformation Bit Masks ---
#define CONF_CAT_MASK  0xE0
#define CONF_SUB_MASK  0x1F

#define TYPE_OTHER     0x00
#define CAT_HELIX      0x00
#define CAT_STRAND     0x20
#define CAT_TURN       0x40
#define CAT_BEND       0x60

// --- Binary Metadata Section Flags ---
#define FLAG_AF_MODE    (1ULL << 0)
#define FLAG_AFDB       (1ULL << 1)
#define FLAG_HYDROGENS  (1ULL << 2)
#define FLAG_CIF        (1ULL << 3)
#define FLAG_AF2        (1ULL << 4)
#define FLAG_FRAGMENTED (1ULL << 5)
#define FLAG_RESOLUTION (1ULL << 20)
#define FLAG_ENTITY     (1ULL << 34)
#define FLAG_STRUCT     (1ULL << 48)
#define FLAG_TEMPLATES  (1ULL << 49)
#define FLAG_DBREF      (1ULL << 50)
#define FLAG_AUDIT      (1ULL << 51)
#define FLAG_NEG_ZERO   (1ULL << 52)

// --- AFDB Audit Revision Preset Bitmask (used inside FLAG_AUDIT) ---
#define AUDIT_MASK_V1   (1U << 0) // 2021-07-01 v1.0 Initial release
#define AUDIT_MASK_V2   (1U << 1) // 2021-12-09 v2.0 Remediation (SwissProt)
#define AUDIT_MASK_V3   (1U << 2) // 2022-06-01 v3.0 Remediation (UniProt)
#define AUDIT_MASK_V4   (1U << 3) // 2022-09-30 v4.0 Remediation (accuracy)
#define AUDIT_MASK_V5   (1U << 4) // 2025-03-31 v5.0 Release (UniProt 2024_06)
#define AUDIT_MASK_V6   (1U << 5) // 2025-08-01 v6.0 Release (UniProt 2025_03)
#define AUDIT_MASK_RAW  0xFFU     // Non-standard revisions follow inline

// --- Template TLV Tags (used inside FLAG_TEMPLATES) ---
enum {
    TPL_TAG_END       = 0x00, // 1-byte record terminator
    TPL_TAG_ID        = 0x01, // uint16_t (omitted when template_id == idx + 1)
    TPL_TAG_ACCESSION = 0x02, // ASCII bytes (len = 1..255)
    TPL_TAG_DB_NAME   = 0x03, // ASCII bytes (omitted when "PDB")
    TPL_TAG_ASYM_ID   = 0x04  // ASCII bytes (omitted when "A")
};

// --- DBRef TLV Tags (used inside FLAG_DBREF) ---
enum {
    DBREF_TAG_END          = 0x00, // 1-byte record terminator
    DBREF_TAG_DB_NAME      = 0x01, // ASCII string (omitted when "UNP")
    DBREF_TAG_ACCESSION    = 0x02, // ASCII string (omitted when == hdr.accession)
    DBREF_TAG_ID_CODE      = 0x03, // ASCII string (e.g. "L2HDH_ACCRE")
    DBREF_TAG_DB_ALIGN     = 0x04, // int32_t[2]: db_seq_begin, db_seq_end (omitted when == chain bounds)
    DBREF_TAG_SEQ_ALIGN    = 0x05, // int32_t[2]: seq_begin, seq_end (omitted when 0,0)
    DBREF_TAG_CHAIN        = 0x06, // char (omitted when default 'A')
    DBREF_TAG_ICODES       = 0x07, // char[4]: iCode_beg/end, db_iCode_beg/end (omitted when "????")
    DBREF_TAG_CHECKSUM     = 0x08, // uint64_t (omitted when 0)
    DBREF_TAG_TAX_ID       = 0x09, // uint32_t (omitted when 0)
    DBREF_TAG_VERSION_DATE = 0x0A, // 3-byte packed {yr-1973, mo, day} (or raw ASCII if len != 3)
    DBREF_TAG_GENE         = 0x0B, // ASCII string (omitted if absent)
    DBREF_TAG_ORGANISM     = 0x0C, // ASCII string (omitted if absent)
    DBREF_TAG_ISOFORM      = 0x0D, // ASCII string (omitted if "?")
    DBREF_TAG_ID_SUFFIX    = 0x0E  // ASCII suffix after "<accession>_" (TrEMBL ID_CODE compaction)
};

// --- Wire-Format Structs ---

typedef struct {
    char     chain;
    char     iCode_begin, iCode_end;
    char     db_iCode_begin, db_iCode_end;
    int32_t  seq_begin, seq_end;
    int32_t  db_seq_begin, db_seq_end;
    uint64_t checksum;
    uint32_t ncbi_tax_id;
    char*    db_name;
    char*    db_accession;
    char*    db_id_code;
    char*    version_date;
    char*    gene;
    char*    organism;
    char*    isoform;
} DbRefEntry;

typedef struct {
    uint8_t  resType;
    char     chain[2];
    char     iCode;
    uint32_t resSeq;
} PdbResId;

typedef struct {
    PdbResId start;
    PdbResId end;
    uint16_t length;
    uint8_t  type_id;
    uint8_t  parent_id;
} MetaStructConf;

typedef struct {
    uint8_t  year;
    uint8_t  month;
    uint8_t  day;
    uint8_t  major;
    uint8_t  minor;
    uint8_t  type_len;
    uint16_t desc_len;
} __attribute__((packed)) MetaAuditRevision;

#define BPDB_PARENT_PREV_C 29

// --- Chemistry Dictionary Entry ---
typedef struct {
    uint8_t code;
    char    name[4];
    char    oneLetter;
    uint8_t isNucleic;
    uint8_t numHeavyAtoms;
    uint8_t numHeavyAtomsTerm;
    char    atomNames[25][5];
    char    elements[25][3];
    uint8_t format_group;
    uint8_t numHydrogens;
    uint8_t numHydrogensTerm;
    char    hydroNames[25][5];
    char    hydroElements[25][3];
    uint8_t parent[30];
} ResidueDict;

// --- In-Memory Intermediate Representations ---
typedef struct {
    uint16_t template_id;
    char*    accession;
    char*    db_name;
    char*    asym_id;
} TemplateEntry;

typedef struct {
    MetaAuditRevision meta;
    char* type;
    char* description;
} AuditRevisionEntry;

typedef struct {
    char     magic[4];      //  0..3 : "bPDB"
    uint8_t  version;       //  4    : 2
    uint8_t  chain_id;      //  5    : 'A'
    uint8_t  c_term_oxt;    //  6    : 1 if C-terminal residue has OXT
    uint8_t  cold_codec;    //  7    : 0 = raw Cold Metadata
    uint64_t flags;         //  8..15: FLAG_* bitmask

    uint32_t total_res;     // 16..19
    uint32_t total_atoms;   // 20..23
    int32_t  first_seq_id;  // 24..27
    uint32_t escape_count;  // 28..31: count of 12-byte EscapeAtom entries

    int32_t  origin[3];     // 32..43: atom 0 exact int32 (x, y, z) * 1000
    uint32_t cold_offset;   // 44..47: byte offset where Cold Metadata starts

    char     accession[12]; // 48..59
    uint16_t frag_num;      // 60..61
    uint16_t reserved;      // 62..63
} BpdbV2Header;

typedef char bpdb_v2_hdr_size_check[(sizeof(BpdbV2Header) == 64) ? 1 : -1];

typedef struct {
    int32_t x, y, z;
} EscapeAtom;

typedef struct {
    char     accession[13];
    int      frag_num;
    uint64_t flags;

    uint16_t global_plddt; // scaled x100 (e.g. 93.67 -> 9367)
    char*    entity_desc;

    MetaStructConf* confs;
    size_t          conf_count;
    size_t          conf_cap;

    TemplateEntry* templates;
    size_t         template_count;

    DbRefEntry dbref;

    uint8_t             audit_mask; // 0x01..0x3F = AFDB v1..v6 bitmask, 0xFF = raw revisions[]
    AuditRevisionEntry* revisions;
    uint32_t            revision_count;

    uint32_t* neg_zeros;      // sorted (stream_atom_idx << 2) | axis (1=X, 2=Y, 3=Z)
    uint32_t  neg_zero_count;
    uint32_t  neg_zero_cap;
} AfdbMetadata;

typedef struct {
    uint8_t  code;
    uint8_t  xyz_lsb[25];  // bit 0=X, bit 1=Y, bit 2=Z (AF2: -0.000; AF3: ion charge)
    uint16_t plddt;        // per-residue pLDDT scaled x100
    uint32_t present_mask; // bit a set if heavy atom a is present
    int32_t  xyz[25][3];   // milli-Angstrom coordinates
} ParsedResidue;

typedef struct {
    ParsedResidue* residues;
    int            res_count;
    int            res_cap;
    int            first_seq_id;
    char           chain_id;
    uint8_t        fmt_widths[8];
} ParsedChain;

typedef struct {
    uint32_t          total_res;
    uint32_t          total_atoms;
    uint32_t          escape_count;
    int32_t           first_seq_id;
    int32_t           origin[3];
    char              chain_id;
    uint8_t           c_term_oxt;
    uint8_t           fmt_widths[8];
    const uint8_t*    seq;
    const uint8_t*    plddt_bytes;  // total_res * 2 bytes (uint16_t scaled x100)
    const uint8_t*    coord_stream; // 12-bit bonded-parent deltas (+ >=7B zero pad)
    size_t            coord_len;
    const EscapeAtom* escapes;
} ChainView;

void afdb_meta_free(AfdbMetadata* meta);
void afdb_meta_reset_reuse(AfdbMetadata* meta);
void parsed_chain_free(ParsedChain* chain);

#endif // BPDB_INTERNAL_H
