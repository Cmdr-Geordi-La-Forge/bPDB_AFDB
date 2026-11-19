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
#define FLAG_RESOLUTION (1ULL << 20)
#define FLAG_ENTITY     (1ULL << 34)
#define FLAG_STRUCT     (1ULL << 48)
#define FLAG_TEMPLATES  (1ULL << 49)
#define FLAG_DBREF      (1ULL << 50)
#define FLAG_AUDIT      (1ULL << 51)

// --- Template TLV Tags (used inside FLAG_TEMPLATES) ---
enum {
    TPL_TAG_END       = 0x00, // 1-byte record terminator
    TPL_TAG_ID        = 0x01, // uint16_t (len = 2)
    TPL_TAG_ACCESSION = 0x02, // ASCII bytes (len = 1..255)
    TPL_TAG_DB_NAME   = 0x03, // ASCII bytes (len = 1..255)
    TPL_TAG_ASYM_ID   = 0x04  // ASCII bytes (len = 1..255, e.g. "A", "Ag", "AAA")
};

// --- DBRef TLV Tags (used inside FLAG_DBREF) ---
enum {
    DBREF_TAG_END          = 0x00, // 1-byte record terminator
    DBREF_TAG_DB_NAME      = 0x01, // ASCII string (e.g. "UNP")
    DBREF_TAG_ACCESSION    = 0x02, // ASCII string (e.g. "A0A024RCN7")
    DBREF_TAG_ID_CODE      = 0x03, // ASCII string (e.g. "A0A024RCN7_HUMAN")
    DBREF_TAG_DB_ALIGN     = 0x04, // int32_t[2]: db_seq_begin, db_seq_end
    DBREF_TAG_SEQ_ALIGN    = 0x05, // int32_t[2]: seq_begin, seq_end (omitted when 0,0)
    DBREF_TAG_CHAIN        = 0x06, // char (omitted when default 'A')
    DBREF_TAG_ICODES       = 0x07, // char[4]: iCode_beg/end, db_iCode_beg/end (omitted when "????")
    DBREF_TAG_CHECKSUM     = 0x08, // uint64_t
    DBREF_TAG_TAX_ID       = 0x09, // uint32_t
    DBREF_TAG_VERSION_DATE = 0x0A, // ASCII string (e.g. "2014-05-14")
    DBREF_TAG_GENE         = 0x0B, // ASCII string (omitted if absent)
    DBREF_TAG_ORGANISM     = 0x0C, // ASCII string (omitted if absent)
    DBREF_TAG_ISOFORM      = 0x0D  // ASCII string (omitted if "?")
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
    char     accession[13];
    int      frag_num;
    uint64_t flags;

    float global_plddt;
    char* entity_desc;

    MetaStructConf* confs;
    size_t          conf_count;

    TemplateEntry* templates;
    size_t         template_count;

    DbRefEntry dbref;

    AuditRevisionEntry* revisions;
    uint32_t            revision_count;
} AfdbMetadata;

typedef struct {
    uint8_t code;
    int32_t coords[30][4];
    uint8_t present[30];
    uint8_t xyz_lsb[30]; // bit 0=X, bit 1=Y, bit 2=Z (AF2: -0.000; AF3: ion charge)
} ParsedResidue;

typedef struct {
    ParsedResidue* residues;
    int            res_count;
    int            first_seq_id;
    char           chain_id;
    uint8_t        fmt_widths[8];
} ParsedChain;

typedef struct {
    uint32_t       total_res;
    char           chain_id;
    int            first_seq_id;
    uint8_t        fmt_widths[8];
    const uint8_t* seq;
    const uint8_t* coord_stream;
    const uint8_t* stream_end;
} ChainView;

void afdb_meta_free(AfdbMetadata* meta);
void parsed_chain_free(ParsedChain* chain);

#endif // BPDB_INTERNAL_H
