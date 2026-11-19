// src/cif/afdb_strings.h
#ifndef BPDB_AFDB_STRINGS_H
#define BPDB_AFDB_STRINGS_H

#include <stddef.h>

typedef struct {
    char        code;
    const char* formula;
    const char* weight;
    const char* id;
    const char* name;
    const char* type;
} ChemCompRow;

extern const ChemCompRow CHEM_COMP_ROWS[20];
extern const size_t      NUM_CHEM_COMP_ROWS;

// --- printf Format Macros (enables compile-time -Wformat checking) ---
#define AFDB_ENTRY_ID_FORMAT \
    "data_AF-%s-F%d\n" \
    "#\n" \
    "_entry.id AF-%s-F%d\n" \
    "#\n"

#define AFDB_DATABASE_2_FORMAT \
    "_database_2.database_code AF-%s-F%d\n" \
    "_database_2.database_id   AlphaFoldDB\n" \
    "#\n"

#define AFDB_QA_METRIC_GLOBAL_FORMAT \
    "_ma_qa_metric_global.metric_id    1\n" \
    "_ma_qa_metric_global.metric_value %.2f\n" \
    "_ma_qa_metric_global.model_id     1\n" \
    "_ma_qa_metric_global.ordinal_id   1\n" \
    "#\n"

#define AFDB_DATABASE_STATUS_FORMAT \
    "_pdbx_database_status.entry_id                      AF-%s-F%d\n" \
    "_pdbx_database_status.recvd_initial_deposition_date %s\n" \
    "_pdbx_database_status.status_code                   REL\n" \
    "#\n"

#define AFDB_STRUCT_ASYM_FORMAT \
    "_struct_asym.entity_id 1\n" \
    "_struct_asym.id        %c\n" \
    "#\n"

#define AFDB_STRUCT_REF_FORMAT \
    "_struct_ref.db_code                  %s\n" \
    "_struct_ref.db_name                  %s\n" \
    "_struct_ref.entity_id                1\n" \
    "_struct_ref.id                       1\n" \
    "_struct_ref.pdbx_align_begin         %d\n" \
    "_struct_ref.pdbx_align_end           %d\n" \
    "_struct_ref.pdbx_db_accession        %s\n" \
    "_struct_ref.pdbx_db_isoform          %s\n" \
    "_struct_ref.pdbx_seq_one_letter_code "

#define AFDB_STRUCT_REF_SEQ_FORMAT \
    "_struct_ref_seq.align_id                    1\n" \
    "_struct_ref_seq.db_align_beg                %d\n" \
    "_struct_ref_seq.db_align_end                %d\n" \
    "_struct_ref_seq.pdbx_PDB_id_code            AF-%s-F%d\n" \
    "_struct_ref_seq.pdbx_auth_seq_align_beg     %d\n" \
    "_struct_ref_seq.pdbx_auth_seq_align_end     %d\n" \
    "_struct_ref_seq.pdbx_db_accession           %s\n" \
    "_struct_ref_seq.pdbx_db_align_beg_ins_code  ?\n" \
    "_struct_ref_seq.pdbx_db_align_end_ins_code  ?\n" \
    "_struct_ref_seq.pdbx_seq_align_beg_ins_code ?\n" \
    "_struct_ref_seq.pdbx_seq_align_end_ins_code ?\n" \
    "_struct_ref_seq.pdbx_strand_id              %c\n" \
    "_struct_ref_seq.ref_id                      1\n" \
    "_struct_ref_seq.seq_align_beg               %d\n" \
    "_struct_ref_seq.seq_align_end               %d\n" \
    "#\n"

// --- Static mmCIF Headers & Boilerplate Blocks ---
extern const char* AFDB_ATOM_TYPE_HEADER;
extern const char* AFDB_CHEM_COMP_HEADER;
extern const char* AFDB_ENTITY_BLOCK_1;
extern const char* AFDB_ENTITY_BLOCK_2;
extern const char* AFDB_ENTITY_POLY_SEQ_HEADER;
extern const char* AFDB_MA_DATA_MULTI;
extern const char* AFDB_MA_DATA_SINGLE;
extern const char* AFDB_QA_METRIC_LOCAL_HEADER;
extern const char* AFDB_AUDIT_REVISION_DETAILS_HEADER;
extern const char* AFDB_AUDIT_REVISION_HISTORY_HEADER;
extern const char* AFDB_POLY_SEQ_SCHEME_HEADER;

extern const char* AFDB_AUDIT_BLOCK;
extern const char* AFDB_CITATION_BLOCK;
extern const char* AFDB_MA_MODEL_LIST;
extern const char* AFDB_MA_POST_BLOCK;
extern const char* AFDB_MA_TEMPLATE_DETAILS_HEADER;
extern const char* AFDB_MA_TEMPLATE_REF_DB_DETAILS_HEADER;
extern const char* AFDB_MA_TEMPLATE_TRANS_MATRIX;
extern const char* AFDB_PDBX_DATA_USAGE;
extern const char* AFDB_SOFTWARE_ASYM;

extern const char* AFDB_STRUCT_CONF_HEADER;
extern const char* AFDB_STRUCT_CONF_TYPE_HEADER;
extern const char* AFDB_ATOM_SITE_HEADER;

#endif // BPDB_AFDB_STRINGS_H
