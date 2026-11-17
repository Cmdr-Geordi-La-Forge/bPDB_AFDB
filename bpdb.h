#ifndef BPDB_H
#define BPDB_H

#include <stdint.h>
#include <stddef.h>
#include <stdio.h> 

// --- Meta-Section Structs & Masks ---
#define CONF_CAT_MASK  0xE0
#define CONF_SUB_MASK  0x1F

#define TYPE_OTHER     0x00
#define CAT_HELIX      0x00
#define CAT_STRAND     0x20
#define CAT_TURN       0x40
#define CAT_BEND       0x60

#define FLAG_AF_MODE    (1ULL << 0)
#define FLAG_AFDB       (1ULL << 1)
#define FLAG_HYDROGENS  (1ULL << 2)
#define FLAG_CIF        (1ULL << 3)
#define FLAG_ELEMENTS   (1ULL << 34) // Placed in Section 3 (> Bit 16)
#define FLAG_STRUCT     (1ULL << 40)

typedef struct {
    uint8_t resType;
    char chain[2];
    char iCode;
    uint32_t resSeq;
} PdbResId;

typedef struct {
    PdbResId start;
    PdbResId end;
    uint16_t length;
    uint8_t type_id;
    uint8_t parent_id;
} MetaStructConf;

typedef struct {
    uint8_t code;
    char name[4];
    char oneLetter;
    uint8_t isNucleic;
    uint8_t numHeavyAtoms;
    uint8_t numHeavyAtomsTerm;
    char atomNames[25][5];
    char elements[25][3];
    uint8_t format_group; 
    uint8_t numHydrogens;
    uint8_t numHydrogensTerm;
    char hydroNames[25][5];
    char hydroElements[25][3];
} ResidueDict;

extern const ResidueDict DICT[];
extern const size_t NUM_DICT_RES;

// --- End-to-end pipeline definitions ---
size_t encode_afdb_to_binary(const char* cif_data, size_t cif_len, uint8_t** out_buf);
void decode_binary_to_afdb(const uint8_t* buf, size_t buf_len, FILE* out_f);
void write_meta_struct_conf(FILE* out_f, const MetaStructConf* confs, size_t count);

#endif // BPDB_H
