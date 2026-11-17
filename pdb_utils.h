#ifndef BPDB_PDB_UTILS_H
#define BPDB_PDB_UTILS_H

#define AFDB_STATIC_META_STRING \
    "loop_\n" \
    "_audit_author.name\n" \
    "_audit_author.pdbx_ordinal\n" \
    "\"Jumper, John\"               1  \n" \
    "\"Evans, Richard\"             2  \n" \
    "\"Pritzel, Alexander\"         3  \n" \
    "\"Green, Tim\"                 4  \n" \
    "\"Figurnov, Michael\"          5  \n" \
    "\"Ronneberger, Olaf\"          6  \n" \
    "\"Tunyasuvunakool, Kathryn\"   7  \n" \
    "\"Bates, Russ\"                8  \n" \
    "\"Zidek, Augustin\"            9  \n" \
    "\"Potapenko, Anna\"            10 \n" \
    "\"Bridgland, Alex\"            11 \n" \
    "\"Meyer, Clemens\"             12 \n" \
    "\"Kohl, Simon A. A.\"          13 \n" \
    "\"Ballard, Andrew J.\"         14 \n" \
    "\"Cowie, Andrew\"              15 \n" \
    "\"Romera-Paredes, Bernardino\" 16 \n" \
    "\"Nikolov, Stanislav\"         17 \n" \
    "\"Jain, Rishub\"               18 \n" \
    "\"Adler, Jonas\"               19 \n" \
    "\"Back, Trevor\"               20 \n" \
    "\"Petersen, Stig\"             21 \n" \
    "\"Reiman, David\"              22 \n" \
    "\"Clancy, Ellen\"              23 \n" \
    "\"Zielinski, Michal\"          24 \n" \
    "\"Steinegger, Martin\"         25 \n" \
    "\"Pacholska, Michalina\"       26 \n" \
    "\"Berghammer, Tamas\"          27 \n" \
    "\"Silver, David\"              28 \n" \
    "\"Vinyals, Oriol\"             29 \n" \
    "\"Senior, Andrew W.\"          30 \n" \
    "\"Kavukcuoglu, Koray\"         31 \n" \
    "\"Kohli, Pushmeet\"            32 \n" \
    "\"Hassabis, Demis\"            33 \n" \
    "#\n" \
    "_audit_conform.dict_location https://raw.githubusercontent.com/ihmwg/ModelCIF/master/dist/mmcif_ma.dic\n" \
    "_audit_conform.dict_name     mmcif_ma.dic\n" \
    "_audit_conform.dict_version  1.3.9\n" \
    "#\n"

#include "bpdb.h"
#include <stdbool.h>
#include <stdio.h>

uint8_t res3_to_type_id(const char* src, int width);
const char* type_id_to_res3(uint8_t resType);

uint8_t encode_struct_conf_type(const char* type_str);
const char* decode_struct_conf_type(uint8_t type_id);

int get_atom_idx(uint8_t one_letter_code, const char* atom_name, bool is_terminal);

// The missing IO routines:
size_t parse_meta_struct_conf(const char* cif_data, const char* end, MetaStructConf** out_confs);
void write_meta_struct_conf(FILE* out_f, const MetaStructConf* confs, size_t count);

#endif // BPDB_PDB_UTILS_H
