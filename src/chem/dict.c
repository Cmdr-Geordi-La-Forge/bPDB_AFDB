#include "chem/dict.h"
#include "util/str_util.h"
#include <string.h>

const ResidueDict DICT[] = {
    // 0x01 - 0x1A: Amino Acids (Strict 1-Letter Alphabetical Mapping)
    {0x01, "ALA", 'A', 0,  5,  6, {" N  "," CA "," C  "," O  "," CB "," OXT"}, {" N"," C"," C"," O"," C"," O"}, 3,  8,  9, {" H1 "," H2 "," H3 "," H  "," HA "," HB1"," HB2"," HB3"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x02, "ASX", 'B', 0,  8,  9, {" N  "," CA "," C  "," O  "," CB "," CG "," XD1"," XD2"," OXT"}, {" N"," C"," C"," O"," C"," C"," X"," X"," O"}, 3,  7,  8, {" H1 "," H2 "," H3 "," H  "," HA "," HB1"," HB2"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x03, "CYS", 'C', 0,  6,  7, {" N  "," CA "," C  "," O  "," CB "," SG "," OXT"}, {" N"," C"," C"," O"," C"," S"," O"}, 3,  8,  9, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3"," HG "," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x04, "ASP", 'D', 0,  8,  9, {" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," OD2"," OXT"}, {" N"," C"," C"," O"," C"," C"," O"," O"," O"}, 3,  8,  9, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3"," HD2"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x05, "GLU", 'E', 0,  9, 10, {" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," OE2"," OXT"}, {" N"," C"," C"," O"," C"," C"," C"," O"," O"," O"}, 3, 10, 11, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3"," HG2"," HG3"," HE2"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x06, "PHE", 'F', 0, 11, 12, {" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," OXT"}, {" N"," C"," C"," O"," C"," C"," C"," C"," C"," C"," C"," O"}, 3, 12, 13, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3"," HD1"," HD2"," HE1"," HE2"," HZ "," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x07, "GLY", 'G', 0,  4,  5, {" N  "," CA "," C  "," O  "," OXT"}, {" N"," C"," C"," O"," O"}, 3,  6,  7, {" H1 "," H2 "," H3 "," H  "," HA2"," HA3"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"}},
    {0x08, "HIS", 'H', 0, 10, 11, {" N  "," CA "," C  "," O  "," CB "," CG "," ND1"," CD2"," CE1"," NE2"," OXT"}, {" N"," C"," C"," O"," C"," C"," N"," C"," C"," N"," O"}, 3, 11, 12, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3"," HD1"," HD2"," HE1"," HE2"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x09, "ILE", 'I', 0,  8,  9, {" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"," CD1"," OXT"}, {" N"," C"," C"," O"," C"," C"," C"," C"," O"}, 3, 14, 15, {" H1 "," H2 "," H3 "," H  "," HA "," HB ","HG12","HG13","HG21","HG22","HG23","HD11","HD12","HD13"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x0B, "LYS", 'K', 0,  9, 10, {" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ "," OXT"}, {" N"," C"," C"," O"," C"," C"," C"," C"," N"," O"}, 3, 16, 17, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3"," HG2"," HG3"," HD2"," HD3"," HE2"," HE3"," HZ1"," HZ2"," HZ3"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x0C, "LEU", 'L', 0,  8,  9, {" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," OXT"}, {" N"," C"," C"," O"," C"," C"," C"," C"," O"}, 3, 14, 15, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3"," HG ","HD11","HD12","HD13","HD21","HD22","HD23"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x0D, "MET", 'M', 0,  8,  9, {" N  "," CA "," C  "," O  "," CB "," CG "," SD "," CE "," OXT"}, {" N"," C"," C"," O"," C"," C"," S"," C"," O"}, 3, 12, 13, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3"," HG2"," HG3"," HE1"," HE2"," HE3"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x0E, "ASN", 'N', 0,  8,  9, {" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," ND2"," OXT"}, {" N"," C"," C"," O"," C"," C"," O"," N"," O"}, 3,  9, 10, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3","HD21","HD22"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x0F, "PYL", 'O', 0, 17, 18, {" CB2"," CG2"," CD2"," CE2"," N2 "," CA2"," C2 "," O2 "," NZ "," CE "," CD "," CG "," CB "," CA "," C  "," O  "," N  "," OXT"}, {" C"," C"," C"," C"," N"," C"," C"," O"," N"," C"," C"," C"," C"," C"," C"," O"," N"," O"}, 3, 22, 23, {" H1 "," H2 "," H3 "," HA "," H  "," HB3"," HB2"," HG3"," HG2"," HD3"," HD2"," HE3"," HE2"," HZ "," HA2","HE22","HD32","HD22","HG22","HB12","HB22","HB32"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x10, "PRO", 'P', 0,  7,  8, {" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OXT"}, {" N"," C"," C"," O"," C"," C"," C"," O"}, 3, 11, 12, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3"," HG2"," HG3"," HD2"," HD3"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x11, "GLN", 'Q', 0,  9, 10, {" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," NE2"," OXT"}, {" N"," C"," C"," O"," C"," C"," C"," O"," N"," O"}, 3, 11, 12, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3"," HG2"," HG3","HE21","HE22"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x12, "ARG", 'R', 0, 11, 12, {" N  "," CA "," C  "," O  "," CB "," CG "," CD "," NE "," CZ "," NH1"," NH2"," OXT"}, {" N"," C"," C"," O"," C"," C"," C"," N"," C"," N"," N"," O"}, 3, 16, 17, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3"," HG2"," HG3"," HD2"," HD3"," HE ","HH11","HH12","HH21","HH22"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x13, "SER", 'S', 0,  6,  7, {" N  "," CA "," C  "," O  "," CB "," OG "," OXT"}, {" N"," C"," C"," O"," C"," O"," O"}, 3,  8,  9, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3"," HG "," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x14, "THR", 'T', 0,  7,  8, {" N  "," CA "," C  "," O  "," CB "," OG1"," CG2"," OXT"}, {" N"," C"," C"," O"," C"," O"," C"," O"}, 3, 10, 11, {" H1 "," H2 "," H3 "," H  "," HA "," HB "," HG1","HG21","HG22","HG23"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x15, "SEC", 'U', 0,  6,  7, {" N  "," CA "," CB ","SE  "," C  "," O  "," OXT"}, {" N"," C"," C","SE"," C"," O"," O"}, 3,  8,  9, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3"," HE "," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x16, "VAL", 'V', 0,  7,  8, {" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"," OXT"}, {" N"," C"," C"," O"," C"," C"," C"," O"}, 3, 12, 13, {" H1 "," H2 "," H3 "," H  "," HA "," HB ","HG11","HG12","HG13","HG21","HG22","HG23"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x17, "TRP", 'W', 0, 14, 15, {" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," NE1"," CE2"," CE3"," CZ2"," CZ3"," CH2"," OXT"}, {" N"," C"," C"," O"," C"," C"," C"," C"," N"," C"," C"," C"," C"," C"," O"}, 3, 13, 14, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3"," HD1"," HE1"," HE3"," HZ2"," HZ3"," HH2"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x18, "UNK", 'X', 0,  6,  7, {" N  "," CA "," C  "," O  "," CB "," CG "," OXT"}, {" N"," C"," C"," O"," C"," C"," O"}, 3, 10, 11, {" H1 "," H2 "," H3 "," H  "," HA "," HB1"," HB2"," HG1"," HG2"," HG3"," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x19, "TYR", 'Y', 0, 12, 13, {" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," OH "," OXT"}, {" N"," C"," C"," O"," C"," C"," C"," C"," C"," C"," C"," O"," O"}, 3, 12, 13, {" H1 "," H2 "," H3 "," H  "," HA "," HB2"," HB3"," HD1"," HD2"," HE1"," HE2"," HH "," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x1A, "GLX", 'Z', 0,  9, 10, {" N  "," CA "," C  "," O  "," CB "," CG "," CD "," XE1"," XE2"," OXT"}, {" N"," C"," C"," O"," C"," C"," C"," X"," X"," O"}, 3,  9, 10, {" H1 "," H2 "," H3 "," HA "," HB1"," HB2"," HG1"," HG2"," H  "," HXT"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},

    // 0x21 - 0x2C: Nucleic Acids
    {0x21, "  A", 'A', 1, 23, 24, {" OP3"," P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," N6 "," N1 "," C2 "," N3 "," C4 "," O3T"}, {" O"," P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," O"," C"," N"," C"," N"," C"," C"," N"," N"," C"," N"," C"," O"}, 3, 14, 15, {"HOP3","HOP2","HO5'"," H5'","H5''"," H4'"," H3'"," H2'","HO2'"," H1'"," H8 "," H61"," H62"," H2 ","HO3'"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x22, "  C", 'C', 1, 21, 22, {" OP3"," P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N1 "," C2 "," O2 "," N3 "," C4 "," N4 "," C5 "," C6 "," O3T"}, {" O"," P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," O"," C"," N"," C"," O"," N"," C"," N"," C"," C"," O"}, 3, 14, 15, {"HOP3","HOP2","HO5'"," H5'","H5''"," H4'"," H3'"," H2'","HO2'"," H1'"," H41"," H42"," H5 "," H6 ","HO3'"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x23, "  G", 'G', 1, 24, 25, {" OP3"," P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," O6 "," N1 "," C2 "," N2 "," N3 "," C4 "," O3T"}, {" O"," P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," O"," C"," N"," C"," N"," C"," C"," O"," N"," C"," N"," N"," C"," O"}, 3, 14, 15, {"HOP3","HOP2","HO5'"," H5'","H5''"," H4'"," H3'"," H2'","HO2'"," H1'"," H8 "," H1 "," H21"," H22","HO3'"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x24, "  U", 'U', 1, 21, 22, {" OP3"," P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N1 "," C2 "," O2 "," N3 "," C4 "," O4 "," C5 "," C6 "," O3T"}, {" O"," P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," O"," C"," N"," C"," O"," N"," C"," O"," C"," C"," O"}, 3, 13, 14, {"HOP3","HOP2","HO5'"," H5'","H5''"," H4'"," H3'"," H2'","HO2'"," H1'"," H3 "," H5 "," H6 ","HO3'"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x25, "  I", 'I', 1, 23, 24, {" OP3"," P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," O6 "," N1 "," C2 "," N3 "," C4 "," O3T"}, {" O"," P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," O"," C"," N"," C"," N"," C"," C"," O"," N"," C"," N"," C"," O"}, 3, 13, 14, {"HOP3","HOP2","HO5'"," H5'","H5''"," H4'"," H3'"," H2'","HO2'"," H1'"," H8 "," H1 "," H2 ","HO3'"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x26, " DA", 'A', 1, 22, 23, {" OP3"," P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," N6 "," N1 "," C2 "," N3 "," C4 "," O3T"}, {" O"," P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," C"," N"," C"," N"," C"," C"," N"," N"," C"," N"," C"," O"}, 3, 14, 15, {"HOP3","HOP2","HO5'"," H5'","H5''"," H4'"," H3'"," H2'","H2''"," H1'"," H8 "," H61"," H62"," H2 ","HO3'"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x27, " DC", 'C', 1, 20, 21, {" OP3"," P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N1 "," C2 "," O2 "," N3 "," C4 "," N4 "," C5 "," C6 "," O3T"}, {" O"," P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," C"," N"," C"," O"," N"," C"," N"," C"," C"," O"}, 3, 14, 15, {"HOP3","HOP2","HO5'"," H5'","H5''"," H4'"," H3'"," H2'","H2''"," H1'"," H41"," H42"," H5 "," H6 ","HO3'"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x28, " DG", 'G', 1, 23, 24, {" OP3"," P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," O6 "," N1 "," C2 "," N2 "," N3 "," C4 "," O3T"}, {" O"," P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," C"," N"," C"," N"," C"," C"," O"," N"," C"," N"," N"," C"," O"}, 3, 14, 15, {"HOP3","HOP2","HO5'"," H5'","H5''"," H4'"," H3'"," H2'","H2''"," H1'"," H8 "," H1 "," H21"," H22","HO3'"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x29, " DT", 'T', 1, 21, 22, {" OP3"," P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N1 "," C2 "," O2 "," N3 "," C4 "," O4 "," C5 "," C7 "," C6 "," O3T"}, {" O"," P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," C"," N"," C"," O"," N"," C"," O"," C"," C"," C"," O"}, 3, 15, 16, {"HOP3","HOP2","HO5'"," H5'","H5''"," H4'"," H3'"," H2'","H2''"," H1'"," H3 "," H71"," H72"," H73"," H6 ","HO3'"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x2A, " DU", 'U', 1, 20, 21, {" OP3"," P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N1 "," C2 "," O2 "," N3 "," C4 "," O4 "," C5 "," C6 "," O3T"}, {" O"," P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," C"," N"," C"," O"," N"," C"," O"," C"," C"," O"}, 3, 13, 14, {"HOP3","HOP2","HO5'"," H5'","H5''"," H4'"," H3'"," H2'","H2''"," H1'"," H3 "," H5 "," H6 ","HO3'"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x2B, " DI", 'I', 1, 22, 23, {" OP3"," P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"," N9 "," C8 "," N7 "," C5 "," C6 "," O6 "," N1 "," C2 "," N3 "," C4 "," O3T"}, {" O"," P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," C"," N"," C"," N"," C"," C"," O"," N"," C"," N"," C"," O"}, 3, 13, 14, {"HOP2","HOP3","HO5'"," H5'","H5''"," H4'"," H3'"," H2'","H2''"," H1'"," H8 "," H1 "," H2 ","HO3'"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},
    {0x2C, "  N", 'N', 1, 13, 14, {" OP3"," P  "," OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'"," O3T"}, {" O"," P"," O"," O"," O"," C"," C"," O"," C"," O"," C"," O"," C"," O"}, 3, 11, 12, {"HOP3","HOP2","HO5'"," H5'","H5''"," H4'"," H3'"," H2'","HO2'"," H1'","H1'2","HO3'"}, {" H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"," H"}},

    // 0x61: Solvents
    {0x61, "HOH", 'O', 0,  1,  1, {" O  "}, {" O"}, 0,  2,  2, {" H1 "," H2 "}, {" H"," H"}},
};

const size_t NUM_DICT_RES = sizeof(DICT) / sizeof(DICT[0]);

const char* const ATOM37_ORDER[37] = {
    "N", "CA", "C", "CB", "O", "CG", "CG1", "CG2", "OG", "OG1", "SG", "CD",
    "CD1", "CD2", "ND1", "ND2", "OD1", "OD2", "SD", "CE", "CE1", "CE2", "CE3",
    "NE", "NE1", "NE2", "OE1", "OE2", "CH2", "NH1", "NH2", "OH", "CZ", "CZ2",
    "CZ3", "NZ", "OXT"
};

// Direct O(1) index into the amino acid prefix of DICT[] ('A'..'Z', 'J' unused)
static const int8_t AA_DICT_INDEX[26] = {
     0,  1,  2,  3,  4,  5,  6,  7,  8, -1,
     9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
    19, 20, 21, 22, 23, 24
};

const ResidueDict* dict_find_aa(uint8_t one_letter) {
    if (one_letter < 'A' || one_letter > 'Z') return NULL;
    int8_t idx = AA_DICT_INDEX[one_letter - 'A'];
    return (idx >= 0) ? &DICT[idx] : NULL;
}

int dict_get_atom_idx(uint8_t one_letter, const char* atom_name, bool is_terminal) {
    const ResidueDict* rd = dict_find_aa(one_letter);
    if (!rd) return -1;

    int limit = is_terminal ? rd->numHeavyAtomsTerm : rd->numHeavyAtoms;
    for (int a = 0; a < limit; a++) {
        char d_clean[5];
        strip_whitespace_fixed(rd->atomNames[a], 4, d_clean, 4);
        if (strcmp(atom_name, d_clean) == 0) {
            return a;
        }
    }
    return -1;
}

uint8_t res3_to_type_id(const char* src, int width) {
    char buf[3] = {' ', ' ', ' '};
    int d = 0;
    for (int i = 0; i < width && d < 3; i++) {
        if (src[i] != ' ' && src[i] != '\n' && src[i] != '\r' && src[i] != '\0') {
            buf[d++] = src[i];
        }
    }

    uint32_t val = (uint8_t)buf[0] | ((uint32_t)(uint8_t)buf[1] << 8) | ((uint32_t)(uint8_t)buf[2] << 16);
    char one_letter;

    switch (val) {
        case 'A' | ('L' << 8) | ('A' << 16): one_letter = 'A'; break;
        case 'A' | ('R' << 8) | ('G' << 16): one_letter = 'R'; break;
        case 'A' | ('S' << 8) | ('N' << 16): one_letter = 'N'; break;
        case 'A' | ('S' << 8) | ('P' << 16): one_letter = 'D'; break;
        case 'C' | ('Y' << 8) | ('S' << 16): one_letter = 'C'; break;
        case 'G' | ('L' << 8) | ('N' << 16): one_letter = 'Q'; break;
        case 'G' | ('L' << 8) | ('U' << 16): one_letter = 'E'; break;
        case 'G' | ('L' << 8) | ('Y' << 16): one_letter = 'G'; break;
        case 'H' | ('I' << 8) | ('S' << 16): one_letter = 'H'; break;
        case 'I' | ('L' << 8) | ('E' << 16): one_letter = 'I'; break;
        case 'L' | ('E' << 8) | ('U' << 16): one_letter = 'L'; break;
        case 'L' | ('Y' << 8) | ('S' << 16): one_letter = 'K'; break;
        case 'M' | ('E' << 8) | ('T' << 16): one_letter = 'M'; break;
        case 'P' | ('H' << 8) | ('E' << 16): one_letter = 'F'; break;
        case 'P' | ('R' << 8) | ('O' << 16): one_letter = 'P'; break;
        case 'S' | ('E' << 8) | ('R' << 16): one_letter = 'S'; break;
        case 'T' | ('H' << 8) | ('R' << 16): one_letter = 'T'; break;
        case 'T' | ('R' << 8) | ('P' << 16): one_letter = 'W'; break;
        case 'T' | ('Y' << 8) | ('R' << 16): one_letter = 'Y'; break;
        case 'V' | ('A' << 8) | ('L' << 16): one_letter = 'V'; break;
        default:                             one_letter = 'X'; break;
    }
    return (uint8_t)(one_letter - 0x40);
}

const char* type_id_to_res3(uint8_t res_type) {
    switch ((char)(res_type + 0x40)) {
        case 'A': return "ALA";
        case 'R': return "ARG";
        case 'N': return "ASN";
        case 'D': return "ASP";
        case 'C': return "CYS";
        case 'Q': return "GLN";
        case 'E': return "GLU";
        case 'G': return "GLY";
        case 'H': return "HIS";
        case 'I': return "ILE";
        case 'L': return "LEU";
        case 'K': return "LYS";
        case 'M': return "MET";
        case 'F': return "PHE";
        case 'P': return "PRO";
        case 'S': return "SER";
        case 'T': return "THR";
        case 'W': return "TRP";
        case 'Y': return "TYR";
        case 'V': return "VAL";
        default:  return "UNK";
    }
}
