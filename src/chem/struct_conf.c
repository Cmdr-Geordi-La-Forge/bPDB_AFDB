#include "chem/struct_conf.h"
#include "bpdb_internal.h"
#include <string.h>

uint8_t encode_struct_conf_type(const char* type_str) {
    if (type_str[0] == 'O') return TYPE_OTHER;
    if (type_str[0] == 'B') return CAT_BEND;
    if (type_str[0] == 'S') return CAT_STRAND;

    if (type_str[0] == 'T') {
        if (strcmp(type_str, "TURN_TY1_P")  == 0) return CAT_TURN | 1;
        if (strcmp(type_str, "TURN_TY1P_P") == 0) return CAT_TURN | 2;
        if (strcmp(type_str, "TURN_TY2_P")  == 0) return CAT_TURN | 3;
        if (strcmp(type_str, "TURN_TY2P_P") == 0) return CAT_TURN | 4;
        if (strcmp(type_str, "TURN_TY3_P")  == 0) return CAT_TURN | 5;
        if (strcmp(type_str, "TURN_TY3P_P") == 0) return CAT_TURN | 6;
        return CAT_TURN | 7;
    }

    if (type_str[0] == 'H') {
        if (strcmp(type_str, "HELX_RH_AL_P") == 0) return CAT_HELIX | (1 << 1) | 0;
        if (strcmp(type_str, "HELX_LH_AL_P") == 0) return CAT_HELIX | (1 << 1) | 1;
        if (strcmp(type_str, "HELX_RH_OM_P") == 0) return CAT_HELIX | (2 << 1) | 0;
        if (strcmp(type_str, "HELX_LH_OM_P") == 0) return CAT_HELIX | (2 << 1) | 1;
        if (strcmp(type_str, "HELX_RH_PI_P") == 0) return CAT_HELIX | (3 << 1) | 0;
        if (strcmp(type_str, "HELX_LH_PI_P") == 0) return CAT_HELIX | (3 << 1) | 1;
        if (strcmp(type_str, "HELX_RH_GA_P") == 0) return CAT_HELIX | (4 << 1) | 0;
        if (strcmp(type_str, "HELX_LH_GA_P") == 0) return CAT_HELIX | (4 << 1) | 1;
        if (strcmp(type_str, "HELX_RH_3T_P") == 0) return CAT_HELIX | (5 << 1) | 0;
        if (strcmp(type_str, "HELX_LH_3T_P") == 0) return CAT_HELIX | (5 << 1) | 1;
        if (strcmp(type_str, "HELX_RH_27_P") == 0) return CAT_HELIX | (9 << 1) | 0;
        if (strcmp(type_str, "HELX_LH_27_P") == 0) return CAT_HELIX | (9 << 1) | 1;
        if (strcmp(type_str, "HELX_RH_PP_P") == 0) return CAT_HELIX | (10 << 1) | 0;
        if (strcmp(type_str, "HELX_LH_PP_P") == 0) return CAT_HELIX | (10 << 1) | 1;
        return CAT_HELIX | (15 << 1);
    }

    return TYPE_OTHER;
}

const char* decode_struct_conf_type(uint8_t type_id) {
    if (type_id == TYPE_OTHER) return "OTHER";

    uint8_t cat = type_id & CONF_CAT_MASK;
    uint8_t sub = type_id & CONF_SUB_MASK;

    if (cat == CAT_BEND)   return "BEND";
    if (cat == CAT_STRAND) return "STRN";

    if (cat == CAT_TURN) {
        switch (sub) {
            case 1:  return "TURN_TY1_P";
            case 2:  return "TURN_TY1P_P";
            case 3:  return "TURN_TY2_P";
            case 4:  return "TURN_TY2P_P";
            case 5:  return "TURN_TY3_P";
            case 6:  return "TURN_TY3P_P";
            default: return "TURN_P";
        }
    }

    if (cat == CAT_HELIX) {
        uint8_t base = sub >> 1;
        uint8_t lh   = sub & 1;

        switch (base) {
            case 1:  return lh ? "HELX_LH_AL_P" : "HELX_RH_AL_P";
            case 2:  return lh ? "HELX_LH_OM_P" : "HELX_RH_OM_P";
            case 3:  return lh ? "HELX_LH_PI_P" : "HELX_RH_PI_P";
            case 4:  return lh ? "HELX_LH_GA_P" : "HELX_RH_GA_P";
            case 5:  return lh ? "HELX_LH_3T_P" : "HELX_RH_3T_P";
            case 9:  return lh ? "HELX_LH_27_P" : "HELX_RH_27_P";
            case 10: return lh ? "HELX_LH_PP_P" : "HELX_RH_PP_P";
            default: return "HELX_P";
        }
    }

    return "OTHER";
}
