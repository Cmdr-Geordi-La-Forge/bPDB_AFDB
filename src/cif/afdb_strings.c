#include "cif/afdb_strings.h"

const ChemCompRow CHEM_COMP_ROWS[20] = {
    {'A', "\"C3 H7 N O2\"",   "89.093",  "ALA", "ALANINE",           "\"L-PEPTIDE LINKING\""},
    {'R', "\"C6 H15 N4 O2\"", "175.209", "ARG", "ARGININE",          "\"L-PEPTIDE LINKING\""},
    {'N', "\"C4 H8 N2 O3\"",  "132.118", "ASN", "ASPARAGINE",        "\"L-PEPTIDE LINKING\""},
    {'D', "\"C4 H7 N O4\"",   "133.103", "ASP", "\"ASPARTIC ACID\"", "\"L-PEPTIDE LINKING\""},
    {'C', "\"C3 H7 N O2 S\"", "121.158", "CYS", "CYSTEINE",          "\"L-PEPTIDE LINKING\""},
    {'Q', "\"C5 H10 N2 O3\"", "146.144", "GLN", "GLUTAMINE",         "\"L-PEPTIDE LINKING\""},
    {'E', "\"C5 H9 N O4\"",   "147.129", "GLU", "\"GLUTAMIC ACID\"", "\"L-PEPTIDE LINKING\""},
    {'G', "\"C2 H5 N O2\"",   "75.067",  "GLY", "GLYCINE",           "\"PEPTIDE LINKING\""},
    {'H', "\"C6 H10 N3 O2\"", "156.162", "HIS", "HISTIDINE",         "\"L-PEPTIDE LINKING\""},
    {'I', "\"C6 H13 N O2\"",  "131.173", "ILE", "ISOLEUCINE",        "\"L-PEPTIDE LINKING\""},
    {'L', "\"C6 H13 N O2\"",  "131.173", "LEU", "LEUCINE",           "\"L-PEPTIDE LINKING\""},
    {'K', "\"C6 H15 N2 O2\"", "147.195", "LYS", "LYSINE",            "\"L-PEPTIDE LINKING\""},
    {'M', "\"C5 H11 N O2 S\"","149.211", "MET", "METHIONINE",        "\"L-PEPTIDE LINKING\""},
    {'F', "\"C9 H11 N O2\"",  "165.189", "PHE", "PHENYLALANINE",     "\"L-PEPTIDE LINKING\""},
    {'P', "\"C5 H9 N O2\"",   "115.130", "PRO", "PROLINE",           "\"L-PEPTIDE LINKING\""},
    {'S', "\"C3 H7 N O3\"",   "105.093", "SER", "SERINE",            "\"L-PEPTIDE LINKING\""},
    {'T', "\"C4 H9 N O3\"",   "119.119", "THR", "THREONINE",         "\"L-PEPTIDE LINKING\""},
    {'W', "\"C11 H12 N2 O2\"","204.225", "TRP", "TRYPTOPHAN",        "\"L-PEPTIDE LINKING\""},
    {'Y', "\"C9 H11 N O3\"",  "181.189", "TYR", "TYROSINE",          "\"L-PEPTIDE LINKING\""},
    {'V', "\"C5 H11 N O2\"",  "117.146", "VAL", "VALINE",            "\"L-PEPTIDE LINKING\""}
};

const size_t NUM_CHEM_COMP_ROWS = sizeof(CHEM_COMP_ROWS) / sizeof(CHEM_COMP_ROWS[0]);

const char* AFDB_ATOM_TYPE_HEADER = 
    "loop_\n"
    "_atom_type.symbol\n"
    "C \n"
    "N \n"
    "O \n";

const char* AFDB_CHEM_COMP_HEADER = 
    "loop_\n"
    "_chem_comp.formula\n"
    "_chem_comp.formula_weight\n"
    "_chem_comp.id\n"
    "_chem_comp.mon_nstd_flag\n"
    "_chem_comp.name\n"
    "_chem_comp.pdbx_synonyms\n"
    "_chem_comp.type\n";

const char* AFDB_ENTITY_BLOCK_1 = 
    "_entity.details                  ?\n"
    "_entity.formula_weight           ?\n"
    "_entity.id                       1\n";

const char* AFDB_ENTITY_BLOCK_2 = 
    "_entity.pdbx_ec                  ?\n"
    "_entity.pdbx_fragment            ?\n"
    "_entity.pdbx_mutation            ?\n"
    "_entity.pdbx_number_of_molecules 1\n"
    "_entity.src_method               man\n"
    "_entity.type                     polymer\n"
    "#\n"
    "_entity_poly.entity_id                    1\n"
    "_entity_poly.nstd_linkage                 no\n"
    "_entity_poly.nstd_monomer                 no\n";

const char* AFDB_ENTITY_POLY_SEQ_HEADER = 
    "loop_\n"
    "_entity_poly_seq.entity_id\n"
    "_entity_poly_seq.hetero\n"
    "_entity_poly_seq.mon_id\n"
    "_entity_poly_seq.num\n";

const char* AFDB_MA_DATA_MULTI = 
    "loop_\n"
    "_ma_data.content_type\n"
    "_ma_data.id\n"
    "_ma_data.name\n"
    "\"model coordinates\" 1 Model             \n"
    "\"input structure\"   2 \"Input structure\" \n"
    "#\n";

const char* AFDB_MA_DATA_SINGLE = 
    "_ma_data.content_type \"model coordinates\"\n"
    "_ma_data.id           1\n"
    "_ma_data.name         Model\n"
    "#\n";

const char* AFDB_QA_METRIC_LOCAL_HEADER = 
    "loop_\n"
    "_ma_qa_metric_local.label_asym_id\n"
    "_ma_qa_metric_local.label_comp_id\n"
    "_ma_qa_metric_local.label_seq_id\n"
    "_ma_qa_metric_local.metric_id\n"
    "_ma_qa_metric_local.metric_value\n"
    "_ma_qa_metric_local.model_id\n"
    "_ma_qa_metric_local.ordinal_id\n";

const char* AFDB_AUDIT_REVISION_DETAILS_HEADER = 
    "loop_\n"
    "_pdbx_audit_revision_details.data_content_type\n"
    "_pdbx_audit_revision_details.description\n"
    "_pdbx_audit_revision_details.ordinal\n"
    "_pdbx_audit_revision_details.provider\n"
    "_pdbx_audit_revision_details.revision_ordinal\n"
    "_pdbx_audit_revision_details.type\n";

const char* AFDB_AUDIT_REVISION_HISTORY_HEADER = 
    "loop_\n"
    "_pdbx_audit_revision_history.data_content_type\n"
    "_pdbx_audit_revision_history.major_revision\n"
    "_pdbx_audit_revision_history.minor_revision\n"
    "_pdbx_audit_revision_history.ordinal\n"
    "_pdbx_audit_revision_history.revision_date\n";

const char* AFDB_POLY_SEQ_SCHEME_HEADER = 
    "loop_\n"
    "_pdbx_poly_seq_scheme.asym_id\n"
    "_pdbx_poly_seq_scheme.auth_seq_num\n"
    "_pdbx_poly_seq_scheme.entity_id\n"
    "_pdbx_poly_seq_scheme.hetero\n"
    "_pdbx_poly_seq_scheme.mon_id\n"
    "_pdbx_poly_seq_scheme.pdb_ins_code\n"
    "_pdbx_poly_seq_scheme.pdb_mon_id\n"
    "_pdbx_poly_seq_scheme.pdb_seq_num\n"
    "_pdbx_poly_seq_scheme.pdb_strand_id\n"
    "_pdbx_poly_seq_scheme.seq_id\n";

const char* AFDB_AUDIT_BLOCK = 
    "loop_\n"
    "_audit_author.name\n"
    "_audit_author.pdbx_ordinal\n"
    "\"Jumper, John\"               1  \n"
    "\"Evans, Richard\"             2  \n"
    "\"Pritzel, Alexander\"         3  \n"
    "\"Green, Tim\"                 4  \n"
    "\"Figurnov, Michael\"          5  \n"
    "\"Ronneberger, Olaf\"          6  \n"
    "\"Tunyasuvunakool, Kathryn\"   7  \n"
    "\"Bates, Russ\"                8  \n"
    "\"Zidek, Augustin\"            9  \n"
    "\"Potapenko, Anna\"            10 \n"
    "\"Bridgland, Alex\"            11 \n"
    "\"Meyer, Clemens\"             12 \n"
    "\"Kohl, Simon A. A.\"          13 \n"
    "\"Ballard, Andrew J.\"         14 \n"
    "\"Cowie, Andrew\"              15 \n"
    "\"Romera-Paredes, Bernardino\" 16 \n"
    "\"Nikolov, Stanislav\"         17 \n"
    "\"Jain, Rishub\"               18 \n"
    "\"Adler, Jonas\"               19 \n"
    "\"Back, Trevor\"               20 \n"
    "\"Petersen, Stig\"             21 \n"
    "\"Reiman, David\"              22 \n"
    "\"Clancy, Ellen\"              23 \n"
    "\"Zielinski, Michal\"          24 \n"
    "\"Steinegger, Martin\"         25 \n"
    "\"Pacholska, Michalina\"       26 \n"
    "\"Berghammer, Tamas\"          27 \n"
    "\"Silver, David\"              28 \n"
    "\"Vinyals, Oriol\"             29 \n"
    "\"Senior, Andrew W.\"          30 \n"
    "\"Kavukcuoglu, Koray\"         31 \n"
    "\"Kohli, Pushmeet\"            32 \n"
    "\"Hassabis, Demis\"            33 \n"
    "#\n"
    "_audit_conform.dict_location https://raw.githubusercontent.com/ihmwg/ModelCIF/master/dist/mmcif_ma.dic\n"
    "_audit_conform.dict_name     mmcif_ma.dic\n"
    "_audit_conform.dict_version  1.3.9\n"
    "#\n";

const char* AFDB_CITATION_BLOCK = 
    "_citation.book_publisher          ?\n"
    "_citation.country                 UK\n"
    "_citation.id                      primary\n"
    "_citation.journal_full            Nature\n"
    "_citation.journal_id_ASTM         NATUAS\n"
    "_citation.journal_id_CSD          0006\n"
    "_citation.journal_id_ISSN         0028-0836\n"
    "_citation.journal_volume          596\n"
    "_citation.page_first              583\n"
    "_citation.page_last               589\n"
    "_citation.pdbx_database_id_DOI    10.1038/s41586-021-03819-2\n"
    "_citation.pdbx_database_id_PubMed 34265844\n"
    "_citation.title                   \"Highly accurate protein structure prediction with AlphaFold\"\n"
    "_citation.year                    2021\n"
    "#\n"
    "loop_\n"
    "_citation_author.citation_id\n"
    "_citation_author.name\n"
    "_citation_author.ordinal\n"
    "primary \"Jumper, John\"               1  \n"
    "primary \"Evans, Richard\"             2  \n"
    "primary \"Pritzel, Alexander\"         3  \n"
    "primary \"Green, Tim\"                 4  \n"
    "primary \"Figurnov, Michael\"          5  \n"
    "primary \"Ronneberger, Olaf\"          6  \n"
    "primary \"Tunyasuvunakool, Kathryn\"   7  \n"
    "primary \"Bates, Russ\"                8  \n"
    "primary \"Zidek, Augustin\"            9  \n"
    "primary \"Potapenko, Anna\"            10 \n"
    "primary \"Bridgland, Alex\"            11 \n"
    "primary \"Meyer, Clemens\"             12 \n"
    "primary \"Kohl, Simon A. A.\"          13 \n"
    "primary \"Ballard, Andrew J.\"         14 \n"
    "primary \"Cowie, Andrew\"              15 \n"
    "primary \"Romera-Paredes, Bernardino\" 16 \n"
    "primary \"Nikolov, Stanislav\"         17 \n"
    "primary \"Jain, Rishub\"               18 \n"
    "primary \"Adler, Jonas\"               19 \n"
    "primary \"Back, Trevor\"               20 \n"
    "primary \"Petersen, Stig\"             21 \n"
    "primary \"Reiman, David\"              22 \n"
    "primary \"Clancy, Ellen\"              23 \n"
    "primary \"Zielinski, Michal\"          24 \n"
    "primary \"Steinegger, Martin\"         25 \n"
    "primary \"Pacholska, Michalina\"       26 \n"
    "primary \"Berghammer, Tamas\"          27 \n"
    "primary \"Silver, David\"              28 \n"
    "primary \"Vinyals, Oriol\"             29 \n"
    "primary \"Senior, Andrew W.\"          30 \n"
    "primary \"Kavukcuoglu, Koray\"         31 \n"
    "primary \"Kohli, Pushmeet\"            32 \n"
    "primary \"Hassabis, Demis\"            33 \n"
    "#\n";

const char* AFDB_MA_MODEL_LIST = 
    "_ma_model_list.data_id          1\n"
    "_ma_model_list.model_group_id   1\n"
    "_ma_model_list.model_group_name \"AlphaFold Monomer v2.0 model\"\n"
    "_ma_model_list.model_id         1\n"
    "_ma_model_list.model_name       \"Top ranked model\"\n"
    "_ma_model_list.model_type       \"Ab initio model\"\n"
    "_ma_model_list.ordinal_id       1\n"
    "#\n"
    "loop_\n"
    "_ma_protocol_step.method_type\n"
    "_ma_protocol_step.ordinal_id\n"
    "_ma_protocol_step.protocol_id\n"
    "_ma_protocol_step.step_id\n"
    "\"coevolution MSA\" 1 1 1 \n"
    "\"template search\" 2 1 2 \n"
    "modeling          3 1 3 \n"
    "#\n"
    "loop_\n"
    "_ma_qa_metric.id\n"
    "_ma_qa_metric.mode\n"
    "_ma_qa_metric.name\n"
    "_ma_qa_metric.software_group_id\n"
    "_ma_qa_metric.type\n"
    "1 global pLDDT 1 pLDDT \n"
    "2 local  pLDDT 1 pLDDT \n"
    "#\n";

const char* AFDB_MA_POST_BLOCK = 
    "_ma_software_group.group_id    1\n"
    "_ma_software_group.ordinal_id  1\n"
    "_ma_software_group.software_id 1\n"
    "#\n"
    "_ma_target_entity.data_id   1\n"
    "_ma_target_entity.entity_id 1\n"
    "_ma_target_entity.origin    \"reference database\"\n"
    "#\n"
    "_ma_target_entity_instance.asym_id   A\n"
    "_ma_target_entity_instance.details   .\n"
    "_ma_target_entity_instance.entity_id 1\n"
    "#\n";

const char* AFDB_MA_TEMPLATE_DETAILS_HEADER = 
    "loop_\n"
    "_ma_template_details.ordinal_id\n"
    "_ma_template_details.target_asym_id\n"
    "_ma_template_details.template_auth_asym_id\n"
    "_ma_template_details.template_data_id\n"
    "_ma_template_details.template_entity_type\n"
    "_ma_template_details.template_id\n"
    "_ma_template_details.template_model_num\n"
    "_ma_template_details.template_origin\n"
    "_ma_template_details.template_trans_matrix_id\n";

const char* AFDB_MA_TEMPLATE_REF_DB_DETAILS_HEADER = 
    "loop_\n"
    "_ma_template_ref_db_details.db_accession_code\n"
    "_ma_template_ref_db_details.db_name\n"
    "_ma_template_ref_db_details.template_id\n";

const char* AFDB_MA_TEMPLATE_TRANS_MATRIX = 
    "_ma_template_trans_matrix.id               1\n"
    "_ma_template_trans_matrix.rot_matrix[1][1] 1.0\n"
    "_ma_template_trans_matrix.rot_matrix[1][2] 0.0\n"
    "_ma_template_trans_matrix.rot_matrix[1][3] 0.0\n"
    "_ma_template_trans_matrix.rot_matrix[2][1] 0.0\n"
    "_ma_template_trans_matrix.rot_matrix[2][2] 1.0\n"
    "_ma_template_trans_matrix.rot_matrix[2][3] 0.0\n"
    "_ma_template_trans_matrix.rot_matrix[3][1] 0.0\n"
    "_ma_template_trans_matrix.rot_matrix[3][2] 0.0\n"
    "_ma_template_trans_matrix.rot_matrix[3][3] 1.0\n"
    "_ma_template_trans_matrix.tr_vector[1]     0.0\n"
    "_ma_template_trans_matrix.tr_vector[2]     0.0\n"
    "_ma_template_trans_matrix.tr_vector[3]     0.0\n"
    "#\n";

const char* AFDB_PDBX_DATA_USAGE = 
    "loop_\n"
    "_pdbx_data_usage.details\n"
    "_pdbx_data_usage.id\n"
    "_pdbx_data_usage.name\n"
    "_pdbx_data_usage.type\n"
    "_pdbx_data_usage.url\n"
    "\"Data in this file is available under a CC-BY-4.0 license.\" 1 CC-BY-4.0 license    https://creativecommons.org/licenses/by/4.0/ \n"
    ";ALPHAFOLD DATA, COPYRIGHT (2021) DEEPMIND TECHNOLOGIES LIMITED. THE INFORMATION\n"
    "PROVIDED IS THEORETICAL MODELLING ONLY AND CAUTION SHOULD BE EXERCISED IN ITS\n"
    "USE. IT IS PROVIDED \"AS-IS\" WITHOUT ANY WARRANTY OF ANY KIND, WHETHER EXPRESSED\n"
    "OR IMPLIED. NO WARRANTY IS GIVEN THAT USE OF THE INFORMATION SHALL NOT INFRINGE\n"
    "THE RIGHTS OF ANY THIRD PARTY. DISCLAIMER: THE INFORMATION IS NOT INTENDED TO BE\n"
    "A SUBSTITUTE FOR PROFESSIONAL MEDICAL ADVICE, DIAGNOSIS, OR TREATMENT, AND DOES\n"
    "NOT CONSTITUTE MEDICAL OR OTHER PROFESSIONAL ADVICE. IT IS AVAILABLE FOR\n"
    "ACADEMIC AND COMMERCIAL PURPOSES, UNDER CC-BY 4.0 LICENCE.\n"
    ";\n"
    "2 ?         disclaimer ?                                            \n"
    "#\n";

const char* AFDB_SOFTWARE_ASYM = 
    "loop_\n"
    "_software.classification\n"
    "_software.date\n"
    "_software.description\n"
    "_software.name\n"
    "_software.pdbx_ordinal\n"
    "_software.type\n"
    "_software.version\n"
    "other ? \"Structure prediction\" AlphaFold 1 package v2.0 \n"
    "other ? \"Secondary structure\"  dssp      2 library 4    \n"
    "#\n";

const char* AFDB_STRUCT_CONF_HEADER =
    "loop_\n"
    "_struct_conf.beg_auth_asym_id\n"
    "_struct_conf.beg_auth_comp_id\n"
    "_struct_conf.beg_auth_seq_id\n"
    "_struct_conf.beg_label_asym_id\n"
    "_struct_conf.beg_label_comp_id\n"
    "_struct_conf.beg_label_seq_id\n"
    "_struct_conf.conf_type_id\n"
    "_struct_conf.end_auth_asym_id\n"
    "_struct_conf.end_auth_comp_id\n"
    "_struct_conf.end_auth_seq_id\n"
    "_struct_conf.end_label_asym_id\n"
    "_struct_conf.end_label_comp_id\n"
    "_struct_conf.end_label_seq_id\n"
    "_struct_conf.id\n"
    "_struct_conf.pdbx_beg_PDB_ins_code\n"
    "_struct_conf.pdbx_end_PDB_ins_code\n";

const char* AFDB_STRUCT_CONF_TYPE_HEADER =
    "loop_\n"
    "_struct_conf_type.criteria\n"
    "_struct_conf_type.id\n";

const char* AFDB_ATOM_SITE_HEADER =
    "loop_\n"
    "_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n"
    "_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n"
    "_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n"
    "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n"
    "_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n"
    "_atom_site.pdbx_formal_charge\n_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n"
    "_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n_atom_site.pdbx_PDB_model_num\n"
    "_atom_site.pdbx_sifts_xref_db_acc\n_atom_site.pdbx_sifts_xref_db_name\n"
    "_atom_site.pdbx_sifts_xref_db_num\n_atom_site.pdbx_sifts_xref_db_res\n";
