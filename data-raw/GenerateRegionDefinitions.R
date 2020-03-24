# Generates region definitions

# IMGT numbering for V segment minus CDR3; i.e. nucleotide positions 1-312.
# CDR, FWR
#
# Lefranc MP, Pommie C, Ruiz M, Giudicelli V, Foulquier E, Truong L, 
# Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin 
# and T cell receptor variable domains and Ig superfamily V-like domains.
# Developmental and comparative immunology. 2003;27:55-77.
IMGT_V <- createRegionDefinition(name="IMGT_V",
                                 boundaries=factor(c(rep("fwr", 78), 
                                                     rep("cdr", 36),  
                                                     rep("fwr", 51), 
                                                     rep("cdr", 30), 
                                                     rep("fwr", 117)),
                                                   levels = c("cdr","fwr")),
                                 description="IMGT numbering scheme defining the V segment, excluding CDR3.",
                                 citation="Lefranc MP, Pommie C, Ruiz M, Giudicelli V, Foulquier E, Truong L, Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains. Developmental and comparative immunology. 2003;27:55-77.") 
usethis::use_data(IMGT_V, overwrite=TRUE)

# IMGT numbering for V segment, broken down by the individual CDR and FWR regions.
# FWR1, CDR1, FWR2, CDR2, FWR3
IMGT_V_BY_REGIONS <- createRegionDefinition(name="IMGT_V_BY_REGIONS",
                                            boundaries=factor(c(rep("fwr1", 78), 
                                                                 rep("cdr1", 36),  
                                                                 rep("fwr2", 51), 
                                                                 rep("cdr2", 30), 
                                                                 rep("fwr3", 117)),
                                                              levels = c(paste0("cdr",1:2), paste0("fwr",1:3))),
                                            description="IMGT numbering scheme defining the V segment by individual CDR and FWR regions, excluding CDR3",
                                            citation="Lefranc MP, Pommie C, Ruiz M, Giudicelli V, Foulquier E, Truong L, Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains. Developmental and comparative immunology. 2003;27:55-77.") 
usethis::use_data(IMGT_V_BY_REGIONS, overwrite=TRUE)

# IMGT numbering for V segment broken down by individual codons
# 
IMGT_V_BY_CODONS <- createRegionDefinition(name="IMGT_V_BY_CODONS",
                                           boundaries=factor(rep(as.character(1:104), each=3), 
                                                             levels=as.character(1:104)),
                                           description="IMGT numbering scheme defining the V segment by individual codons, excluding CDR3.",
                                           citation="Lefranc MP, Pommie C, Ruiz M, Giudicelli V, Foulquier E, Truong L, Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains. Developmental and comparative immunology. 2003;27:55-77.") 
usethis::use_data(IMGT_V_BY_CODONS, overwrite=TRUE)

# IMGT numbering for V segment with a single region for the entire V.
# 
IMGT_V_BY_SEGMENTS <- createRegionDefinition(name="IMGT_V_BY_SEGMENTS",
                                             boundaries=factor(rep("v", 312), 
                                                               levels="v"),
                                             description="IMGT numbering scheme defining the V segment by individual codons, excluding CDR3.",
                                             citation="Lefranc MP, Pommie C, Ruiz M, Giudicelli V, Foulquier E, Truong L, Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains. Developmental and comparative immunology. 2003;27:55-77.") 
usethis::use_data(IMGT_V_BY_SEGMENTS, overwrite=TRUE)
