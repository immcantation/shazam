# Generates and loads region definitions


# IMGT Numbering for V gene (minus CDR3, i.e. nucleotide positions 1-312)
# 
# Lefranc MP, Pommie C, Ruiz M, Giudicelli V, Foulquier E, Truong L, 
# Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin 
# and T cell receptor variable domains and Ig superfamily V-like domains.
# Developmental and comparative immunology. 2003;27:55-77.
IMGT_V_NO_CDR3 <- createRegionDefinition(name="IMGT_V_NO_CDR3",
                                         boundaries=factor( c( rep("FWR", 78), 
                                                               rep("CDR", 36),  
                                                               rep("FWR", 51), 
                                                               rep("CDR", 30), 
                                                               rep("FWR", 117) ),
                                                            levels = c("CDR","FWR")),
                                         description="IMGT_Numbering scheme defining the V gene up till but not including CDR3",
                                         citation="Lefranc MP, Pommie C, Ruiz M, Giudicelli V, Foulquier E, Truong L, Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains. Developmental and comparative immunology. 2003;27:55-77.") 
devtools::use_data(IMGT_V_NO_CDR3, overwrite=TRUE)


# IMGT Numbering for V gene.
#
# Lefranc MP, Pommie C, Ruiz M, Giudicelli V, Foulquier E, Truong L, 
# Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin 
# and T cell receptor variable domains and Ig superfamily V-like domains.
# Developmental and comparative immunology. 2003;27:55-77.
IMGT_V <- createRegionDefinition(name="IMGT_V_NO_CDR3",
                                 boundaries=factor( c( rep("FWR", 78), 
                                                       rep("CDR", 36),  
                                                       rep("FWR", 51), 
                                                       rep("CDR", 30), 
                                                       rep("FWR", 117),
                                                       rep("CDR", 36)),
                                                    levels = c("CDR","FWR")),
                                 description="IMGT_Numbering scheme defining the V gene up till but not including CDR3",
                                 citation="Lefranc MP, Pommie C, Ruiz M, Giudicelli V, Foulquier E, Truong L, Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains. Developmental and comparative immunology. 2003;27:55-77.") 
devtools::use_data(IMGT_V, overwrite=TRUE)


# IMGT Numbering for V gene (minus CDR3), broken down by the individual regions.
# I.e. FWR1, CDR1, FWR2, CDR2, FWR3, CDR3.
#
# Lefranc MP, Pommie C, Ruiz M, Giudicelli V, Foulquier E, Truong L, 
# Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin 
# and T cell receptor variable domains and Ig superfamily V-like domains.
# Developmental and comparative immunology. 2003;27:55-77.
IMGT_V_BY_REGIONS_NO_CDR3 <- createRegionDefinition(name="IMGT_V_NO_CDR3",
                                                    boundaries=factor( c( rep("FWR1", 78), 
                                                                          rep("CDR1", 36),  
                                                                          rep("FWR2", 51), 
                                                                          rep("CDR2", 30), 
                                                                          rep("FWR3", 117)),
                                                                       levels = c(paste0("CDR",1:2),paste0("FWR",1:3)) ),
                                                    description="IMGT_Numbering scheme defining the V gene up till but not including CDR3",
                                                    citation="Lefranc MP, Pommie C, Ruiz M, Giudicelli V, Foulquier E, Truong L, Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains. Developmental and comparative immunology. 2003;27:55-77.") 
devtools::use_data(IMGT_V_BY_REGIONS_NO_CDR3, overwrite=TRUE)


# IMGT Numbering for V gene, broken down by the individual regions.
# I.e. FWR1, CDR1, FWR2, CDR2, FWR3, CDR3.
#
# Lefranc MP, Pommie C, Ruiz M, Giudicelli V, Foulquier E, Truong L, 
# Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin 
# and T cell receptor variable domains and Ig superfamily V-like domains.
# Developmental and comparative immunology. 2003;27:55-77.
IMGT_V_BY_REGIONS <- createRegionDefinition(name="IMGT_V_NO_CDR3",
                                            boundaries=factor( c( rep("FWR1", 78), 
                                                                  rep("CDR1", 36),  
                                                                  rep("FWR2", 51), 
                                                                  rep("CDR2", 30), 
                                                                  rep("FWR3", 117),
                                                                  rep("CDR3", 36)),
                                                               levels = c(paste0("CDR",1:3),paste0("FWR",1:3)) ),
                                            description="IMGT_Numbering scheme defining the V gene up till but not including CDR3",
                                            citation="Lefranc MP, Pommie C, Ruiz M, Giudicelli V, Foulquier E, Truong L, Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains. Developmental and comparative immunology. 2003;27:55-77.") 
devtools::use_data(IMGT_V_BY_REGIONS, overwrite=TRUE)
