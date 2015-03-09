# Project documentation for shm
# 
# @author     Mohamed Uduman, Gur Yaari, Namita Gupta, Jason Anthony Vander Heiden
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.1
# @date       2015.03.07


#' The shm package
#'
#' Provides tools for advanced anaylisis of Ig Somatic HyperMutation. Includes BASELINe,
#' a novel method for quantifying selection in high-throughput Immunoglobulin sequencing
#' data sets.
#' 
#' @section  Selection analysis:
#' \itemize{
#'   \item  \code{\link{addClonalSequence}}:        Build clonal consensus sequence.
#'   \item  \code{\link{addExpectedFrequencies}}:   Compute expected mutation frequencies.
#'   \item  \code{\link{addObservedMutations}}:     Compute observed mutation counts.
#'   \item  \code{\link{computeBaselinePDF}}:       Compute selection strength.
#'   \item  \code{\link{plotSelection}}:            Plot selection strength.
#' }
#'
#' @section  Targeting models:
#' \itemize{
#'   \item  \code{\link{createMutabilityModel}}:    Builds a mutability model.
#'   \item  \code{\link{createSubstitutionModel}}:  Builds a substitution model.
#'   \item  \code{\link{plotHedgehog}}:             Plots a targeting model.
#' }
#'
#' @section  Distance profiling:
#' \itemize{
#'   \item  \code{\link{distToNearest}}:            Calculate distances to nearest-neighbors.
#'   \item  \code{\link{getPairwiseDistances}}:     Calculate a matrix of pairwise distances.
#' }
#'
#' @references
#' \enumerate{
#'   \item  Uduman M, et al. Detecting selection in immunoglobulin sequences. 
#'            Nucleic Acids Res. 2011 39(Web Server issue):W499–504.
#'   \item  Yaari G, et al. Quantifying selection in high-throughput Immunoglobulin 
#'            sequencing data sets. 
#'            Nucleic Acids Res. 2012 40(17):e134.
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution 
#'            based on synonymous mutations from high-throughput immunoglobulin sequencing 
#'            data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#'
#' @seealso
#' The Change-O suite of tools includes three separate R packages: \code{\link{alakazam}}, 
#' \code{\link{tigger}}, and \code{\link{shm}}.
#' 
#' @name     shm
#' @docType  package
#'
#' @import   alakazam
#' @import   data.table
#' @import   doSNOW
#' @import   foreach
#' @import   ggplot2
#' @import   grid
#' @import   plyr
#' @import   seqinr
#' @import   SDMTools
#' @import   zoo
NULL


#### Constants ####

# IMGT V-region definitions
VREGIONS <- factor(c(rep("FWR", 78), 
                     rep("CDR", 36), 
                     rep("FWR", 51), 
                     rep("CDR", 30), 
                     rep("FWR", 117)), 
                   levels=c("FWR", "CDR"))

# IMGT V-region length
#readEnd <- 312
VLENGTH <- length(VREGIONS)

# Nucleotide characters
NUCLEOTIDES <- c("A", "C", "G", "T")

# Amino acid characters
AMINO_ACIDS <- c("F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y", "*", "*", "C", "C", "*", "W", "L", "L", "L", "L", "P", "P", "P", "P", "H", "H", "Q", "Q", "R", "R", "R", "R", "I", "I", "I", "M", "T", "T", "T", "T", "N", "N", "K", "K", "S", "S", "R", "R", "V", "V", "V", "V", "A", "A", "A", "A", "D", "D", "E", "E", "G", "G", "G", "G")
names(AMINO_ACIDS) <- c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")
names(AMINO_ACIDS) <- names(AMINO_ACIDS)

#Amino Acid Traits
#"*" "A" "C" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y"
#B = "Hydrophobic/Burried"  N = "Intermediate/Neutral"  S="Hydrophilic/Surface")
TRAITS_AMINO_ACIDS_CHOTHIA98 <- c("*","N","B","S","S","B","N","N","B","S","B","B","S","N","S","S","N","N","B","B","N")
names(TRAITS_AMINO_ACIDS_CHOTHIA98) <- sort(unique(AMINO_ACIDS))
TRAITS_AMINO_ACIDS <- array(NA, 21)

#### Sysdata ####

# 5x312 logical matrix of CDR positions
# CDR_Nuc_Mat

# 5x312 logical matrix of FWR positions
# FWR_Nuc_Mat

# Vector of codon amino acid translations
# CODON_AA_TABLE

# 12x216 matrix of replacement and silent mutation permutations
# CODON_TABLE