# Project documentation for shm
# 
# @author     Mohamed Uduman, Gur Yaari, Namita Gupta, Jason Anthony Vander Heiden
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.1
# @date       2015.03.07


#' The shm package
#'
#' Provides tools for advanced anaylisis of Immunoglobulin (Ig) Somatic HyperMutation 
#' (SHM), including BASELINe, a novel method for quantifying antigen-driven selection in 
#' high-throughput Ig sequencing data.
#' 
#' Dramatic improvements in high-throughput sequencing technologies now enable 
#' large-scale characterization of Ig repertoires, defined as the collection of transmembrane 
#' antigen-receptor proteins located on the surface of T and B lymphocytes.
#' 
#' The \code{shm} package provides tools for advanced analysis of Ig sequences following 
#' germline segment assignment. Namely, the analysis of SHM. 
#' Which includes:
#'  \itemize{
#'      \item   Statistical analysis of SHM patterns \cr
#'              Computational models and analyses of SHM have separated the process 
#'              into two independent components: 
#'              \enumerate{
#'                  \item  A mutability model that defines where mutations occur.
#'                  \item  A nucleotide substitution model that defines the resulting mutation.
#'              }
#'              Collectively these are what form the targeting model of SHM. \code{shm} 
#'              provides tools to build these mutability and substitution (i.e. targeting) 
#'              models.
#'                  
#'      \item   BASELINe \cr
#'              Bayesian Estimation of Antigen-driven Selection in Ig Sequences is a 
#'              novel method for quantifying antigen-driven selection in high-throughput
#'              Ig sequence data. The targeting model created using \code{shm} is used 
#'              to estimate the null distribution of expected mutation frequencies in 
#'              BASELINe.
#'              
#'      \item   Distance calculations \cr
#'              Based on the underlying SHM targeting (calculated using \code{shm}) one 
#'              can compute evolutionary distances between sequences or groups of 
#'              sequences. This information is particularly useful in understanding and 
#'              defining clonal relationships.
#'  }
#' 
#' Below are the functions in \code{shm} broken down by the three main tasks described
#' above:
#' 
#' @section  Targeting models:
#' \itemize{
#'   \item  \link{createTargetingModel}:     Build a 5-mer targeting model.
#'   \item  \link{plotMutability}:           Plot 5-mer mutability rates.
#' }
#' 
#' @section  Mutational profiling:
#' \itemize{
#'   \item  \link{collapseByClone}:    Build clonal consensus sequence.
#'   \item  \link{calcDBObservedMutations}:  Compute observed mutation counts.
#'   \item  \link{calcDBExpectedMutations}:  Compute expected mutation frequencies.
#' }
#'
#' @section  Selection analysis:
#' \itemize{
#'   \item  \link{calcBaseline}:             Calculate the BASELINe probability
#'                                           density functions (PDFs).
#'   \item  \link{groupBaseline}:            Combine PDFs from sequences grouped
#'                                           by biological or experimental relevance.
#'   \item  \link{summarizeBaseline}:        Compute summary statistics from BASELINe PDFs.
#'   \item  \link{plotBaselineDensity}:      Plot the probability density functions
#'                                           resulting from selection analysis.
#'   \item  \link{plotBaselineSummary}:      Plot summary stastistics resulting from 
#'                                           selection analysis.
#' }
#'
#' @section  Distance profiling:
#' \itemize{
#'   \item  \link{distToNearest}:            Tune clonal assignment thresholds by calculating 
#'                                           distances to nearest-neighbors.
#'   \item  \link{calcTargetingDistance}:    Construct a nucleotide distance matrix from a 
#'                                           5-mer targeting model.
#' }
#'
#' @references
#' \enumerate{
#'   \item  Hershberg U, et al. Improved methods for detecting selection by mutation 
#'            analysis of Ig V region sequences. 
#'            Int Immunol. 2008 20(5):683-94.
#'   \item  Uduman M, et al. Detecting selection in immunoglobulin sequences. 
#'            Nucleic Acids Res. 2011 39(Web Server issue):W499-504.
#'   \item  Yaari G, et al. Quantifying selection in high-throughput immunoglobulin 
#'            sequencing data sets. 
#'            Nucleic Acids Res. 2012 40(17):e134.
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4:358.
#'  }
#'
#' @seealso
#' The Change-O suite of tools includes three separate R packages: \link[alakazam]{alakazam}, 
#' \link[tigger]{tigger}, and \link[shm]{shm}.
#' 
#' @name     shm
#' @docType  package
#' @import   alakazam
#' @import   data.table
#' @import   doParallel
#' @import   foreach
#' @import   grid
#' @import   ggplot2
#' @import   methods
#' @import   plyr
#' @import   reshape2
#' @import   scales
#' @importFrom  iterators  icount
#' @importFrom  SDMTools   wt.sd
#' @importFrom  seqinr     c2s
#' @importFrom  seqinr     s2c
#' @importFrom  seqinr     words
#' @importFrom  seqinr     translate
NULL


#### Constants ####

# IMGT V-region definitions
#VREGIONS <- factor(c(rep("FWR", 78), 
#                     rep("CDR", 36), 
#                     rep("FWR", 51), 
#                     rep("CDR", 30), 
#                     rep("FWR", 117)), 
#                   levels=c("FWR", "CDR"))

# IMGT V-region length
#readEnd <- 312
#VLENGTH <- length(VREGIONS)
VLENGTH <- 312

# Nucleotide characters
NUCLEOTIDES <- c("A", "C", "G", "T", "N", "-", ".")

# Amino acid characters
#AMINO_ACIDS <- c("F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y", "*", "*", "C", "C", "*", "W", "L", "L", "L", "L", "P", "P", "P", "P", "H", "H", "Q", "Q", "R", "R", "R", "R", "I", "I", "I", "M", "T", "T", "T", "T", "N", "N", "K", "K", "S", "S", "R", "R", "V", "V", "V", "V", "A", "A", "A", "A", "D", "D", "E", "E", "G", "G", "G", "G")
#names(AMINO_ACIDS) <- c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")

#Amino Acid Traits
#"*" "A" "C" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y"
#B = "Hydrophobic/Burried"  N = "Intermediate/Neutral"  S="Hydrophilic/Surface")
#TRAITS_AMINO_ACIDS_CHOTHIA98 <- c("*","N","B","S","S","B","N","N","B","S","B","B","S","N","S","S","N","N","B","B","N")
#names(TRAITS_AMINO_ACIDS_CHOTHIA98) <- sort(unique(AMINO_ACIDS))
#TRAITS_AMINO_ACIDS <- array(NA, 21)

#### Sysdata ####

# 5x312 logical matrix of CDR positions
# CDR_Nuc_Mat

# 5x312 logical matrix of FWR positions
# FWR_Nuc_Mat

# Vector of codon amino acid translations
# CODON_AA_TABLE

# 12x216 matrix of replacement and silent mutation permutations
# CODON_TABLE