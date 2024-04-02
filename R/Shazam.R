#' @keywords internal
#' @aliases shazam-package
"_PACKAGE"

# Shazam package documentation and import directives

#' The shazam package
#'
#' Dramatic improvements in high-throughput sequencing technologies now enable 
#' large-scale characterization of Ig repertoires, defined as the collection of transmembrane 
#' antigen-receptor proteins located on the surface of T and B lymphocytes. The \code{shazam}
#' package provides tools for advanced analysis of somatic hypermutation (SHM) in
#' immunoglobulin (Ig) sequences. The key functions in \code{shazam}, broken down topic, are 
#' described below.
#' 
#' @section  Mutational profiling:
#' \code{shazam} provides tools to quantify the extent and nature of SHM within
#' full length V(D)J sequences as well as sub-regions (eg, FWR and CDR).
#' Quantification of expected mutational loaded, under specific SHM targeting 
#' models, can also be performed along with model driven simulations of SHM.
#' 
#' \itemize{
#'   \item  \link{collapseClones}:           Build clonal consensus sequences.
#'   \item  \link{consensusSequence}:        Build a single consensus sequence.
#'   \item  \link{observedMutations}:        Compute observed mutation counts and frequencies.
#'   \item  \link{expectedMutations}:        Compute expected mutation frequencies.
#'   \item  \link{shmulateSeq}:              Simulate mutations in a single sequence.
#'   \item  \link{shmulateTree}:             Simulate mutations over a lineage tree.
#'   \item  \link{setRegionBoundaries}:      Extends a region definition to include CDR3 and FWR4.
#' }
#' 
#' @section  SHM targeting models:
#' Computational models and analyses of SHM have separated the process 
#' into two independent components: 
#' \enumerate{
#'   \item  A mutability model that defines where mutations occur.
#'   \item  A nucleotide substitution model that defines the resulting mutation.
#' }
#' Collectively these are what form the targeting model of SHM. \code{shazam} 
#' provides empirically derived targeting models for both humans and mice,
#' along with tools to build these mutability and substitution models from data.
#' 
#' \itemize{
#'   \item  \link{createTargetingModel}:     Build a 5-mer targeting model.
#'   \item  \link{plotMutability}:           Plot 5-mer mutability rates.
#'   \item  \link{HH_S5F}:                   Human 5-mer SHM targeting model.
#'   \item  \link{MK_RS5NF}:                 Mouse 5-mer SHM targeting model.
#' }
#'
#' @section  Quantification of selection pressure:
#' Bayesian Estimation of Antigen-driven Selection in Ig Sequences is a 
#' novel method for quantifying antigen-driven selection in high-throughput
#' Ig sequence data. Targeting models created using \code{shazam} can be used 
#' to estimate the null distribution of expected mutation frequencies used
#' by BASELINe, providing measures of selection pressure informed by known 
#' AID targeting biases.
#'              
#' \itemize{
#'   \item  \link{calcBaseline}:             Calculate the BASELINe probability
#'                                           density functions (PDFs).
#'   \item  \link{groupBaseline}:            Combine PDFs from sequences grouped
#'                                           by biological or experimental relevance.
#'   \item  \link{summarizeBaseline}:        Compute summary statistics from BASELINe PDFs.
#'   \item  \link{testBaseline}:             Perform significance testing for the difference
#'                                           between BASELINe PDFs.
#'   \item  \link{plotBaselineDensity}:      Plot the probability density functions
#'                                           resulting from selection analysis.
#'   \item  \link{plotBaselineSummary}:      Plot summary statistics resulting from 
#'                                           selection analysis.
#' }
#'
#' @section  Mutational distance calculation:
#' \code{shazam} provides tools to compute evolutionary distances between 
#' sequences or groups of sequences, which can leverage SHM targeting 
#' models. This information is particularly useful in understanding and 
#' defining clonal relationships.
#'              
#' \itemize{
#'   \item  \link{findThreshold}:            Identify clonal assignment threshold based on 
#'                                           distances to nearest neighbors.
#'   \item  \link{distToNearest}:            Tune clonal assignment thresholds by calculating 
#'                                           distances to nearest neighbors.
#'   \item  \link{calcTargetingDistance}:    Construct a nucleotide distance matrix from a 
#'                                           5-mer targeting model.
#' }
#'
#' @name     shazam
#' @docType  _PACKAGE
#' @references
#' \enumerate{
#'   \item  Hershberg U, et al. Improved methods for detecting selection by mutation 
#'            analysis of Ig V region sequences. 
#'            Int Immunol. 2008 20(5):683-94.
#'   \item  Uduman M, et al. Detecting selection in immunoglobulin sequences. 
#'            Nucleic Acids Res. 2011 39(Web Server issue):W499-504. (Corrections at 
#'            http://selection.med.yale.edu/baseline/correction/) 
#'   \item  Yaari G, et al. Quantifying selection in high-throughput immunoglobulin 
#'            sequencing data sets. 
#'            Nucleic Acids Res. 2012 40(17):e134.
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4:358.
#'   \item  Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O'Connor K,
#'          Vigneault F, Shlomchik M and Kleinstein S (2016). A Model of Somatic Hypermutation 
#'          Targeting in Mice Based on High-Throughput Ig Sequencing Data. The Journal of 
#'          Immunology, 197(9), 3566-3574.
#'  }
#' 
#' @import   ggplot2
#' @import   graphics
#' @import   methods
#' @import   utils
#' @importFrom  alakazam    getAllele getGene getFamily getSegment getLocus groupGenes
#'                          getAAMatrix getDNAMatrix IUPAC_DNA
#'                          pairwiseDist nonsquareDist pairwiseEqual 
#'                          seqDist seqEqual
#'                          isValidAASeq translateStrings gridPlot
#'                          getMRCA getPathLengths tableEdges
#'                          progressBar baseTheme checkColumns cpuCount
#'                          makeChangeoClone summarizeSubtrees buildPhylipLineage
#' @importFrom  ape         mst
#' @importFrom  diptest     dip.test
#' @importFrom  doParallel  registerDoParallel
#' @importFrom  dplyr       do n desc %>%
#'                          bind_cols bind_rows combine
#'                          filter select arrange 
#'                          group_by ungroup group_indices
#'                          mutate summarize
#'                          mutate_at summarize_at
#'                          rename transmute
#'                          left_join recode pull
#' @importFrom  foreach     foreach %dopar% registerDoSEQ
#' @importFrom  igraph      V E as_adjacency_matrix graph_from_data_frame
#'                          vertex_attr set_vertex_attr 
#'                          layout_as_tree V<-
#' @importFrom  iterators   icount
#' @importFrom  KernSmooth  bkde
#' @importFrom  lazyeval    interp
#' @importFrom  MASS        fitdistr
#' @importFrom  progress    progress_bar
#' @importFrom  rlang       sym syms .data
#' @importFrom  scales      log2_trans log10_trans trans_breaks trans_format
#'                          math_format percent scientific pretty_breaks
#' @importFrom  seqinr      c2s s2c words translate
#' @importFrom  stats       na.omit setNames ecdf sd cor cov median mad
#'                          approx convolve weighted.mean p.adjust
#'                          dbeta pbeta qbeta rbeta optim optimize
#'                          dnorm pnorm runif dgamma pgamma uniroot na.exclude
#'                          as.dist cutree integrate
#' @importFrom  stringi     stri_dup stri_flatten stri_join stri_length 
#'                          stri_sub stri_sub_replace stri_detect_regex
#'                          stri_count_boundaries stri_count_regex 
#'                          stri_extract_all_regex stri_extract_first_regex  
#'                          stri_replace_all_regex stri_replace_first_regex
#' @importFrom  tidyr       gather spread pivot_wider
#' @importFrom  tidyselect  all_of any_of
NULL

# Package loading actions
.onAttach <- function(libname, pkgname) {
    msg <- citation(pkgname)
    msg <-paste(c(format(msg,"citation")),collapse="\n\n")
    packageStartupMessage(msg)
}

#### Sysdata ####

# Deprecated (v0.1.4) mouse single nucleotide distance matrix
#
# Single nucleotide distance matrix of somatic hypermutation targeting based on 
# Mus musculus Ig sequence data.
#
# @format   A symmetric matrix of nucleotide substitution distance scores with 
#           row names and column names definition the specific substitution.
# 
# @references
# \enumerate{
#   \item  Smith DS, et al. Di- and trinucleotide target preferences of somatic 
#            mutagenesis in normal and autoreactive B cells. 
#            J Immunol. 1996 156:2642-52. 
# }
#
# M1N_Compat


# Deprecated (v0.1.4) Human single nucleotide distance matrix.
#
# Single nucleotide distance matrix of somatic hypermutation targeting based on 
# human Ig sequence data.
#
# @format   A symmetric matrix of nucleotide substitution distance scores with 
#           row names and column names definition the specific substitution.
# 
# @references
# \enumerate{
#   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#            Front Immunol. 2013 4(November):358.
# }
#
# HS1F_Compat

# Ordered nucleotide character set
# NUCLEOTIDES <- c("A", "C", "G", "T", "N", "-", ".")

# IMGT V segment length
# VLENGTH <- 312

# 5x312 logical matrix of CDR positions
# CDR_Nuc_Mat

# 5x312 logical matrix of FWR positions
# FWR_Nuc_Mat

# 12x216 matrix of replacement and silent mutation permutations
# CODON_TABLE

# 1x24 vector of amino acid charge classes
# AMINO_ACIDS_CHARGE

# 1x24 vector of amino acid hydropathy classes
# AMINO_ACIDS_HYDROPATHY

# 1x24 vector of amino acid polarity classes
# AMINO_ACIDS_POLARITY

# TODO: What is this?
# CONST_I

# TODO: And what is this?
# BAYESIAN_FITTED

# IMGT-KABAT numbering mapping
# As described hare http://www.imgt.org/IMGTScientificChart/Numbering/IMGT-Kabat_part1.html
# CONVERT_NUM_REF

# Add built-in variables to global variables environment
utils::globalVariables(c("HH_S1F", "HKL_S1F", "MK_RS1NF",
                         "HH_S5F", "HKL_S5F", "MK_RS5NF", "U5N",
                         "IMGT_V_BY_REGIONS"), package="shazam")