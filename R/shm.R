#' shm
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
#' @section  Mutability profiling:
#' \itemize{
#'   \item  \code{\link{createMutabilityModel}}:    Builds a targeting model.
#'   \item  \code{\link{createSubstitutionModel}}:  Builds a mutability model.
#'   \item  \code{\link{plotHedgehog}}:             Plots a mutability model.
#' }
#'
#' @section  Distance profiling:
#' \itemize{
#'   \item  \code{\link{distToNearest}}:            Calculate distances to nearest-neighbors.
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
#' @name     shm
#' @docType  package
#'
#' @import   alakazam
#' @import   doSNOW
#' @import   ggplot2
#' @import   grid
#' @import   plyr
#' @import   seqinr
#' @import   SDMTools
#' @import   zoo
NULL
