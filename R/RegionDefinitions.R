# S4 Class defining RegionDefinitions
# @author     Mohamed Uduman, Gur Yaari
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.1
# @date       2015.03.08

#' @include shm.R
NULL

#### Classes ####

#' S4 class defining a region definition
#' 
#' \code{RegionDefinition} defines a common data structure for defining the regions
#' (boundaries) of the Ig sequence.
#' 
#' @slot    name            Name of the RegionDefinition
#' @slot    description     Description of the model and its source
#' @slot    boundaries      \code{factor} defining the regions (boundaries) of the 
#'                              sequence. The levels and values of the \code{factor} 
#'                              determine the nubmer of regions (e.g. CDR & FWR).
#' @slot    seqLength      the length of the sequence                        
#' @slot    regions         the levels of the boundaries (e.g. CDR & FWR)
#' @slot    labels          the labels for the boundary/mutations combinations
#'                          e.g. CDR_R CDR_S FWR_R, FWR_S.
#' @slot    citation       Publication source.
#'    
#' @name RegionDefinition
#' @export
setClass("RegionDefinition", 
         slots=c(name="character",
                 description="character",
                 boundaries="factor",
                 seqLength="numeric",
                 regions="character",
                 labels="character",
                 citation="character"),
         prototype=c(name="IMGT_V_NO_CDR3",
                     description="IMGT_Numbering scheme defining the V gene up till but not including CDR3",
                     boundaries=factor( c( rep("FWR", 78), 
                                           rep("CDR", 36),  
                                           rep("FWR", 51), 
                                           rep("CDR", 30), 
                                           rep("FWR", 117) ),
                                        levels = c("CDR","FWR")),
                     seqLength=312,
                     regions=c("CDR","FWR"),
                     labels=c("CDR_R", "CDR_S","FWR_R","FWR_S"),
                     citation="Lefranc MP et al. (2003)"))

#### RegionDefinition building functions #####

#' Creates a RegionDefinition
#' 
#' \code{createRegionDefinition} creates a \code{RegionDefinition}.
#'
#' @param    name           name of the region definition.
#' @param    boundaries     \code{factor} defining the regions (boundaries) of the sequence.
#'                          the levels and values of the \code{factor} determine the 
#'                          nubmer of regions (e.g. CDR & FWR).
#' @param    description    description of the region definition and its source data.
#' @param    citation       publication source.
#' 
#' @return   A \code{RegionDefinition} object.
#' 
#' @seealso  See \code{\link{RegionDefinition}} for the return object.
#' 
#' @examples
#' library(shm)
#' 
#' 
#' @export
createRegionDefinition <- function(name=NULL,
                                   boundaries=NULL,
                                   description=NULL,
                                   citation=NULL
                                   ) {
    #Extract information from 'boundaries'
    # Determine the number of levels (e.g. CDR, FWR)
    regions <- levels(boundaries)
    # Determine the length of the boundaries
    seqLength <- length(boundaries)
    
    # Determine the combinations of levels_regionDefinition and R/S
    # e.g. CDR_R CDR_S FWR_R, FWR_S
    labels <- paste( rep(regions, each=2), 
                           rep(c("R", "S"), length(regions)), 
                           sep="_" )
    
    # Define RegionDefinition object
    regionDefinition <- new("RegionDefinition",
                            name=name,
                            description=description,
                            boundaries=boundaries,
                            seqLength=seqLength,
                            regions=regions,
                            labels=labels,
                            citation=citation)
    
    return(regionDefinition)
}

#### Data ####

#' IMGT unique numbering for V-DOMAIN
#'
#' Defines the CDR and FWR, according to the IMGT unique numbering scheme, for V segments
#' uptill the begining of (but not including CDR3)
#'
#' @format \code{\link{RegionDefinition}} object.
#' 
#' @family IMGT unique numbering schemes
#' @references
#' \enumerate{
#'   \item  Lefranc MP, Pommié C, Ruiz M, Giudicelli V, Foulquier E, Truong L, 
#'              Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin 
#'              and T cell receptor variable domains and Ig superfamily V-like domains. 
#'              Developmental and comparative immunology. 2003;27:55–77.
#' }
"IMGT_V_NO_CDR3"




#' IMGT unique numbering for V-DOMAIN
#'
#' Defines the CDR and FWR, according to the IMGT unique numbering scheme, for V segments
#' (including) CDR3.
#'
#' @format \code{\link{RegionDefinition}} object.
#' 
#' @family IMGT unique numbering schemes
#' 
#' @references
#' \enumerate{
#'   \item  Lefranc MP, Pommié C, Ruiz M, Giudicelli V, Foulquier E, Truong L, 
#'              Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin 
#'              and T cell receptor variable domains and Ig superfamily V-like domains. 
#'              Developmental and comparative immunology. 2003;27:55–77.
#' }
#'
"IMGT_V"



#' IMGT unique numbering for V-DOMAIN.
#'
#' Defines the indiviudal CDRs and FWRs, according to the IMGT unique numbering scheme, 
#' for V segments uptill the begining of (but not including CDR3).
#'
#' @format \code{\link{RegionDefinition}} object.
#' 
#' @family IMGT unique numbering schemes
#' 
#' @references
#' \enumerate{
#'   \item  Lefranc MP, Pommié C, Ruiz M, Giudicelli V, Foulquier E, Truong L, 
#'              Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin 
#'              and T cell receptor variable domains and Ig superfamily V-like domains. 
#'              Developmental and comparative immunology. 2003;27:55–77.
#' }
"IMGT_V_BY_REGIONS_NO_CDR3"


#' IMGT unique numbering for V-DOMAIN.
#'
#' Defines the indiviudal CDRs and FWRs, according to the IMGT unique numbering scheme, 
#' for V segments (including CDR3).
#'
#' @format \code{\link{RegionDefinition}} object.
#' 
#' @family IMGT unique numbering schemes
#' 
#' @references
#' \enumerate{
#'   \item  Lefranc MP, Pommié C, Ruiz M, Giudicelli V, Foulquier E, Truong L, 
#'              Thouvenin-Contet V, Lefranc G. IMGT unique numbering for immunoglobulin 
#'              and T cell receptor variable domains and Ig superfamily V-like domains. 
#'              Developmental and comparative immunology. 2003;27:55–77.
#' }
"IMGT_V_BY_REGIONS"
