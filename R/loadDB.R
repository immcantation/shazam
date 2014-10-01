#' Load DB file
#'
#' This function loads the basic (CLIP) DB file.\cr
#' This file should be tab delimited and contain at least the following three columns: \cr
#' 1. SEQUENCE_ID \cr
#' 2. SEQUENCE_GAP \cr
#' 3. GERMLINE_GAP (or GERMLINE_GAP_D_MASK) \cr
#'
#' Note that the sequences need to gapped according to the IMGT numbering scheme.
#'
#' @param   input_file   Full path to tab-delimitted input file
#' @return  a data.frame of the DB file
#' @export
loadDB <- function(input_file) {
  # Read file
  db <- read.delim(input_file, colClasses='character', fill=F,
                   na.strings=c('', 'NA', 'None', 'No Results'))
  # Fix column data type
  col_names <- names(db)
  if ('CLONE' %in% col_names) { db$CLONE <- as.numeric(db$CLONE) }
  return(db)
}
