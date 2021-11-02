#### Convert Numbering ####

#' convertNumbering: IMGT-Kabat number conversion
#' 
#' Converts numbering systems like Kabat or IMGT using these conventions:
#' http://www.imgt.org/IMGTScientificChart/Numbering/IMGT-Kabat_part1.html
#' with Gaps (unoccupied positions) shown by "G" and Asterisks (*) shown by "S": 
#' arbitrary mappings (multiple possible "to" values) represented with "NA"
#'
#' @param   locus     string indicating heavy ("IGH") or light chains ("IGK" or "IGL)
#' @param   from      string indicating numbering system to convert to ("IMGT" or "KABAT")
#' @param   to        string indicating original numbering system ("IMGT" or "KABAT")
#' @param   calls     vector of strings representing original numbering
#' @return  A vector of string indicating the corresponding numbering
#' 
#' @examples
#' convertNumbering("IGH", "IMGT", "KABAT", c("51", "23", "110"))
#' convertNumbering("IGH", "KABAT", "IMGT", c("51", "23", "G"))
#' @export
convertNumbering <- function(locus, from, to, calls) {
    # Generate mapping from references
    from_map <- pull(CONVERT_NUM_REF, paste(locus, from, sep='_')) 
    to_map <- pull(CONVERT_NUM_REF, paste(locus, to, sep='_')) 
    
    # Check for compatibility of reference with input calls
    if (!all(calls %in% from_map)) {
        stop(paste("Formatting of following characters does not match reference: ", 
            toString(calls[!calls %in% from_map])))
    }
    
    out_calls <- dplyr::recode(as.character(calls), 
                                     !!! setNames(to_map, from_map))
    
    # Separate check for ambiguous values to convert to NA, not arbitrary call
    from_counts <- table(from_map)
    
    # Generate output
    for(i in calls){
        if (from_counts[i] > 1) {
            out_calls[which(calls==i)] <- "NA"
        }
    }
    
    return(as.character(out_calls))
}
