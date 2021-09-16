#### Convert Numbering ####

# Converts numbering systems like Kabat or IMGT using these conventions:
# http://www.imgt.org/IMGTScientificChart/Numbering/IMGT-Kabat_part1.html
# with Gaps (unoccupied positions) shown by "G" and Asterisks (*) shown by "S": 
# arbitrary mappings (multiple possible "to" values) represented with "NA"
#
# @param   locus   string indicating heavy ("IGH") or light chains ("IGK" or "IGL)
# @param   from    string indicating numbering system to convert to ("IMGT" or "KABAT")
# @param   to      string indicating original numbering system ("IMGT" or "KABAT")
# @param   x       vector of strings representing original numbering
# @return  A vector of string indicating the corresponding numbering
#
# @examples
# convertNumbering("IGH", "IMGT", "KABAT", c("51", "23", "110"))
# convertNumbering("IGH", "KABAT", "IMGT", c("51", "23", "G"))

library(dplyr)

convertNumbering <- function(locus, from, to, arrayOfSeqs) {
  from_map <- pull(CONVERT_NUM_REF, paste(locus, from, sep='_')) 
  to_map <- pull(CONVERT_NUM_REF, paste(locus, to, sep='_')) 

  # Recode has issues with alternative types
  if (!all(arrayOfSeqs %in% from_map)) {
    stop(paste("Formatting of following characters does not match reference: ", 
               toString(arrayOfSeqs[!arrayOfSeqs %in% from_map])))
  }
  
  out_arrayOfSeqs <- dplyr::recode(as.character(arrayOfSeqs), !!! setNames(to_map, from_map))
  
  # Separate check for ambiguous values to convert to NA, not arbitrary
  from_counts <- table(from_map)
  
  for(i in arrayOfSeqs){
    if (from_counts[i] > 1) {
      out_arrayOfSeqs[which(arrayOfSeqs==i)[1]] <- "NA"
    }
  }
  
  return(as.character(out_arrayOfSeqs))
}
