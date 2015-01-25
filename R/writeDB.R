#' Writes DB to file
#'
#' Write the contents of DB file to a tab delimited text file.
#' 
#' @param   db  a data.frame of the DB file.
#' @param   fileName  The filename including path to write to
#' @export
writeDB <- function(db, fileName)  {
  write.table(db, file=fileName, quote=FALSE, sep="\t", row.names=FALSE)
}
