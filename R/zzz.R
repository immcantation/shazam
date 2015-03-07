.onLoad <- function(libname, pkgname) {
  options(warn=-1)
  #REGIONS <<- factor( c( rep("FWR",78), rep("CDR",36), rep("FWR",51), rep("CDR",30), rep("FWR",117)), levels=c("FWR","CDR") )
  #LENGTH_REGIONS <<- length(REGIONS)
  #testID <<- 1
  modelName = "HS5F_Targeting.RData"
  load(system.file("extdata", modelName , package="shm"))
  substitution <<- Targeting[["Substitution"]]
  mutability <<- Targeting[["Mutability"]]
  #NUCLEOTIDES_FAC <<- factor( 1:5, labels=c("A","C","G","T", "N") )
  #readEnd <<- length(REGIONS)
  #readEnd <<- 312
}
