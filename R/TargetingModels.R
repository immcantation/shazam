# Targeting models
# 
# @author     Gur Yaari, Mohamed Uduman, Jason Anthony Vander Heiden
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.1
# @date       2015.04.07

#' @include shm.R
NULL

#### Data ####

#' 3-mer targeting model.
#'
#' 3-mer model of somatic hypermutation targeting based on Mus musculus Ig sequence data.
#'
#' @format \code{\link{TargetingModel}} object.
#' 
#' @references
#' \enumerate{
#'   \item  Smith DS, et al. Di- and trinucleotide target preferences of somatic 
#'            mutagenesis in normal and autoreactive B cells. 
#'            J Immunol. 1996 156:2642–52. 
#' }
#'
#' @seealso  See \code{\link{HS5FModel}} for the 5-mer model.
"M3NModel"


#' 5-mer targeting model.
#'
#' 5-mer model of somatic hypermutation targeting based on analysis of silent mutations
#' in functional Ig sequences from Homo sapiens.
#'
#' @format \code{\link{TargetingModel}} object.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#'  
#' @seealso  See \code{\link{M3NModel}} for the 3-mer model.
"HS5FModel"


#### Classes ####

#' S4 class defining a targeting model
#' 
#' \code{TargetingModel} defines a common data structure for mutability, substitution and
#' targeting of immunoglobulin (Ig) sequencing data in a 5-mer microsequence context.
#' 
#' @slot     name          Name of the model.
#' @slot     description   Description of the model and its source data.
#' @slot     species       Genus and species of the source sequencing data.
#' @slot     date          Date the model was built.
#' @slot     citation      Publication source.
#' @slot     substitution  Normalized probability of the center nucleotide of a given 5-mer 
#'                         mutating to a different nucleotide. The substitution model 
#'                         is stored as a 5x3125 matrix of rates. Rows define
#'                         the mutated nucleotide at the center of each 5-mer, one of 
#'                         c(A, C, G, T, N), and columns define the complete 5-mer 
#'                         of the unmutated nucleotide sequence.
#' @slot     mutability    Normalized probability of a given 5-mer being mutated. The 
#'                         mutability model is stored as a numeric vector of length 3125 
#'                         with mutability rates for each 5-mer.
#' @slot     targeting     Probability of a given mutation ocurring, defined as 
#'                         \eqn{mutability * substitution}. The targeting model 
#'                         is stored as a 5x3125 matrix of rates. Rows define
#'                         the mutated nucleotide at the center of each 5-mer, one of 
#'                         c(A, C, G, T, N), and columns define the complete 5-mer 
#'                         of the unmutated nucleotide sequence.
#' 
#' @seealso  See \code{\link{createMutabilityModel}} and \code{\link{createSubstitutionModel}} 
#'           for building models from sequencing data.
#'           
#' @name TargetingModel
#' @export
setClass("TargetingModel", 
         slots=c(name="character",
                 description="character",
                 species="character",
                 date="character",
                 citation="character",
                 mutability="numeric",
                 substitution="matrix",
                 targeting="matrix"),
         prototype=c(name="name",
                     description="description",
                     species="species",
                     date="2000-01-01",
                     citation="citation",
                     mutability=numeric(3125),
                     substitution=matrix(0, 5, 3125),
                     targeting=matrix(0, 5, 3125)))

#### Model building functions #####

#' Builds a substitution model
#'
#' \code{createSubstitutionModel} builds a 5-mer nucleotide substitution model by counting 
#' the number of substitution mutations occuring in the center position for all 5-mer 
#' motifs.
#'
#' @param    db                data.frame containing sequence data.
#' @param    model             type of model to create. The default model, "RS", creates 
#'                             a model by counting both replacement and silent mutations.
#'                             The "S" specification builds a model by counting only 
#'                             silent mutations.
#' @param    sequenceColumn    name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn    name of the column containing IMGT-gapped germline sequences.
#' @param    vCallColumn       name of the column containing the V-segment allele call.
#' @param    multipleMutation  string specifying how to handle multiple mutations occuring 
#'                             within the same 5-mer. If \code{"independent"} then multiple 
#'                             mutations within the same 5-mer are counted indepedently. 
#'                             If \code{"ignore"} then 5-mers with multiple mutations are 
#'                             excluded from the total mutation tally.
#' 
#' @return   A 4x1024 matrix of substitution counts for each 5-mer motif with rows names 
#'           defining the center nucleotide, one of c(A, C, G, T), and column names 
#'           defining the 5-mer nucleotide sequence.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#'
#' @seealso  See \code{\link{createMutabilityModel}} for building a mutability model.
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' db <- readChangeoDb(file)
#' 
#' # Create model using all mutations
#' sub_model <- createSubstitutionModel(db)
#' 
#' # Create model using only silent mutations
#' sub_model <- createSubstitutionModel(db, model="S")
#' 
#' @export
createSubstitutionModel <- function(db, model=c("RS", "S"), sequenceColumn="SEQUENCE_IMGT",
                                    germlineColumn="GERMLINE_IMGT_D_MASK",
                                    vCallColumn="V_CALL",
                                    multipleMutation=c("independent", "ignore"))  {
    # Evaluate argument choices
    model <- match.arg(model)
    multipleMutation <- match.arg(multipleMutation)
    
    # Make sure the columns specified exist
    if (!(germlineColumn %in% names(db))) {
        stop("Please specify the germline column name")
    } 
    if (!(sequenceColumn %in% names(db))) {
        stop("Please specify the sequence column name")
    } 
    if (!(vCallColumn %in% names(db))) {
        stop("Please specify the V call column name")
    }   
    
    # Setup
    nuc_chars <- NUCLEOTIDES[1:4]
    nuc_words <- words(4, nuc_chars)
    vh_families <- paste(rep("VH", 7), 1:7, sep="")
    
    # Define empty return list of lists
    substitutionMatrix <- matrix(0, ncol=4, nrow=4, dimnames=list(nuc_chars, nuc_chars))
    substitutionList <- list()    
    for(vh_fam in vh_families){
        substitutionList[[vh_fam]] = list()
        for(word in nuc_words){
            substitutionList[[vh_fam]][[word]] = substitutionMatrix
        }
    }
    
    
    # Redefine vh_families to only those found in the data
    vh_families <- substr(db[, vCallColumn],
                          gregexpr("IGH", db[, vCallColumn])[[1]][1] + 4,
                          gregexpr("IGH", db[, vCallColumn])[[1]][1] + 4)
    vh_families <- paste("VH", vh_families, sep="")
    
    mutations <- listObservedMutations(db, sequenceColumn=sequenceColumn, 
                                       germlineColumn=germlineColumn)
    
    
    if (model == "S") { # Silent model
        for(index in 1:length(mutations)) {
            cSeq <-  s2c(db[index,sequenceColumn])
            cGL  <-  s2c(db[index,germlineColumn])
            indexMutation <- mutations[[index]]
            vh_fam <- vh_families[index]
            
            positions <- as.numeric(names(indexMutation))
            positions <- positions[positions<=VLENGTH]
            positions <- positions[!is.na(positions)]
            for( position in  positions){
                wrd <-  c2s(c(cGL[(position-2):(position-1)],cGL[(position+1):(position+2)]))
                codonNucs = getCodonPos(position)
                codonGL = cGL[codonNucs]
                codonSeq = cSeq[codonNucs]
                muCodonPos = {position-1}%%3+1
                seqAtMutation <- codonSeq[muCodonPos]
                glAtMutation <- codonGL[muCodonPos]
                if( !any(codonGL=="N") & !any(codonSeq=="N") ){
                    if(multipleMutation=="independent"){ # if independent mutations are included
                        codonSeq <- codonGL
                        codonSeq[muCodonPos] <- seqAtMutation
                    }
                    codonPermutate <- matrix(rep(codonGL,3),ncol=3,byrow=T)
                    codonPermutate[,muCodonPos] <- canMutateTo(glAtMutation)[-4]
                    codonPermutate <- apply(codonPermutate,1,paste,collapse="")
                    codonPermutate <- matrix( c( codonPermutate, rep(c2s(codonGL),3) ), ncol=2, byrow=F)
                    muType <- mutationTypeOptimized(codonPermutate)
                    if(!length(grep("N",wrd))){
                        if( sum(muType=="S") == length(muType) ){
                            substitutionList[[vh_fam]][[wrd]][glAtMutation,seqAtMutation] <- (substitutionList[[vh_fam]][[wrd]][glAtMutation,seqAtMutation] + 1)
                        }
                    }
                }
            }
        }
    } else if (model == "RS") { # RS model (All mutations)
        for (index in 1:length(mutations)) {
            cSeq <-  s2c(db[index,sequenceColumn])
            cGL  <-  s2c(db[index,germlineColumn])
            indexMutation <- mutations[[index]]
            vh_fam <- vh_families[index]
            
            positions <- as.numeric(names(indexMutation))
            positions <- positions[positions<=VLENGTH]
            positions <- positions[!is.na(positions)]
            for( position in  positions){
                wrd <-  c2s(c(cGL[(position-2):(position-1)],cGL[(position+1):(position+2)]))
                codonNucs = getCodonPos(position)
                codonGL = cGL[codonNucs]
                codonSeq = cSeq[codonNucs]
                muCodonPos = {position-1}%%3+1
                seqAtMutation <- codonSeq[muCodonPos]
                glAtMutation <- codonGL[muCodonPos]
                if( !any(codonGL=="N") & !any(codonSeq=="N") ){
                    if(multipleMutation=="independent"){ # if independent mutations are included
                        codonSeq <- codonGL
                        codonSeq[muCodonPos] <- seqAtMutation
                    }
                    if(!length(grep("N",wrd))){
                        substitutionList[[vh_fam]][[wrd]][glAtMutation,seqAtMutation] <- (substitutionList[[vh_fam]][[wrd]][glAtMutation,seqAtMutation] + 1)
                    }
                }
            }
        }
    }
    
    # TODO: this is redundant code. can probably move the whole listSubstitution definition to the top of the function.
    nuc_words = words(4, nuc_chars)
    vh_families <- paste(rep("VH", 7), 1:7, sep="")
    
    X1 <- rep(nuc_words, 7)
    X2 <- rep(vh_families, each=256)
    arrNames <- paste(X2, X1, sep="_")
    rm(X2, X1)
    
    listSubstitution <- array(0, dim=c(256*7,4,4), dimnames=list(arrNames, nuc_chars, nuc_chars))
    
    
    # substitutionList
    for(vh_fam in vh_families){
        listSubstitution[paste(vh_fam,nuc_words,sep="_"),,]<-t(sapply(nuc_words,function(word){substitutionList[[vh_fam]][[word]]}))
    }
    
    M<-list()
    NAMES<-sapply(dimnames(listSubstitution)[[1]],function(x)strsplit(x,"_",fixed=TRUE)[[1]])
    CTS<-list()
    X1<-list()
    
    CT_FiveMers <-
        c( paste(substring(words(3,nuc_chars),1,1), "C",substring(words(3,nuc_chars),2,3),sep=""),
           paste(substring(words(3,nuc_chars),1,1),"G",substring(words(3,nuc_chars),2,3),sep=""),
           paste(substring(words(3,nuc_chars),1,1),"T",substring(words(3,nuc_chars),2,3),sep="")
        )
    
    for(CT in CT_FiveMers){
        CTS[[CT]]<-nuc_words[sapply(nuc_words,function(x)substring(x,1,4))== CT]
        INDEx<-NAMES[2,]%in%CTS[[CT]]
        M[[CT]]<- t(sapply(1:4,function(i)apply(listSubstitution[INDEx,i,],2,sum)))
        rownames(M[[CT]]) <- nuc_chars
    }
    M[["E"]]<-M[[paste(substring(words(3,nuc_chars),1,1),"C",substring(words(3,nuc_chars),2,3),sep="")[1]]]
    for(CT in c(paste(substring(words(3,nuc_chars),1,1),"C",substring(words(3,nuc_chars),2,3),sep=""),paste(substring(words(3,nuc_chars),1,1),"G",substring(words(3,nuc_chars),2,3),sep=""),paste(substring(words(3,nuc_chars),1,1),"T",substring(words(3,nuc_chars),2,3),sep=""))[-1]){
        M[["E"]]<- M[["E"]]+M[[CT]]
    }
    rownames(M[["E"]]) <- nuc_chars
    
    
    # fivemer=M; FIVEMER="CCATT"
    .simplifivemer <- function(fivemer, FIVEMER, Thresh=20) {
        Nuc=substr(FIVEMER,3,3)
        Nei=paste(substr(FIVEMER,1,2),substr(FIVEMER,4,5),collapse="",sep="")
        
        if(sum(fivemer[[Nei]][Nuc,])>Thresh & sum(fivemer[[Nei]][Nuc,]==0)==1){
            return(fivemer[[Nei]][Nuc,]);
        }
        else{
            if(substring(Nei,2,2)!="A"){
                FIVE=fivemer[[Nei]][Nuc,]
                for(i in 1:4){
                    for(j in 1:4){
                        MutatedNeighbor=paste(nuc_chars[i],substring(Nei,2,3),nuc_chars[j],collapse="",sep="")
                        FIVE=FIVE+fivemer[[MutatedNeighbor]][Nuc,]
                    }
                }
                return(FIVE)
            }
            if(substring(Nei,2,2)=="A"){
                FIVE=fivemer[[paste(substring(Nei,1,1),canMutateTo(substring(Nei,2,2))[1],substring(Nei,3,4),collapse="",sep="")]][Nuc,]
                for(i in 2){
                    for(j in 2:3){
                        MutatedNeighbor=paste(substring(Nei,1,(i-1)),canMutateTo(substring(Nei,i,i))[j],substring(Nei,(i+1),4),collapse="",sep="")
                        FIVE=FIVE+fivemer[[MutatedNeighbor]][Nuc,]
                    }
                }
                if(sum(FIVE)>Thresh & sum(FIVE==0)==1){
                    return(FIVE)
                }
                else{
                    for(i in 1:3){
                        for(j in 1:3){
                            MutatedNeighbor=paste(canMutateTo(substring(Nei,1,1))[i],canMutateTo(substring(Nei,2,2))[j],substring(Nei,3,4),collapse="",sep="")
                            FIVE=FIVE+fivemer[[MutatedNeighbor]][Nuc,]
                            MutatedNeighbor=paste(substring(Nei,1,1),canMutateTo(substring(Nei,2,2))[j],canMutateTo(substring(Nei,3,3))[i],substring(Nei,4,4),collapse="",sep="")
                            FIVE=FIVE+fivemer[[MutatedNeighbor]][Nuc,]
                            MutatedNeighbor=paste(substring(Nei,1,1),canMutateTo(substring(Nei,2,2))[j],substring(Nei,3,3),canMutateTo(substring(Nei,4,4))[i],collapse="",sep="")
                            FIVE=FIVE+fivemer[[MutatedNeighbor]][Nuc,]
                        }
                    }
                }
                return(FIVE)
            }
        }
    }
    
    substitutionModel <- sapply(words(5, nuc_chars), function(x) .simplifivemer(M, x))

    return(substitutionModel)
}


#' Builds a mutability model
#'
#' \code{createMutabilityModel} builds a 5-mer nucleotide mutability model by counting 
#' the number of mutations occuring in the center position for all 5-mer motifs.
#'
#' @param    db                 data.frame containing sequence data.
#' @param    substitutionModel  matrix of 5-mers substitution counts built by 
#'                              \code{\link{createSubstitutionModel}}.
#' @param    model              type of model to create. The default model, "RS", creates 
#'                              a model by counting both replacement and silent mutations.
#'                              The "S" specification builds a model by counting only 
#'                              silent mutations.
#' @param    sequenceColumn     name of the column containing IMGT-gapped sample sequences.
#' @param    germlineColumn     name of the column containing IMGT-gapped germline sequences.
#' @param    vCallColumn        name of the column containing the V-segment allele call.
#' @param    multipleMutation   string specifying how to handle multiple mutations occuring 
#'                              within the same 5-mer. If \code{"independent"} then multiple 
#'                              mutations within the same 5-mer are counted indepedently. 
#'                              If \code{"ignore"} then 5-mers with multiple mutations are 
#'                              excluded from the total mutation tally.
#' 
#' @return  A list with one element for each row in the input \code{db} where each entry
#'          contains a list itself with:
#'          \itemize{
#'            \item \code{[[1]]} = numeric vector of 5-mer mutation counts normalized by 
#'                                 the total number of mutations in the sequence.
#'            \item \code{[[2]]} = count 5-mers that were mutated in the sequence.
#'            \item \code{[[3]]} = total number of mutations in the sequence.
#'          }
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution 
#'            based on synonymous mutations from high-throughput immunoglobulin sequencing 
#'            data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#' 
#' @seealso  See \code{\link{createSubstitutionModel}} for building the substitution model.
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' db <- readChangeoDb(file)
#' 
#' # Create model using all mutations
#' sub_model <- createSubstitutionModel(db)
#' mut_model <- createMutabilityModel(db, sub_model)
#'
#' # Create model using only silent mutations
#' sub_model <- createSubstitutionModel(db, model="S")
#' mut_model <- createMutabilityModel(db, sub_model, model="S")
#' 
#' @export
createMutabilityModel <- function(db, substitutionModel, model=c("RS", "S"),
                                  sequenceColumn="SEQUENCE_IMGT", 
                                  germlineColumn="GERMLINE_IMGT_D_MASK",
                                  vCallColumn="V_CALL",
                                  multipleMutation=c("independent", "ignore")) {
    # model="RS"; sequenceColumn="SEQUENCE_IMGT"; germlineColumn="GERMLINE_IMGT_D_MASK"; vCallColumn="V_CALL"; multipleMutation="independent"
    # Evaluate argument choices
    model <- match.arg(model)
    multipleMutation <- match.arg(multipleMutation)
    
    # Make sure the columns specified exist
    if (!(germlineColumn %in% names(db))) {
        stop("The germline column", germlineColumn, "was not found.")
    } 
    if (!(sequenceColumn %in% names(db))) {
        stop("The sequence column", sequenceColumn, "was not found.")
    } 
    if (!(vCallColumn %in% names(db))) {
        stop("The V-segment allele call column", vCallColumn, "was not found.")
    } 
    
    nuc_chars <- NUCLEOTIDES[1:4]
    mutations <- listObservedMutations(db, sequenceColumn=sequenceColumn, 
                                       germlineColumn=germlineColumn)
    
    substitutionModel <- apply(substitutionModel,2,function(x)x/sum(x))
    template <- rep(0, 1024)
    names(template) <- words(5,nuc_chars)
    
    COUNT<-list()
    for(index in 1:length(mutations)){
        COUNT[[index]]<-template
        indexMutation <- mutations[[index]]
        if(!sum(is.na(indexMutation))){
            cSeq <-  s2c(db[index,sequenceColumn])
            cGL  <-  s2c(db[index,germlineColumn])
            positions <- as.numeric(names(indexMutation))
            positions <- positions[positions<=VLENGTH]
            positions <- positions[!is.na(positions)]
            for( position in  positions){
                wrd5<-substr(db[index,germlineColumn],position-2,position+2)
                if( !grepl("[^ACGT]", wrd5) & nchar(wrd5)==5  ){
                    codonNucs = getCodonPos(position)
                    codonGL = cGL[codonNucs]
                    codonSeq = cSeq[codonNucs]
                    muCodonPos = {position-1}%%3+1
                    seqAtMutation <- codonSeq[muCodonPos]
                    glAtMutation <- codonGL[muCodonPos]
                    if( !any(codonGL%in%c("N","-", ".")) & !any(codonSeq%in%c("N","-", ".")) ){
                        if(multipleMutation=="independent"){ # if independent mutations are included
                            codonSeq <- codonGL
                            codonSeq[muCodonPos] <- seqAtMutation
                        }
                        if(!length(grep("N",wrd5))){
                            COUNT[[index]][wrd5]<- COUNT[[index]][wrd5]+1;
                        }
                    }
                }
            }
        }
    }
    
    BG_COUNT<-list()
    cat("Progress: 0%      50%     100%\n")
    cat("          ")
    pb <- txtProgressBar(min=1,max=length(mutations),width=20)
    for(index in 1:length(mutations)){
        setTxtProgressBar(pb, index)
        BG_COUNT[[index]]<-template
        sSeq <- gsub("\\.","",db[index,sequenceColumn])
        sGL <- gsub("\\.","",db[index,germlineColumn])
        cSeq <-  s2c(sSeq)
        cGL  <-  s2c(sGL)[1:VLENGTH]
        positions <- 3:(length(cGL)-2)
        for( position in  positions){
            wrd5<-substr(sGL,position-2,position+2)
            if( !grepl("[^ACGT]", wrd5) & nchar(wrd5)==5 ){
                codonNucs = getCodonPos(position)
                codonGL = cGL[codonNucs]
                codonSeq = cSeq[codonNucs]
                muCodonPos = {position-1}%%3+1
                seqAtMutation <- codonSeq[muCodonPos]
                glAtMutation <- codonGL[muCodonPos]
                if( !any(codonGL%in%c("N","-")) & !any(codonSeq%in%c("N","-")) ){
                    if(multipleMutation=="independent"){ # if independent mutations are included
                        codonSeq <- codonGL
                        codonSeq[muCodonPos] <- seqAtMutation
                    }
                    codonPermutate <- matrix(rep(codonGL,3),ncol=3,byrow=T)
                    codonPermutate[,muCodonPos] <- canMutateTo(glAtMutation)[-4]
                    codonPermutate <- apply(codonPermutate,1,paste,collapse="")
                    codonPermutate <- matrix( c( codonPermutate, rep(c2s(codonGL),3) ), ncol=2, byrow=F)
                    muType <- mutationTypeOptimized(codonPermutate)
                    if(!length(grep("N",wrd5)) & !length(grep("-",wrd5))){
                        for(m in 1:3){
                            if(muType[m]=="S"){
                                BG_COUNT[[index]][wrd5]<- BG_COUNT[[index]][wrd5] + substitutionModel[substr(codonPermutate[m,1],muCodonPos,muCodonPos),wrd5];
                            }
                        }
                    }
                }
            }
        }
        BG_COUNT[[index]][BG_COUNT[[index]]==0]<-NA
    }
    close(pb)
    cat("\n")
    
    Mutability<-list()
    for(i in 1:length(mutations)){
        Mutability[[i]]<-list()
        Mutability[[i]][[1]]<-COUNT[[i]]/BG_COUNT[[i]]
        Mutability[[i]][[1]]<-Mutability[[i]][[1]]/sum(Mutability[[i]][[1]],na.rm=TRUE)
        Mutability[[i]][[2]]<-length( mutations[[i]])
        Mutability[[i]][[3]]<-sum(COUNT[[i]],na.rm=TRUE)
        
    }
    
    # TODO:  what is this?
    # Aggregate mutability
    #MutabilityMatrix <- sapply(Mutability, function(x) x[[1]][1:1024])
    #MutabilityWeights <- sapply(Mutability, function(x) x[[2]][1:1024])
    #MutabilityMatrixNorm <- apply(MutabilityMatrix, 2, function(x) x/sum(x, na.rm=TRUE) * sum(!is.na(x)))
    #Mutability_Mean <- sapply(1:1024, function(i) weighted.mean(MutabilityMatrix[i,], MutabilityWeights[i,], na.rm=TRUE))
    #Mutability_SD<-sapply(1:1024, function(i) SDMTools::wt.sd(MutabilityMatrixNorm[i,], MutabilityWeights[i,]))
    
    return(Mutability)
}


#' Normalized substitution
#'
createSymmetricSubstitution <- function() {
    library("seqinr")
    library("dclone")
    library("gdata")
    
    load("FiveS_Substitution.RData")
    FiveS_Substitution <- rbind(FiveS_Substitution,"N"=rep(NA,ncol(FiveS_Substitution)))
    
    NUCS <- c("A", "C", "G", "T","N")
    
    
    withN_5mer <- words(5, alphabet=s2c("ACGTN"))
    withN_5mer_mat <- matrix(NA, nrow=5, ncol=length(withN_5mer), dimnames=list(NUCS,withN_5mer))
    for(mer in withN_5mer){
        #cat(mer,"\n")
        if( is.finite(match(mer,colnames(FiveS_Substitution))) ){
            withN_5mer_mat[,mer] <- FiveS_Substitution[,mer]
        }else{
            merAsChar <- s2c(mer)
            N_Index <- grep("[N]",merAsChar)
            if( any(N_Index==3) ){
                withN_5mer_mat[,mer] <- NA
            }else{
                merAsChar[N_Index] <- "."
                merAsStr <- c2s(merAsChar)
                merIndex <- grep(merAsStr,colnames(FiveS_Substitution))
                withN_5mer_mat[,mer] <- apply(FiveS_Substitution[,merIndex],1,mean,na.rm=TRUE)
            }
        }
    }
    S5F_Substitution <- withN_5mer_mat
    #save(S5F_Substitution,file="S5F_Substitution.RData")
    
    #Make Subs an array
    #S5F_Substitution_Array <- unmatrix(S5F_Substitution)
    
    #Mutability
    #load("FiveS_Mutability.RData")
    load("S5F_Mutability.RData")
    withN_5mer <- words(5, alphabet=s2c("ACGTN"))
    withN_5mer_mat <- array(NA, dim=length(withN_5mer), dimnames=list(withN_5mer))
    FiveS_Mutability <- as.vector(S5F_Mutability)
    FiveS_Mutability <- FiveS_Mutability/sum(FiveS_Mutability,na.rm=TRUE)
    for(mer in withN_5mer){
        #cat(mer,"\n")
        if( is.finite(match(mer,rownames(FiveS_Mutability))) ){
            withN_5mer_mat[mer] <- FiveS_Mutability[mer,]
        }else{
            merAsChar <- s2c(mer)
            N_Index <- grep("[N]",merAsChar)
            if( any(N_Index==3) ){
                withN_5mer_mat[mer] <- NA
            }else{
                merAsChar[N_Index] <- "."
                merAsStr <- c2s(merAsChar)
                merIndex <- grep(merAsStr,rownames(FiveS_Mutability))
                withN_5mer_mat[mer] <- mean(FiveS_Mutability[merIndex,],na.rm=TRUE)
            }
        }
    }
    
    
    S5F_Mutability <- withN_5mer_mat
    S5F_Mutability <- (S5F_Mutability/sum(S5F_Mutability,na.rm=TRUE)) * sum(!is.na(S5F_Mutability))
    #save(S5F_Mutability,file="S5F_Mutability.RData")
    
    
    
    Targeting <- list()
    Targeting[["Name"]] <- "HS5F"
    Targeting[["Species"]] <- "Human"
    Targeting[["Date"]] <- Sys.time()
    Targeting[["Data"]] <- "Yaari G (2013). Front Immunol."
    Targeting[["Substitution"]] <- S5F_Substitution
    Targeting[["Mutability"]] <- S5F_Mutability
    Targeting[["Targeting"]] <- sweep(S5F_Substitution,2,S5F_Mutability,`*`) 
    save(Targeting,file=paste(Targeting[["Name"]],"_Targeting.RData", sep=""))

}

#### I/O Functions ####

#' Write targeting model to tab-delimited file
#' 
#' \code{writeTargetingModel} writes a five-mer targeting model matrix 
#' of \eqn{mutability * substitution} to a tab-delimited file.
#' 
#' @param    model     \code{\link{TargetingModel}} object with 
#'                     mutation likelihood information.
#' @param    file      name of file to write.
#'                                                
#' @return   NULL
#'           
#' @details
#' \code{writeTargetingModel} takes as input a \code{\link{TargetingModel}} object and
#' writes the \code{targeting} slot, which is the probability of a given mutation 
#' ocurring, defined as \eqn{mutability * substitution}. The targeting model is stored as 
#' a 5x3125 matrix of rates. Rows define the mutated nucleotide at the center of each 5-mer, 
#' one of c(A, C, G, T, N), and columns define the complete 5-mer of the unmutated nucleotide 
#' sequence. It then replaces NAs in this matrix with 0 and writes to a tab-delimited file 
#' with column and row names.
#'    
#' @seealso  Takes as input a \code{\link{TargetingModel}} object.
#' 
#' @examples
#' \dontrun{
#' # Write HS5F targeting model to working directory as hs5f.tab
#' writeTargetingModel(HS5FModel, "hs5f.tab") 
#' }
#' 
#' @export
writeTargetingModel <- function(model, file) {
    to_write <- as.data.frame(model@targeting)
    to_write[is.na(to_write)] <- 0
    write.table(to_write, file, quote=FALSE, sep="\t")
}
