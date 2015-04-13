# Targeting models
# 
# @author     Gur Yaari, Mohamed Uduman, Jason Anthony Vander Heiden
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.1
# @date       2015.04.11

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
#' @seealso  See \code{\link{createMutabilityMatrix}} and \code{\link{createSubstitutionMatrix}} 
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
#' \code{createSubstitutionMatrix} builds a 5-mer nucleotide substitution model by counting 
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
#' @return   A 4x1024 matrix of normalized substitution rate for each 5-mer motif with 
#'           rows names defining the center nucleotide, one of \code{c("A", "C", "G", "T", "N")}, 
#'           and column names defining the 5-mer nucleotide sequence.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
#'            on synonymous mutations from high-throughput immunoglobulin sequencing data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#'
#' @seealso  See \code{\link{createMutabilityMatrix}} for building a mutability model.
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(file)
#' 
#' # Create model using all mutations
#' sub <- createSubstitutionMatrix(db)
#' 
#' # Create model using only silent mutations
#' sub <- createSubstitutionMatrix(db, model="S")
#' 
#' @export
createSubstitutionMatrix <- function(db, model=c("RS", "S"), sequenceColumn="SEQUENCE_IMGT",
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
    nuc_words <- seqinr::words(4, nuc_chars)
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
    nuc_words = seqinr::words(4, nuc_chars)
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
        c( paste(substring(seqinr::words(3,nuc_chars),1,1), "C",substring(seqinr::words(3,nuc_chars),2,3),sep=""),
           paste(substring(seqinr::words(3,nuc_chars),1,1),"G",substring(seqinr::words(3,nuc_chars),2,3),sep=""),
           paste(substring(seqinr::words(3,nuc_chars),1,1),"T",substring(seqinr::words(3,nuc_chars),2,3),sep="")
        )
    
    for(CT in CT_FiveMers){
        CTS[[CT]]<-nuc_words[sapply(nuc_words,function(x)substring(x,1,4))== CT]
        INDEx<-NAMES[2,]%in%CTS[[CT]]
        M[[CT]]<- t(sapply(1:4,function(i)apply(listSubstitution[INDEx,i,],2,sum)))
        rownames(M[[CT]]) <- nuc_chars
    }
    M[["E"]]<-M[[paste(substring(seqinr::words(3,nuc_chars),1,1),"C",substring(seqinr::words(3,nuc_chars),2,3),sep="")[1]]]
    for(CT in c(paste(substring(seqinr::words(3,nuc_chars),1,1),"C",substring(seqinr::words(3,nuc_chars),2,3),sep=""),paste(substring(seqinr::words(3,nuc_chars),1,1),"G",substring(seqinr::words(3,nuc_chars),2,3),sep=""),paste(substring(seqinr::words(3,nuc_chars),1,1),"T",substring(seqinr::words(3,nuc_chars),2,3),sep=""))[-1]){
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
    
    substitutionModel <- sapply(seqinr::words(5, nuc_chars), function(x) { .simplifivemer(M, x) })
    substitutionModel <- apply(substitutionModel, 2, function(x) { x/sum(x, na.rm=TRUE) })
    substitutionModel[!is.finite(substitutionModel)] <- NA
    
    return(substitutionModel)
}


#' Builds a mutability model
#'
#' \code{createMutabilityMatrix} builds a 5-mer nucleotide mutability model by counting 
#' the number of mutations occuring in the center position for all 5-mer motifs.
#'
#' @param    db                 data.frame containing sequence data.
#' @param    substitutionModel  matrix of 5-mer substitution rates built by 
#'                              \code{\link{createSubstitutionMatrix}}.
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
#' @return   A named numeric vector of 1024 normalized mutability rates for each 5-mer 
#'           motif with names 5-mer nucleotide sequence.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution 
#'            based on synonymous mutations from high-throughput immunoglobulin sequencing 
#'            data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#' 
#' @family   targeting model functions
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(file)
#' 
#' # Create model using all mutations
#' sub_model <- createSubstitutionMatrix(db)
#' mut_model <- createMutabilityMatrix(db, sub_model)
#'
#' # Create model using only silent mutations
#' sub_model <- createSubstitutionMatrix(db, model="S")
#' mut_model <- createMutabilityMatrix(db, sub_model, model="S")
#' 
#' @export
createMutabilityMatrix <- function(db, substitutionModel, model=c("RS", "S"),
                                  sequenceColumn="SEQUENCE_IMGT", 
                                  germlineColumn="GERMLINE_IMGT_D_MASK",
                                  vCallColumn="V_CALL",
                                  multipleMutation=c("independent", "ignore")) {
    # model="RS"
    # sequenceColumn="SEQUENCE_IMGT"
    # germlineColumn="GERMLINE_IMGT_D_MASK"
    # vCallColumn="V_CALL"
    # multipleMutation="independent"
    
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
    
    # Count mutations
    nuc_chars <- NUCLEOTIDES[1:4]
    mutations <- listObservedMutations(db, sequenceColumn=sequenceColumn, 
                                       germlineColumn=germlineColumn)
    
    # Do something
    template <- rep(0, 1024)
    names(template) <- seqinr::words(5, nuc_chars)
    COUNT <- list()
    for(index in 1:length(mutations)){
        COUNT[[index]] <- template
        indexMutation <- mutations[[index]]
        if(!sum(is.na(indexMutation))){
            cSeq <-  s2c(db[index, sequenceColumn])
            cGL  <-  s2c(db[index, germlineColumn])
            positions <- as.numeric(names(indexMutation))
            positions <- positions[positions <= VLENGTH]
            positions <- positions[!is.na(positions)]
            for (position in  positions){
                wrd5 <- substr(db[index, germlineColumn], position - 2, position + 2)
                if(!grepl("[^ACGT]", wrd5) & nchar(wrd5)==5){
                    codonNucs = getCodonPos(position)
                    codonGL = cGL[codonNucs]
                    codonSeq = cSeq[codonNucs]
                    muCodonPos = {position - 1} %% 3 + 1
                    seqAtMutation <- codonSeq[muCodonPos]
                    glAtMutation <- codonGL[muCodonPos]
                    if (!any(codonGL %in% c("N", "-", ".")) & !any(codonSeq %in% c("N", "-", "."))) {
                        if (multipleMutation == "independent") { # if independent mutations are included
                            codonSeq <- codonGL
                            codonSeq[muCodonPos] <- seqAtMutation
                        }
                        if (!length(grep("N", wrd5))) {
                            COUNT[[index]][wrd5]<- COUNT[[index]][wrd5] + 1;
                        }
                    }
                }
            }
        }
    }
    
    #cat("Progress: 0%      50%     100%\n")
    #cat("          ")
    #pb <- txtProgressBar(min=1, max=length(mutations), width=20)

    # Do something else
    BG_COUNT <- list()
    for (index in 1:length(mutations)){
        #setTxtProgressBar(pb, index)
        BG_COUNT[[index]] <- template
        sSeq <- gsub("\\.", "", db[index, sequenceColumn])
        sGL <- gsub("\\.", "", db[index, germlineColumn])
        cSeq <-  s2c(sSeq)
        cGL  <-  s2c(sGL)[1:VLENGTH]
        positions <- 3:(length(cGL) - 2)
        for (position in  positions) {
            wrd5 <- substr(sGL, position - 2, position + 2)
            if (!grepl("[^ACGT]", wrd5) & nchar(wrd5) == 5 ) {
                codonNucs = getCodonPos(position)
                codonGL = cGL[codonNucs]
                codonSeq = cSeq[codonNucs]
                muCodonPos = {position - 1} %% 3 + 1
                seqAtMutation <- codonSeq[muCodonPos]
                glAtMutation <- codonGL[muCodonPos]
                if (!any(codonGL %in% c("N", "-")) & !any(codonSeq %in% c("N", "-"))) {
                    if(multipleMutation == "independent") { # if independent mutations are included
                        codonSeq <- codonGL
                        codonSeq[muCodonPos] <- seqAtMutation
                    }
                    codonPermutate <- matrix(rep(codonGL, 3), ncol=3, byrow=TRUE)
                    codonPermutate[, muCodonPos] <- canMutateTo(glAtMutation)[-4]
                    codonPermutate <- apply(codonPermutate, 1, paste,collapse="")
                    codonPermutate <- matrix(c(codonPermutate, rep(c2s(codonGL), 3)), ncol=2, byrow=FALSE)
                    muType <- mutationTypeOptimized(codonPermutate)
                    if (!length(grep("N", wrd5)) & !length(grep("-", wrd5))) {
                        for (m in 1:3) {
                            if (muType[m] == "S") {
                                BG_COUNT[[index]][wrd5]<- BG_COUNT[[index]][wrd5] + 
                                    substitutionModel[substr(codonPermutate[m, 1], muCodonPos, muCodonPos), wrd5];
                            }
                        }
                    }
                }
            }
        }
        BG_COUNT[[index]][BG_COUNT[[index]] == 0] <- NA
    }
    #close(pb)
    #cat("\n")
    
    Mutability<-list()
    for(i in 1:length(mutations)){
        mut_mat <- COUNT[[i]] / BG_COUNT[[i]]
        mut_mat <- mut_mat / sum(mut_mat, na.rm=TRUE)
        mut_mat[!is.finite(mut_mat)] <- NA
        wgt_mat <- length(mutations[[i]])
        Mutability[[i]] <- list(mut_mat, wgt_mat)
        
        #Mutability[[i]] <- list()
        #Mutability[[i]][[1]]<-COUNT[[i]]/BG_COUNT[[i]]
        #Mutability[[i]][[1]]<-Mutability[[i]][[1]]/sum(Mutability[[i]][[1]],na.rm=TRUE)
        #Mutability[[i]][[2]]<-length( mutations[[i]])
        #Mutability[[i]][[3]]<-sum(COUNT[[i]],na.rm=TRUE)
        #Mutability[[i]][[1]][!is.finite(Mutability[[i]][[1]])] <- NA
    }
    
    
    # Aggregate mutability
    MutabilityMatrix <- sapply(Mutability, function(x) x[[1]])
    MutabilityWeights <- sapply(Mutability, function(x) x[[2]])
    Mutability_Mean <- apply(MutabilityMatrix, 1, weighted.mean, w=MutabilityWeights, na.rm=TRUE)
    Mutability_Mean[!is.finite(Mutability_Mean)] <- NA
    
    # Normalize
    Mutability_Mean <- Mutability_Mean / sum(Mutability_Mean, na.rm=TRUE)
    
    # TODO:  what is this?
    # MutabilityMatrixNorm = Normalized by fivemers that are observed.
    #Mutability_Mean <- sapply(1:1024, function(i) weighted.mean(MutabilityMatrix[i, ], MutabilityWeights[i, ], na.rm=TRUE))
    #MutabilityMatrixNorm <- apply(MutabilityMatrix, 2, function(x) x/sum(x, na.rm=TRUE) * sum(!is.na(x)))
    #Mutability_SD <- sapply(1:1024, function(i) SDMTools::wt.sd(MutabilityMatrixNorm[i, ], MutabilityWeights[i]))    
    
    return(Mutability_Mean)
}


#' Extends a substitution model to include Ns.
#' 
#' \code{extendSubstitutionMatrix} extends a 5-mer nucleotide substitution model 
#' with 5-mers that include Ns by averaging over all corresponding 5-mers without Ns.
#'
#' @param    substitutionModel  matrix of 5-mers substitution counts built by 
#'                              \code{\link{createSubstitutionMatrix}}.
#' 
#' @return   A 5x3125 matrix of normalized substitution rate for each 5-mer motif with 
#'           rows names defining the center nucleotide, one of \code{c("A", "C", "G", "T", "N")}, 
#'           and column names defining the 5-mer nucleotide sequence.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution 
#'            based on synonymous mutations from high-throughput immunoglobulin sequencing 
#'            data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#' 
#' @family   targeting model functions
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(file)
#' 
#' # Create model using all mutations
#' sub_model <- createSubstitutionMatrix(db)
#' ext_model <- extendSubstitutionMatrix(sub_model)
#' 
#' @export
extendSubstitutionMatrix <- function(substitutionModel) {
    # Define old and new column/row names
    input_names <- colnames(substitutionModel)
    nuc_chars <- NUCLEOTIDES[1:5]
    nuc_5mers <- seqinr::words(5, alphabet=nuc_chars)
    
    # Define empty extended matrix with Ns
    extend_mat <- matrix(NA, nrow=length(nuc_chars), ncol=length(nuc_5mers), 
                         dimnames=list(nuc_chars, nuc_5mers))
    
    # Extend matrix with Ns
    for (mer in nuc_5mers) {
        if (mer %in% input_names) {
            extend_mat[, mer] <- c(substitutionModel[, mer], "N"=NA)
        } else {
            mer_char <- s2c(mer)
            n_index <- grep("N", mer_char)
            if (any(n_index == 3)) {
                extend_mat[, mer] <- NA
            } else {
                mer_char[n_index] <- "."
                mer_str <- c2s(mer_char)
                mer_index <- grep(mer_str, input_names)
                extend_mat[, mer] <- c(apply(substitutionModel[, mer_index], 1, mean, na.rm=TRUE), "N"=NA)
            }
        }
    }
    
    # Normalize
    extend_mat <- apply(extend_mat, 2, function(x) { x/sum(x, na.rm=TRUE) })
    extend_mat[!is.finite(extend_mat)] <- NA
    
    #extend_mat <- extend_mat / sum(extend_mat, na.rm=TRUE)
    
    return (extend_mat)
}


#' Extends a mutability model to include Ns.
#' 
#' \code{extendMutabilityMatrix} extends a 5-mer nucleotide mutability model 
#' with 5-mers that include Ns by averaging over all corresponding 5-mers without Ns.
#'
#' @param    mutabilityModel  vector of 5-mer mutability rates built by 
#'                            \code{\link{createMutabilityMatrix}}.
#' 
#' @return   A 3125 vector of normalized mutability rates for each 5-mer motif with 
#'           names defining the 5-mer nucleotide sequence.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution 
#'            based on synonymous mutations from high-throughput immunoglobulin sequencing 
#'            data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#' 
#' @family   targeting model functions
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(file)
#' 
#' # Create model using all mutations
#' sub_model <- createSubstitutionMatrix(db)
#' mut_model <- createMutabilityMatrix(db, sub_model)
#' ext_model <- extendMutabilityMatrix(mut_model)
#' 
#' @export
extendMutabilityMatrix <- function(mutabilityModel) {
    # TODO: fix order so Ns are at the end? (c(input_names, words not in input_names))
    
    # Define old and new column/row names
    input_names <- names(mutabilityModel)
    nuc_chars <- NUCLEOTIDES[1:5]
    nuc_5mers <- seqinr::words(5, alphabet=nuc_chars)
    
    # Define empty extended matrix with Ns
    extend_mat <- array(NA, dim=length(nuc_5mers), dimnames=list(nuc_5mers))
    
    # Extend matrix with Ns
    for(mer in nuc_5mers) {
        #cat(mer,"\n")
        if (mer %in% input_names) {
            extend_mat[mer] <- mutabilityModel[mer]
        } else {
            mer_char <- s2c(mer)
            n_index <- grep("N", mer_char)
            if (any(n_index == 3)) {
                extend_mat[mer] <- NA
            } else {
                mer_char[n_index] <- "."
                mer_str <- c2s(mer_char)
                mer_index <- grep(mer_str, input_names)
                extend_mat[mer] <- mean(mutabilityModel[mer_index], na.rm=TRUE)
            }
        }
    }
    
    # Normalize    
    extend_mat <- extend_mat / sum(extend_mat, na.rm=TRUE)
    extend_mat[!is.finite(extend_mat)] <- NA
    
    return(extend_mat)
}
 

#' Calculates a targeting rate matrix
#' 
#' \code{createTargetingMatrix} calculates the targeting model matrix as the
#' combined probability of mutability and substitution.
#'
#' @param    substitutionModel  matrix of 5-mers substitution rates built by 
#'                              \code{\link{createSubstitutionMatrix}} or 
#'                              \code{\link{extendSubstitutionMatrix}}.
#' @param    mutabilityModel    vector of 5-mers mutability rates built by 
#'                              \code{\link{createMutabilityMatrix}} or 
#'                              \code{\link{extendMutabilityMatrix}}.
#' 
#' @return   A matrix with the same dimensions as the input \code{substitutionModel} 
#'           containing normalized targeting probabilities for each 5-mer motif with 
#'           rows names defining the center nucleotide and column names defining the 
#'           5-mer nucleotide sequence.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution 
#'            based on synonymous mutations from high-throughput immunoglobulin sequencing 
#'            data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#' 
#' @family   targeting model functions
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(file)
#' 
#' # Create 4x1024 model using all mutations
#' sub_model <- createSubstitutionMatrix(db)
#' mut_model <- createMutabilityMatrix(db, sub_model)
#' tar_model <- createTargetingMatrix(sub_model, mut_model)
#' 
#' # Create 5x3125 model including Ns
#' sub_model <- extendSubstitutionMatrix(sub_model)
#' mut_model <- extendMutabilityMatrix(mut_model)
#' tar_model <- createTargetingMatrix(sub_model, mut_model)
#' 
#' @export
createTargetingMatrix <- function(substitutionModel, mutabilityModel) {
    # Calculate targeting
    tar_mat <- sweep(substitutionModel, 2, mutabilityModel, `*`)
    
    # Normalize    
    tar_mat <- tar_mat / sum(tar_mat, na.rm=TRUE)
    tar_mat[!is.finite(tar_mat)] <- NA
    
    return(tar_mat)
}


#' Creates a TargetingModel
#' 
#' \code{createTargetingModel} creates a \code{TargetingModel}.
#'
#' @param    db                 data.frame containing sequence data.
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
#' @param    modelName          name of the model.
#' @param    modelDescription   description of the model and its source data.
#' @param    modelSpecies       genus and species of the source sequencing data.
#' @param    modelDate          date the model was built. If \code{NULL} the current date
#'                              will be used.
#' @param    modelCitation      publication source.
#' 
#' @return   A \code{\link{TargetingModel}} object.
#' 
#' @references
#' \enumerate{
#'   \item  Yaari G, et al. Models of somatic hypermutation targeting and substitution 
#'            based on synonymous mutations from high-throughput immunoglobulin sequencing 
#'            data. 
#'            Front Immunol. 2013 4(November):358.
#'  }
#' 
#' @seealso  See \code{\link{TargetingModel}} for the return object.
#' @family   targeting model functions
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(file)
#' 
#' # Create model using all mutations
#' model <- createTargetingModel(db)
#'
#' # Create model using only silent mutations
#' model <- createTargetingModel(db, model="S")
#' 
#' @export
createTargetingModel <- function(db, model=c("RS", "S"), sequenceColumn="SEQUENCE_IMGT",
                                 germlineColumn="GERMLINE_IMGT_D_MASK",
                                 vCallColumn="V_CALL",
                                 multipleMutation=c("independent", "ignore"),
                                 modelName="", modelDescription="", modelSpecies="", 
                                 modelCitation="", modelDate=NULL) {
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
    
    # Set date
    if (is.null(modelDate)) { modelDate <- format(Sys.time(), "%Y-%m-%d") }

    # Create models
    sub_model <- createSubstitutionMatrix(db)
    mut_model <- createMutabilityMatrix(db, sub_model)

    # Extend 5-mers with Ns
    sub_model <- extendSubstitutionMatrix(sub_model)
    mut_model <- extendMutabilityMatrix(mut_model)
    
    # TODO: this is wrong somehow
    tar_model <- createTargetingMatrix(sub_model, mut_model) 
    
    # Define TargetingModel object
    model_obj <- new("TargetingModel",
                     name=modelName,
                     description=modelDescription,
                     species=modelSpecies,
                     date=modelDate,
                     citation=modelCitation,
                     substitution=sub_model,
                     mutability=mut_model,
                     targeting=tar_model)

    return(model_obj)
}


#' Calculates a 5-mer distance matrix from a TargetingModel object
#' 
#' \code{getTargetingDistance} renormalizes the targeting probabilities
#' in a TargetingModel model and returns a distance matrix.
#' 
#' @param    model     \code{\link{TargetingModel}} object with 
#'                     mutation likelihood information.
#'                                                
#' @return   A matrix of distances for each 5-mer motif with rows names defining 
#'           the center nucleotide and column names defining the 5-mer nucleotide 
#'           sequence.
#'           
#' @details
#' The targeting distance is defined as the targeting probability of a given base change 
#' from each 5-mer (columns) to each center nucleotide (rows), re-centered such that the
#' mean distance for any base change is 0.25.
#'    
#' @seealso  Takes as input a \code{\link{TargetingModel}} object.
#' @family   targeting model functions
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(file)
#' 
#' # Create model and get targeting distance
#' model <- createTargetingModel(db)
#' dist <- getTargetingDistance(model)
#' 
#' @export
getTargetingDistance <- function(model) {
    if (is(model, "TargetingModel")) {
        model <- model@targeting
    } else if (!is(model, "matrix")) {
        stop("Input must be either a targeting matrix or TargetingModel object.")
    }
    
    # TODO:  need to fix normalization in extend functions
    # TODO:  need to make A->A, C->C, G->G, T->T NA.
    
    # Get mean of 5-mers without Ns
    nuc_chars <- NUCLEOTIDES[1:4]
    nuc_5mers <- words(5, nuc_chars)
    index_5mers <- colnames(model) %in% nuc_5mers
    dist_5mers <- 1 - model[nuc_chars, index_5mers]
    mean_5mers <- mean(dist_5mers, na.rm=TRUE)
    
    # Convert to distance
    dist <- 1 - model
    dist <- 1 + dist - mean_5mers    
    dist[!is.finite(dist)] <- NA
    
    return(dist)
}


#' Rescales mutability probabilities from a TargetingModel
#' 
#' \code{getRescaledMutability} renormalizes the mutability probabilities
#' in a TargetingModel model and returns a rescaled matrix of mutability scores.
#' 
#' @param    model     \code{\link{TargetingModel}} object with 
#'                     mutation likelihood information.
#' @param    mean      the mean value for the rescaled mutability scores.
#'                                                
#' @return   A named vector of mutability scores for each 5-mer motif with mean
#'           equal to \code{mean}.
#'           
#' @seealso  Takes as input a \code{\link{TargetingModel}} object.
#' @family   targeting model functions
#' 
#' @examples
#' # Load example data
#' library(alakazam)
#' file <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(file)
#' 
#' # Create model and get targeting distance
#' model <- createTargetingModel(db)
#' mut <- getRescaledMutability(model)
#' 
#' @export
getRescaledMutability <- function(model, mean=1.0) {
    if (is(model, "TargetingModel")) {
        model <- model@mutability
    } else if (!is(model, "vector")) {
        stop("Input must be either a mutability vector or TargetingModel object.")
    }
    
    # Renormalize
    rescaled <- model / sum(model, na.rm=T) * sum(!is.na(model)) * mean
    rescaled[!is.finite(rescaled)] <- NA
    
    return(rescaled)
}


#### I/O Functions ####

#' Write targeting model distances to a file
#' 
#' \code{writeTargetingDistance} writes a 5-mer targeting distance matrix 
#' to a tab-delimited file.
#' 
#' @param    model     \code{\link{TargetingModel}} object with 
#'                     mutation likelihood information.
#' @param    file      name of file to write.
#'                                                
#' @return   NULL
#'           
#' @details
#' Takes as input a \code{\link{TargetingModel}} object and writes the re-centered 
#' \code{targeting} slot, which is the probability of a given mutation ocurring, 
#' defined as \eqn{mutability * substitution}. The targeting distance is stored as 
#' a 5x3125 matrix of rates. Rows define the mutated nucleotide at the center of each 5-mer, 
#' one of \code{c("A", "C", "G", "T", "N")}, and columns define the complete 5-mer of the 
#' unmutated nucleotide sequence. \code{NA} values in the distance matrix are replaced with 
#' distance 0. This writes to a tab-delimited file with column and row names.
#'    
#' @seealso  Takes as input a \code{\link{TargetingModel}} object.
#' @family   targeting model functions
#' 
#' @examples
#' \dontrun{
#' # Write HS5F targeting model to working directory as hs5f.tab
#' writeTargetingModel(HS5FModel, "hs5f.tab") 
#' }
#' 
#' @export
writeTargetingDistance <- function(model, file) {
    to_write <- as.data.frame(getTargetingDistance(model))
    to_write[is.na(to_write)] <- 0
    write.table(to_write, file, quote=FALSE, sep="\t")
}


#### Plotting functions ####

#' Plot mutability probabilities
#' 
#' \code{plotMutability} creates a hedgehog plot of the mutability rates.
#' 
#' @param    model        \code{\link{TargetingModel}} object or mutability matrix
#'                        with mutability probabilities.
#' @param    nucleotides  vector of center nucleotide characters to plot mutability for.
#'                                                
#' @return   NULL
#'    
#' @seealso  Takes as input a \code{\link{TargetingModel}} object.
#' @family   targeting model functions
#' 
#' @examples
#' # Plot all nucleotides
#' plotMutability(HS5FModel)
#' 
#' # Plot one nucleotide
#' plotMutability(HS5FModel, "C")
#' 
#' @export
plotMutability <- function(model, nucleotides=c("A", "C", "G", "T")) {
    # Check input
    if (is(model, "TargetingModel")) {
        model <- model@mutability
    } else if (!is(model, "vector")) {
        stop("Input must be either a mutability vector or TargetingModel object.")
    }
    nucleotides <- toupper(nucleotides)

    # Set base plot settings
    base_theme <- theme_bw() +
        theme(panel.border=element_blank(), 
              axis.text=element_blank(), 
              axis.ticks=element_blank(),
              legend.position="bottom")
        
    # Set guide colors
    motif_colors <- setNames(c("#33a02c", "#e31a1c", "#6a3d9a", "#999999"),
                             c("WA/TW", "WRC/GYW", "SYC/GRS", "Neutral"))
    dna_colors <- setNames(c("#b2df8a", "#fdbf6f", "#fb9a99", "#a6cee3", "#aaaaaa"),
                           c("A", "C", "G", "T", "N"))    
    #motif_colors <- setNames(c("#4daf4a", "#e41a1c", "#377eb8", "#999999"),
    #                         c("WA/TW", "WRC/GYW", "SYC/GRS", "Neutral"))
    #dna_colors <- setNames(c("#ccebc5", "#fed9a6", "#fbb4ae", "#b3cde3"),
    #                       c("A", "C", "G", "T"))
    
    # Build data.frame of mutability scores
    mut_scores <- model[!grepl("N", names(model))]
    mut_scores <- getRescaledMutability(mut_scores)
    mut_scores[!is.finite(mut_scores)] <- 0
    mut_words <- names(mut_scores)
    mut_positions <- as.data.frame(t(sapply(mut_words, s2c)))
    colnames(mut_positions) <- paste0("pos", 1:ncol(mut_positions))
    mut_df <- data.frame(word=mut_words, 
                         score=mut_scores, 
                         mut_positions)
    
    # Add hot and cold-spot motif information
    mut_df$motif <- "Neutral"
    mut_df$motif[grepl("(.[AT]A..)|(..T[AT].)", mut_df$word, perl=TRUE)] <- "WA/TW"
    mut_df$motif[grepl("([AT][GA]C..)|(..G[CT][AT])", mut_df$word, perl=TRUE)] <- "WRC/GYW"
    mut_df$motif[grepl("([CG][CT]C..)|(..G[GA][CG])", mut_df$word, perl=TRUE)] <- "SYC/GRS"
    mut_df$motif <- factor(mut_df$motif, levels=c("WA/TW", "WRC/GYW", "SYC/GRS", "Neutral"))

    # Build plots for each center nucleotide
    plot_list <- list()
    for (center_nuc in nucleotides) {
        # Subset data to center nucleotide
        sub_df <- subset(mut_df, pos3 == center_nuc)
        sub_df$x <- 1:nrow(sub_df)
        
        # Melt 5-mer position data
        sub_melt <- reshape2::melt(sub_df, id.vars=c("x"), 
                                 measure.vars=colnames(mut_positions),
                                 variable.name="pos",
                                 value.name="char")
        sub_melt$pos <- factor(sub_melt$pos, levels=colnames(mut_positions))
        sub_melt$pos <- as.numeric(sub_melt$pos)
        
        # Define text position data
        sub_text <- list()
        for (i in 1:5) {
            tmp_rle <- rle(sub_melt$char[sub_melt$pos == i])
            tmp_sum <- cumsum(tmp_rle$lengths)
            tmp_pos <- tmp_sum - diff(c(0, tmp_sum))/2
            if (length(tmp_pos) == 1) { tmp_pos = 1/2 }
            
            tmp_df <- data.frame(text_x=tmp_pos, 
                                  text_y=i,
                                  text_label=tmp_rle$values)
            
            sub_text[[i]] <- tmp_df
        }
        
        y_limits <- c(-1, max(sub_df$score) * 2 + 5.6)
        # Plot mutability for center nucleotide
        p1 <- ggplot(sub_df, aes(x=x)) + 
            base_theme + 
            ggtitle(paste0("NN", center_nuc, "NN")) +
            xlab("") +
            ylab("") + 
            coord_polar(theta="x") +
            scale_y_continuous(limits=y_limits, expand=c(0, 0)) +
            scale_color_manual(name="Motif", values=motif_colors) +
            scale_fill_manual(name="Nucleotide", values=dna_colors, guide=FALSE) +
            geom_segment(data=sub_df, mapping=aes(x=x, xend=x, y=5.6, yend=5.6 + score * 2, color=motif), 
                         size=2) +
            #geom_bar(data=sub_df, mapping=aes(x=x, y=5 + score * 2, fill=motif), stat="identity", 
            #         position="identity", size=2) +
            #geom_rect(xmin=0, xmax = Inf, ymin=0, ymax=5, fill="white") +
            geom_tile(data=sub_melt, mapping=aes(x=x, y=pos, fill=char), size=0) +
            geom_text(data=sub_text[[1]], mapping=aes(x=text_x, y=text_y, label=text_label), 
                      color="black", hjust=0.5, vjust=0.5, size=3) +
            geom_text(data=sub_text[[2]], mapping=aes(x=text_x, y=text_y, label=text_label), 
                      color="black", hjust=0.5, vjust=0.5, size=3) +
            geom_text(data=sub_text[[3]], mapping=aes(x=text_x, y=text_y, label=text_label), 
                      color="black", hjust=0.5, vjust=0.5, size=4) +
            geom_text(data=sub_text[[4]], mapping=aes(x=text_x, y=text_y, label=text_label), 
                      color="black", hjust=0.5, vjust=0.5, size=2.5)
            #geom_rect(xmin=0, xmax = Inf, ymin=0, ymax=0.5, fill="white")
        
        # Add plot to list
        plot_list[[center_nuc]] <- p1
    }
    
    # Plot
    do.call(multiggplot, args=c(plot_list, ncol=4))
}


# TODO:  put this in alakazam and export it
# Plot multiple ggplot objects
# 
# @param   ...    ggplot objects to plot
# @param   ncol   number of columns in the plot 
# @return  NULL
# 
# @references  http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)
multiggplot <- function(..., ncol=1) {
    p <- list(...)
    n <- length(p)
    layout <- matrix(seq(1, ncol*ceiling(n/ncol)), ncol=ncol, nrow=ceiling(n/ncol))
    
    # Plot
    if (n == 1) {
        plot(p[[1]])
    } else {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout=grid::grid.layout(nrow(layout), ncol(layout))))
        for (i in 1:n) {
            idx <- as.data.frame(which(layout == i, arr.ind=T))
            plot(p[[i]], vp=grid::viewport(layout.pos.row = idx$row, layout.pos.col=idx$col))
        }
    }
}


#### Testing functions ####

# Function to make dummy data for testing targetting functions
# 
# @param   nseq  number of sequences
# @param   nmut  number of mutations per sequence
# @param   nmer  number of 5-mers per sequence (sequence length = 5 * nmer)
#
# @return  a data.frame with columns SEQUENCE_ID, SEQUENCE_IMGT, GERMLINE_IMGT_D_MASK, V_CALL.
makeTargetingTestDb <- function(nseq=10, nmut=40, nmers=50) {
    nuc_chars <- c("A", "C", "G", "T")
    
    .mut <- function(x, n) {
       i <- sample(1:nchar(x), n) 
       y <- seqinr::s2c(x)
       y[i] <- sapply(y[i], function(z) sample(nuc_chars[nuc_chars != z], 1))
       return(seqinr::c2s(y))
    }
    
    seq <- apply(replicate(nseq, sample(seqinr::words(5, nuc_chars), nmers)), 2, paste, collapse="")
    germ <- sapply(seq, .mut, n=nmut)
    db <- data.frame(SEQUENCE_ID=paste0("SEQ", 1:nseq),
                     SEQUENCE_IMGT=seq,
                     GERMLINE_IMGT_D_MASK=germ,
                     V_CALL="Homsap IGHV3-66*02 F", stringsAsFactors=FALSE)
    rownames(db) <- NULL
    
    return(db)
}