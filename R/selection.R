# BASELINe
# 
# @author     Mohamed Uduman
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.1
# @date       2015.03.06

#' @include shm.R
NULL

#### Functions ####



# Translate amino acid to trait change
translateAminoAcidToTraitChange<-function(AminoAcid){
    return(TRAITS_AMINO_ACIDS[AminoAcid])
}

# Initialize Amino Acid Trait Changes
initializeTraitChange <- function(traitChangeModel=1,species=1,traitChangeFileName=NULL){
    if(!is.null(traitChangeFileName)){
        tryCatch(
            traitChange <- read.delim(traitChangeFileName,sep="\t",header=T)
            , error = function(ex){
                cat("Error|Error reading trait changes. Please check file name/path and format.\n")
                q()
            }
        )
    }else{
        traitChange <- TRAITS_AMINO_ACIDS_CHOTHIA98
    }
    TRAITS_AMINO_ACIDS <<- traitChange
}

# Read in formatted nucleotide substitution matrix
initializeSubstitutionMatrix <- function(substitutionModel,species,subsMatFileName=NULL){
    if(!is.null(subsMatFileName)){
        tryCatch(
            subsMat <- read.delim(subsMatFileName,sep="\t",header=T)
            , error = function(ex){
                cat("Error|Error reading substitution matrix. Please check file name/path and format.\n")
                q()
            }
        )
        if(sum(apply(subsMat,1,sum)==1)!=4) subsMat = t(apply(subsMat,1,function(x)x/sum(x)))
    }else{
        if(substitutionModel==1)subsMat <- substitution_Literature_Mouse
        if(substitutionModel==2)subsMat <- substitution_Flu_Human
        if(substitutionModel==3)subsMat <- substitution_Flu25_Human
    }
    
    if(substitutionModel==0){
        subsMat <- matrix(1,4,4)
        subsMat[,] = 1/3
        subsMat[1,1] = 0
        subsMat[2,2] = 0
        subsMat[3,3] = 0
        subsMat[4,4] = 0
    }
    
    
    NUCLEOTIDESN = c(NUCLEOTIDES,"N", "-")
    if(substitutionModel==5){
        subsMat <- S5F_Targeting[["Substitution"]]
        return(subsMat)
    }else{
        subsMat <- rbind(subsMat,rep(NA,4),rep(NA,4))
        return( matrix(data.matrix(subsMat),6,4,dimnames=list(NUCLEOTIDESN,NUCLEOTIDES) ) )
    }
}


# Read in formatted Mutability file
initializeMutabilityMatrix <- function(mutabilityModel=1, species=1,mutabilityMatFileName=NULL){
    if(!is.null(mutabilityMatFileName)){
        tryCatch(
            mutabilityMat <- read.delim(mutabilityMatFileName,sep="\t",header=T)
            , error = function(ex){
                cat("Error|Error reading mutability matrix. Please check file name/path and format.\n")
                q()
            }
        )
    }else{
        mutabilityMat <- triMutability_Literature_Human
        if(species==2) mutabilityMat <- triMutability_Literature_Mouse
    }
    
    if(mutabilityModel==0){ mutabilityMat <- matrix(1,64,3)}
    
    if(mutabilityModel==5){
        mutabilityMat <- HS5FModel@mutability
        return(mutabilityMat)
    }else{
        return( matrix( data.matrix(mutabilityMat), 64, 3, dimnames=list(triMutability_Names,1:3)) )
    }
}


# Read FASTA file formats
# Modified from read.fasta from the seqinR package
baseline.read.fasta <-
    function (file = system.file("sequences/sample.fasta", package = "seqinr"),
              seqtype = c("DNA", "AA"), as.string = FALSE, forceDNAtolower = TRUE,
              set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE,
              strip.desc = FALSE,  sizeof.longlong = .Machine$sizeof.longlong,
              endian = .Platform$endian, apply.mask = TRUE)
    {
        seqtype <- match.arg(seqtype)
        
        lines <- readLines(file)
        
        if (legacy.mode) {
            comments <- grep("^;", lines)
            if (length(comments) > 0)
                lines <- lines[-comments]
        }
        
        
        ind_groups<-which(substr(lines, 1L, 3L) == ">>>")
        lines_mod<-lines
        
        if(!length(ind_groups)){
            lines_mod<-c(">>>All sequences combined",lines)
        }
        
        ind_groups<-which(substr(lines_mod, 1L, 3L) == ">>>")
        
        lines <- array("BLA",dim=(length(ind_groups)+length(lines_mod)))
        id<-sapply(1:length(ind_groups),function(i)ind_groups[i]+i-1)+1
        lines[id] <- "THIS IS A FAKE SEQUENCE"
        lines[-id] <- lines_mod
        rm(lines_mod)
        
        ind <- which(substr(lines, 1L, 1L) == ">")
        nseq <- length(ind)
        if (nseq == 0) {
            stop("no line starting with a > character found")
        }
        start <- ind + 1
        end <- ind - 1
        
        while( any(which(ind%in%end)) ){
            ind=ind[-which(ind%in%end)]
            nseq <- length(ind)
            if (nseq == 0) {
                stop("no line starting with a > character found")
            }
            start <- ind + 1
            end <- ind - 1
        }
        
        end <- c(end[-1], length(lines))
        sequences <- lapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]], collapse = ""))
        if (seqonly)
            return(sequences)
        nomseq <- lapply(seq_len(nseq), function(i) {
            
            #firstword <- strsplit(lines[ind[i]], " ")[[1]][1]
            substr(lines[ind[i]], 2, nchar(lines[ind[i]]))
            
        })
        if (seqtype == "DNA") {
            if (forceDNAtolower) {
                sequences <- as.list(tolower(chartr(".","-",sequences)))
            }else{
                sequences <- as.list(toupper(chartr(".","-",sequences)))
            }
        }
        if (as.string == FALSE)
            sequences <- lapply(sequences, s2c)
        if (set.attributes) {
            for (i in seq_len(nseq)) {
                Annot <- lines[ind[i]]
                if (strip.desc)
                    Annot <- substr(Annot, 2L, nchar(Annot))
                attributes(sequences[[i]]) <- list(name = nomseq[[i]],
                                                   Annot = Annot, class = switch(seqtype, AA = "SeqFastaAA",
                                                                                 DNA = "SeqFastadna"))
            }
        }
        names(sequences) <- nomseq
        return(sequences)
    }


# Replaces non FASTA characters in input files with N
replaceNonFASTAChars <-function(inSeq="ACGTN-AApA"){
    gsub('[^ACGTNacgt[:punct:]-[:punct:].]','N',inSeq,perl=TRUE)
}

# Find the germlines in the FASTA list
germlinesInFile <- function(seqIDs){
    firstChar = sapply(seqIDs,function(x){substr(x,1,1)})
    secondChar = sapply(seqIDs,function(x){substr(x,2,2)})
    return(firstChar==">" & secondChar!=">")
}

# Find the groups in the FASTA list
groupsInFile <- function(seqIDs){
    sapply(seqIDs,function(x){substr(x,1,2)})==">>"
}

# In the process of finding germlines/groups, expand from the start to end of the group
expandTillNext <- function(vecPosToID){
    IDs = names(vecPosToID)
    posOfInterests =  which(vecPosToID)
    
    expandedID = rep(NA,length(IDs))
    expandedIDNames = gsub(">","",IDs[posOfInterests])
    startIndexes = c(1,posOfInterests[-1])
    stopIndexes = c(posOfInterests[-1]-1,length(IDs))
    expandedID  = unlist(sapply(1:length(startIndexes),function(i){
        rep(i,stopIndexes[i]-startIndexes[i]+1)
    }))
    names(expandedID) = unlist(sapply(1:length(startIndexes),function(i){
        rep(expandedIDNames[i],stopIndexes[i]-startIndexes[i]+1)
    }))
    return(expandedID)
}

# Process FASTA (list) to return a matrix[input, germline)
processInputAdvanced <- function(inputFASTA){
    
    seqIDs = names(inputFASTA)
    numbSeqs = length(seqIDs)
    posGermlines1 = germlinesInFile(seqIDs)
    numbGermlines = sum(posGermlines1)
    posGroups1 = groupsInFile(seqIDs)
    numbGroups = sum(posGroups1)
    consDef = NA
    
    if(numbGermlines==0){
        posGermlines = 2
        numbGermlines = 1
    }
    
    glPositionsSum = cumsum(posGermlines1)
    glPositions = table(glPositionsSum)
    #Find the position of the conservation row
    consDefPos = as.numeric(names(glPositions[names(glPositions)!=0 & glPositions==1]))+1
    if( length(consDefPos)> 0 ){
        consDefID =  match(consDefPos, glPositionsSum)
        #The coservation rows need to be pulled out and stores seperately
        consDef =  inputFASTA[consDefID]
        inputFASTA =  inputFASTA[-consDefID]
        
        seqIDs = names(inputFASTA)
        numbSeqs = length(seqIDs)
        posGermlines1 = germlinesInFile(seqIDs)
        numbGermlines = sum(posGermlines1)
        posGroups1 = groupsInFile(seqIDs)
        numbGroups = sum(posGroups1)
        if(numbGermlines==0){
            posGermlines = 2
            numbGermlines = 1
        }
    }
    
    posGroups <- expandTillNext(posGroups1)
    posGermlines <- expandTillNext(posGermlines1)
    posGermlines[posGroups1] = 0
    names(posGermlines)[posGroups1] = names(posGroups)[posGroups1]
    posInput = rep(TRUE,numbSeqs)
    posInput[posGroups1 | posGermlines1] = FALSE
    
    matInput = matrix(NA, nrow=sum(posInput), ncol=2)
    rownames(matInput) = seqIDs[posInput]
    colnames(matInput) = c("Input","Germline")
    
    vecInputFASTA = unlist(inputFASTA)
    matInput[,1] = vecInputFASTA[posInput]
    matInput[,2] = vecInputFASTA[ which( names(inputFASTA)%in%paste(">",names(posGermlines)[posInput],sep="") )[ posGermlines[posInput]] ]
    
    germlines = posGermlines[posInput]
    groups = posGroups[posInput]
    
    return( list("matInput"=matInput, "germlines"=germlines, "groups"=groups, "conservationDefinition"=consDef ))
}


# Replace leading and trailing dashes in the sequence
replaceLeadingTrailingDashes <- function(x,VLENGTH){
    iiGap = unlist(gregexpr("-",x[1]))
    ggGap = unlist(gregexpr("-",x[2]))
    #posToChange = intersect(iiGap,ggGap)
    
    
    seqIn = replaceLeadingTrailingDashesHelper(x[1])
    seqGL = replaceLeadingTrailingDashesHelper(x[2])
    #seqTemplate = rep('N',VLENGTH)
    seqIn <- c(seqIn,array("N",max(VLENGTH-length(seqIn),0)))
    seqGL <- c(seqGL,array("N",max(VLENGTH-length(seqGL),0)))
    seqIn <- seqIn[1:VLENGTH]
    seqGL <- seqGL[1:VLENGTH]
    #    if(posToChange!=-1){
    #      seqIn[posToChange] = "-"
    #      seqGL[posToChange] = "-"
    #    }
    
    
    N_Index <- grep("-",seqIn)
    if(any(N_Index)){
        NposCodons <- sapply(N_Index,getCodonNumb)
        names(NposCodons) <- N_Index
        tblGapPos <- table(NposCodons)
        seqIn[as.numeric(names(NposCodons)[ tblGapPos[tblGapPos<3] ])] = "N"
    }
    N_Index <- grep("-",seqGL)
    if(any(N_Index)){
        NposCodons <- sapply(N_Index,getCodonNumb)
        names(NposCodons) <- N_Index
        tblGapPos <- table(NposCodons)
        seqGL[as.numeric(names(NposCodons)[ tblGapPos[tblGapPos<3] ])] = "N"
    }
    
    seqIn = c2s(seqIn)
    seqGL = c2s(seqGL)
    
    lenGL = nchar(seqGL)
    if(lenGL<VLENGTH){
        seqGL = paste(seqGL,c2s(rep("N",VLENGTH-lenGL)),sep="")
    }
    
    lenInput = nchar(seqIn)
    if(lenInput<VLENGTH){
        seqIn = paste(seqIn,c2s(rep("N",VLENGTH-lenInput)),sep="")
    }
    return( c(seqIn,seqGL) )
}

replaceLeadingTrailingDashesHelper <- function(x){
    grepResults = gregexpr("-*",x)
    grepResultsPos = unlist(grepResults)
    grepResultsLen =  attr(grepResults[[1]],"match.length")
    x = s2c(x)
    if(x[1]=="-"){
        x[1:grepResultsLen[1]] = "N"
    }
    if(x[length(x)]=="-"){
        x[(length(x)-grepResultsLen[length(grepResultsLen)]+1):length(x)] = "N"
    }
    return(x)
}




# Check sequences for indels
checkForInDels <- function(matInputP){
    insPos <- checkInsertion(matInputP)
    delPos <- checkDeletions(matInputP)
    return(list("Insertions"=insPos, "Deletions"=delPos))
}

# Check sequences for insertions
checkInsertion <- function(matInputP){
    insertionCheck = apply( matInputP,1, function(x){
        inputGaps <- as.vector( gregexpr("-",x[1])[[1]] )
        glGaps <- as.vector( gregexpr("-",x[2])[[1]] )
        return( is.finite( match(FALSE, glGaps%in%inputGaps ) ) )
    })
    return(as.vector(insertionCheck))
}
# Fix inserstions
fixInsertions <- function(matInputP){
    insPos <- checkInsertion(matInputP)
    sapply((1:nrow(matInputP))[insPos],function(rowIndex){
        x <- matInputP[rowIndex,]
        inputGaps <- gregexpr("-",x[1])[[1]]
        glGaps <- gregexpr("-",x[2])[[1]]
        posInsertions <- glGaps[!(glGaps%in%inputGaps)]
        inputInsertionToN <- s2c(x[2])
        inputInsertionToN[posInsertions]!="-"
        inputInsertionToN[posInsertions] <- "N"
        inputInsertionToN <- c2s(inputInsertionToN)
        matInput[rowIndex,2] <<- inputInsertionToN
    })
    return(insPos)
}

# Check sequences for deletions
checkDeletions <-function(matInputP){
    deletionCheck = apply( matInputP,1, function(x){
        inputGaps <- as.vector( gregexpr("-",x[1])[[1]] )
        glGaps <- as.vector( gregexpr("-",x[2])[[1]] )
        return( is.finite( match(FALSE, inputGaps%in%glGaps ) ) )
    })
    return(as.vector(deletionCheck))
}
# Fix sequences with deletions
fixDeletions <- function(matInputP){
    delPos <- checkDeletions(matInputP)
    sapply((1:nrow(matInputP))[delPos],function(rowIndex){
        x <- matInputP[rowIndex,]
        inputGaps <- gregexpr("-",x[1])[[1]]
        glGaps <- gregexpr("-",x[2])[[1]]
        posDeletions <- inputGaps[!(inputGaps%in%glGaps)]
        inputDeletionToN <- s2c(x[1])
        inputDeletionToN[posDeletions] <- "N"
        inputDeletionToN <- c2s(inputDeletionToN)
        matInput[rowIndex,1] <<- inputDeletionToN
    })
    return(delPos)
}


# Trim DNA sequence to the last codon
trimToLastCodon <- function(seqToTrim){
    seqLen = nchar(seqToTrim)
    trimmedSeq = s2c(seqToTrim)
    poi = seqLen
    tailLen = 0
    
    while(trimmedSeq[poi]=="-" || trimmedSeq[poi]=="."){
        tailLen = tailLen + 1
        poi = poi - 1
    }
    
    trimmedSeq = c2s(trimmedSeq[1:(seqLen-tailLen)])
    seqLen = nchar(trimmedSeq)
    # Trim sequence to last codon
    if( getCodonPos(seqLen)[3] > seqLen )
        trimmedSeq = substr(seqToTrim,1, ( (getCodonPos(seqLen)[1])-1 ) )
    
    return(trimmedSeq)
}


# Given a nuclotide position, returns the codon number
# e.g. nuc 86  = codon 29
getCodonNumb <- function(nucPos){
    return( ceiling(nucPos/3) )
}

# Given a codon, returns all the nuc positions that make the codon
getCodonNucs <- function(codonNumb){
    getCodonPos(codonNumb*3)
}

# Given a nucleotide postions return the position in the codon
getContextInCodon <- function(nucPos){
    return( {nucPos-1}%%3+1 )
}

computeCodonTable <- function(testID=1){
    
    #if(testID<=5){
    # Pre-compute every codons
    intCounter = 1
    for(pOne in NUCLEOTIDES[1:4]){
        for(pTwo in NUCLEOTIDES[1:4]){
            for(pThree in NUCLEOTIDES[1:4]){
                codon = paste(pOne,pTwo,pThree,sep="")
                colnames(CODON_TABLE)[intCounter] =  codon
                intCounter = intCounter + 1
                CODON_TABLE[,codon] = mutationTypeOptimized(cbind(permutateAllCodon(codon),rep(codon,12)))
            }
        }
    }
    chars = c("N","A","C","G","T", "-")
    for(a in chars){
        for(b in chars){
            for(c in chars){
                if(a=="N" | b=="N" | c=="N"){
                    #cat(paste(a,b,c),sep="","\n")
                    CODON_TABLE[,paste(a,b,c,sep="")] = rep(NA,12)
                }
            }
        }
    }
    
    chars = c("-","A","C","G","T")
    for(a in chars){
        for(b in chars){
            for(c in chars){
                if(a=="-" | b=="-" | c=="-"){
                    #cat(paste(a,b,c),sep="","\n")
                    CODON_TABLE[,paste(a,b,c,sep="")] = rep(NA,12)
                }
            }
        }
    }
    CODON_TABLE <<- as.matrix(CODON_TABLE)
    #}
}

collapseCloneTry <- function(vecInputSeqs,glSeq,VLENGTH,nonTerminalOnly=0){
    tryCatch(return(collapseClone(vecInputSeqs,glSeq,VLENGTH,nonTerminalOnly)), error = function(e) {
        return(NA)
    })
}



# Compute the expected for each sequence-germline pair
getExpectedIndividual <- function(matInput){
    if( any(grep("multicore",search())) ){
        facGL <- factor(matInput[,2])
        facLevels = levels(facGL)
        LisGLs_MutabilityU = mclapply(1:length(facLevels),  function(x){
            computeMutabilities(facLevels[x])
        })
        facIndex = match(facGL,facLevels)
        
        LisGLs_Mutability = mclapply(1:nrow(matInput),  function(x){
            cInput = rep(NA,nchar(matInput[x,1]))
            cInput[s2c(matInput[x,1])!="N"] = 1
            LisGLs_MutabilityU[[facIndex[x]]] * cInput
        })
        
        LisGLs_Targeting =  mclapply(1:dim(matInput)[1],  function(x){
            computeTargeting(matInput[x,2],LisGLs_Mutability[[x]])
        })
        
        LisGLs_MutationTypes  = mclapply(1:length(matInput[,2]),function(x){
            #print(x)
            computeMutationTypes(matInput[x,2])
        })
        
        LisGLs_Exp = mclapply(1:dim(matInput)[1],  function(x){
            computeExpected(LisGLs_Targeting[[x]],LisGLs_MutationTypes[[x]])
        })
        
        ul_LisGLs_Exp =  unlist(LisGLs_Exp)
        return(matrix(ul_LisGLs_Exp,ncol=4,nrow=(length(ul_LisGLs_Exp)/4),byrow=T))
    }else{
        facGL <- factor(matInput[,2])
        facLevels = levels(facGL)
        LisGLs_MutabilityU = lapply(1:length(facLevels),  function(x){
            computeMutabilities(facLevels[x])
        })
        facIndex = match(facGL,facLevels)
        
        LisGLs_Mutability = lapply(1:nrow(matInput),  function(x){
            cInput = rep(NA,nchar(matInput[x,1]))
            cInput[s2c(matInput[x,1])!="N"] = 1
            LisGLs_MutabilityU[[facIndex[x]]] * cInput
        })
        
        LisGLs_Targeting =  lapply(1:dim(matInput)[1],  function(x){
            computeTargeting(matInput[x,2],LisGLs_Mutability[[x]])
        })
        
        LisGLs_MutationTypes  = lapply(1:length(matInput[,2]),function(x){
            #print(x)
            computeMutationTypes(matInput[x,2])
        })
        
        LisGLs_Exp = lapply(1:dim(matInput)[1],  function(x){
            computeExpected(LisGLs_Targeting[[x]],LisGLs_MutationTypes[[x]])
        })
        
        ul_LisGLs_Exp =  unlist(LisGLs_Exp)
        return(matrix(ul_LisGLs_Exp,ncol=4,nrow=(length(ul_LisGLs_Exp)/4),byrow=T))
        
    }
}


# Returns the mutability of a triplet at a given position
getMutability <- function(codon, pos=1:3){
    triplets <- rownames(HS5FModel@mutability)
    HS5FModel@mutability[  match(codon,triplets) ,pos]
}

getMutability5 <- function(fivemer){
    return(HS5FModel@mutability[fivemer])
}

# Returns the substitution probabilty
getTransistionProb <- function(nuc){
    HS5FModel@substitution[nuc,]
}

getTransistionProb5 <- function(fivemer){
    if(any(which(fivemer==colnames(HS5FModel@substitution)))){
        return(HS5FModel@substitution[,fivemer])
    }else{
        return(array(NA,4))
    }
}

# Given a nucleotide, returns the probabilty of other nucleotide it can mutate to
canMutateToProb <- function(nuc){
    HS5FModel@substitution[nuc, canMutateTo(nuc)]
}



# Compute the mutations types
computeMutationTypesFast <- function(param_strSeq){
    matMutationTypes = matrix( CODON_TABLE[,param_strSeq] ,ncol=3,nrow=4, byrow=F)
    if(mutabilityModel==5){
        matMutationTypes <- rbind(  matMutationTypes, matrix(NA,ncol=ncol(matMutationTypes),nrow=2))
    }
    #dimnames( matMutationTypes ) =  list(NUCLEOTIDES,1:(length(vecSeq)))
    return(matMutationTypes)
}


# Returns a vector of codons 1 mutation away from the given codon
permutateAllCodon <- function(codon){
    cCodon = s2c(codon)
    matCodons = t(array(cCodon,dim=c(3,12)))
    matCodons[1:4,1] = NUCLEOTIDES[1:4]
    matCodons[5:8,2] = NUCLEOTIDES[1:4]
    matCodons[9:12,3] = NUCLEOTIDES[1:4]
    apply(matCodons,1,c2s)
}


#given a mat of targeting & it's corresponding mutationtypes returns
#a vector of Exp_R frequecies in FWR1,CDR1,FWR2,CDR3,FWR3
computeExpected_Regions <- function(paramTargeting,paramMutationTypes){
    emptyArray  = array("FALSE",ncol(FWR_Nuc_Mat))
    
    FWR1= emptyArray
    FWR1[1:78]=T
    FWR1 = rbind(FWR1,FWR1,FWR1,FWR1,FWR1,FWR1)
    
    FWR2= emptyArray
    FWR2[115:165]=T
    FWR2 = rbind(FWR2,FWR2,FWR2,FWR2,FWR2,FWR2)
    
    FWR3= emptyArray
    FWR3[196:312]=T
    FWR3 = rbind(FWR3,FWR3,FWR3,FWR3,FWR3,FWR3)
    
    CDR1= emptyArray
    CDR1[79:114]=T
    CDR1 = rbind(CDR1,CDR1,CDR1,CDR1,CDR1,CDR1)
    
    CDR2= emptyArray
    CDR2[166:195]=T
    CDR2 = rbind(CDR2,CDR2,CDR2,CDR2,CDR2,CDR2)
    
    # Replacements
    RPos = which(paramMutationTypes=="R")
    Exp_R_FWR1 = sum(paramTargeting[ RPos[which(FWR1[RPos]==T)] ],na.rm=T)
    Exp_R_FWR2 = sum(paramTargeting[ RPos[which(FWR2[RPos]==T)] ],na.rm=T)
    Exp_R_FWR3 = sum(paramTargeting[ RPos[which(FWR3[RPos]==T)] ],na.rm=T)
    Exp_R_CDR1 = sum(paramTargeting[ RPos[which(CDR1[RPos]==T)] ],na.rm=T)
    Exp_R_CDR2 = sum(paramTargeting[ RPos[which(CDR2[RPos]==T)] ],na.rm=T)
    
    SPos = which(paramMutationTypes!="R")
    Exp_S_FWR1 = sum(paramTargeting[ SPos[which(FWR1[SPos]==T)] ],na.rm=T)
    Exp_S_FWR2 = sum(paramTargeting[ SPos[which(FWR2[SPos]==T)] ],na.rm=T)
    Exp_S_FWR3 = sum(paramTargeting[ SPos[which(FWR3[SPos]==T)] ],na.rm=T)
    Exp_S_CDR1 = sum(paramTargeting[ SPos[which(CDR1[SPos]==T)] ],na.rm=T)
    Exp_S_CDR2 = sum(paramTargeting[ SPos[which(CDR2[SPos]==T)] ],na.rm=T)
    
    return(c(Exp_R_FWR1/(Exp_R_FWR1+Exp_S_FWR1),
             Exp_R_CDR1/(Exp_R_CDR1+Exp_S_CDR1),
             Exp_R_FWR2/(Exp_R_FWR2+Exp_S_FWR2),
             Exp_R_CDR2/(Exp_R_CDR2+Exp_S_CDR2),
             Exp_R_FWR3/(Exp_R_FWR3+Exp_S_FWR3)
    ))
}


processNucMutations2 <- function(mu){
    if(!is.na(mu)){
        #R
        if(any(mu=="R")){
            Rs = mu[mu=="R"]
            nucNumbs = as.numeric(names(Rs))
            R_CDR = sum(as.integer(CDR_Nuc[nucNumbs]),na.rm=T)
            R_FWR = sum(as.integer(FWR_Nuc[nucNumbs]),na.rm=T)
        }else{
            R_CDR = 0
            R_FWR = 0
        }
        
        #S
        if(any(mu=="S")){
            Ss = mu[mu=="S"]
            nucNumbs = as.numeric(names(Ss))
            S_CDR = sum(as.integer(CDR_Nuc[nucNumbs]),na.rm=T)
            S_FWR = sum(as.integer(FWR_Nuc[nucNumbs]),na.rm=T)
        }else{
            S_CDR = 0
            S_FWR = 0
        }
        
        
        retVec = c(R_CDR,S_CDR,R_FWR,S_FWR)
        retVec[is.na(retVec)]=0
        return(retVec)
    }else{
        return(rep(0,4))
    }
}


## Z-score Test
computeZScore <- function(mat, test="Focused"){
    matRes <- matrix(NA,ncol=2,nrow=(nrow(mat)))
    if(test=="Focused"){
        #Z_Focused_CDR
        #P_Denom = sum( mat[1,c(5,6,8)], na.rm=T )
        P = apply(mat[,c(5,6,8)],1,function(x){(x[1]/sum(x))})
        R_mean = apply(cbind(mat[,c(1,2,4)],P),1,function(x){x[4]*(sum(x[1:3]))})
        R_sd=sqrt(R_mean*(1-P))
        matRes[,1] = (mat[,1]-R_mean)/R_sd
        
        #Z_Focused_FWR
        #P_Denom = sum( mat[1,c(7,6,8)], na.rm=T )
        P = apply(mat[,c(7,6,8)],1,function(x){(x[1]/sum(x))})
        R_mean = apply(cbind(mat[,c(3,2,4)],P),1,function(x){x[4]*(sum(x[1:3]))})
        R_sd=sqrt(R_mean*(1-P))
        matRes[,2] = (mat[,3]-R_mean)/R_sd
    }
    
    if(test=="Local"){
        #Z_Focused_CDR
        #P_Denom = sum( mat[1,c(5,6,8)], na.rm=T )
        P = apply(mat[,c(5,6)],1,function(x){(x[1]/sum(x))})
        R_mean = apply(cbind(mat[,c(1,2)],P),1,function(x){x[3]*(sum(x[1:2]))})
        R_sd=sqrt(R_mean*(1-P))
        matRes[,1] = (mat[,1]-R_mean)/R_sd
        
        #Z_Focused_FWR
        #P_Denom = sum( mat[1,c(7,6,8)], na.rm=T )
        P = apply(mat[,c(7,8)],1,function(x){(x[1]/sum(x))})
        R_mean = apply(cbind(mat[,c(3,4)],P),1,function(x){x[3]*(sum(x[1:2]))})
        R_sd=sqrt(R_mean*(1-P))
        matRes[,2] = (mat[,3]-R_mean)/R_sd
    }
    
    if(test=="Imbalanced"){
        #Z_Focused_CDR
        #P_Denom = sum( mat[1,c(5,6,8)], na.rm=T )
        P = apply(mat[,5:8],1,function(x){((x[1]+x[2])/sum(x))})
        R_mean = apply(cbind(mat[,1:4],P),1,function(x){x[5]*(sum(x[1:4]))})
        R_sd=sqrt(R_mean*(1-P))
        matRes[,1] = (mat[,1]-R_mean)/R_sd
        
        #Z_Focused_FWR
        #P_Denom = sum( mat[1,c(7,6,8)], na.rm=T )
        P = apply(mat[,5:8],1,function(x){((x[3]+x[4])/sum(x))})
        R_mean = apply(cbind(mat[,1:4],P),1,function(x){x[5]*(sum(x[1:4]))})
        R_sd=sqrt(R_mean*(1-P))
        matRes[,2] = (mat[,3]-R_mean)/R_sd
    }
    
    matRes[is.nan(matRes)] = NA
    return(matRes)
}

# Return a p-value for a z-score
z2p <- function(z){
    p=NA
    if( !is.nan(z) && !is.na(z)){
        if(z>0){
            p = (1 - pnorm(z,0,1))
        } else if(z<0){
            p = (-1 * pnorm(z,0,1))
        } else{
            p = 0.5
        }
    }else{
        p = NA
    }
    return(p)
}


## Bayesian  Test

# Fitted parameter for the bayesian framework
BAYESIAN_FITTED<-c(0.407277142798302, 0.554007336744485, 0.63777155771234, 0.693989162719009, 0.735450014674917, 0.767972534429806, 0.794557287143399, 0.816906816601605, 0.83606796225341, 0.852729446430296, 0.867370424541641, 0.880339760590323, 0.891900995024999, 0.902259181289864, 0.911577919359,0.919990301665853, 0.927606458124537, 0.934518806350661, 0.940805863754375, 0.946534836475715, 0.951763691199255, 0.95654428191308, 0.960920179487397, 0.964930893680829, 0.968611312149038, 0.971992459313836, 0.975102110004818, 0.977964943023096, 0.980603428208439, 0.983037660179428, 0.985285800977406, 0.987364285326685, 0.989288037855441, 0.991070478823525, 0.992723699729969, 0.994259575477392, 0.995687688867975, 0.997017365051493, 0.998257085153047, 0.999414558305388, 1.00049681357804, 1.00151036237481, 1.00246080204981, 1.00335370751909, 1.0041939329768, 1.0049859393417, 1.00573382091263, 1.00644127217376, 1.00711179729107, 1.00774845526417, 1.00835412715854, 1.00893143010366, 1.00948275846309, 1.01001030293661, 1.01051606798079, 1.01100188771288, 1.01146944044216, 1.01192026195449, 1.01235575766094, 1.01277721370986)
CONST_i <- sort(c(((2^(seq(-39,0,length.out=201)))/2)[1:200],(c(0:11,13:99)+0.5)/100,1-(2^(seq(-39,0,length.out=201)))/2))

# Given x, M & p, returns a pdf
# calculate_bayes <- function ( x=3, N=10, p=0.33,
#                               i=CONST_i,
#                               max_sigma=20,length_sigma=4001
# ){
#   if(!0%in%N){
#     G <- max(length(x),length(N),length(p))
#     x=array(x,dim=G)
#     N=array(N,dim=G)
#     p=array(p,dim=G)
#     sigma_s<-seq(-max_sigma,max_sigma,length.out=length_sigma)
#     sigma_1<-log({i/{1-i}}/{p/{1-p}})
#     index<-min(N,60)
#     y<-dbeta(i,x+BAYESIAN_FITTED[index],N+BAYESIAN_FITTED[index]-x)*(1-p)*p*exp(sigma_1)/({1-p}^2+2*p*{1-p}*exp(sigma_1)+{p^2}*exp(2*sigma_1))
#     if(!sum(is.na(y))){
#       tmp<-approx(sigma_1,y,sigma_s)$y
#       tmp/sum(tmp)/{2*max_sigma/{length_sigma-1}}
#     }else{
#       return(NA)
#     }
#   }else{
#     return(NA)
#   }
# }
calculate_bayes <- function ( x=3, n=10, p=0.33,
                              i=CONST_i,
                              max_sigma=20,length_sigma=4001){
    if(n!=0){
        sigma_s<-seq(-max_sigma,max_sigma,length.out=length_sigma)
        sigma_1<-log({i/{1-i}}/{p/{1-p}})
        index<-min(n,60)
        y<- dbeta(i,x+BAYESIAN_FITTED[index],n+BAYESIAN_FITTED[index]-x)*(1-p)*p*exp(sigma_1)/({1-p}^2+2*p*{1-p}*exp(sigma_1)+{p^2}*exp(2*sigma_1))
        if(!sum(is.na(y))){
            tmp<-approx(sigma_1,y,sigma_s)$y
            tmp/sum(tmp)/{2*max_sigma/{length_sigma-1}}
        }else{
            return(NA)
        }
    }else{
        return(NA)
    }
}
# Given a mat of observed & expected, return a list of CDR & FWR pdf for selection
computeBayesianScore <- function(mat, test="Focused", max_sigma=20,length_sigma=4001){
    flagOneSeq = F
    if(nrow(mat)==1){
        mat=rbind(mat,mat)
        flagOneSeq = T
    }
    if(test=="Focused"){
        #CDR
        P = c(apply(mat[,c(5,6,8)],1,function(x){(x[1]/sum(x))}),0.5)
        N = c(apply(mat[,c(1,2,4)],1,function(x){(sum(x))}),0)
        X = c(mat[,1],0)
        bayesCDR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})
        bayesCDR = bayesCDR[-length(bayesCDR)]
        
        #FWR
        P = c(apply(mat[,c(7,6,8)],1,function(x){(x[1]/sum(x))}),0.5)
        N = c(apply(mat[,c(3,2,4)],1,function(x){(sum(x))}),0)
        X = c(mat[,3],0)
        bayesFWR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})
        bayesFWR = bayesFWR[-length(bayesFWR)]
    }
    
    if(test=="Local"){
        #CDR
        P = c(apply(mat[,c(5,6)],1,function(x){(x[1]/sum(x))}),0.5)
        N = c(apply(mat[,c(1,2)],1,function(x){(sum(x))}),0)
        X = c(mat[,1],0)
        bayesCDR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})
        bayesCDR = bayesCDR[-length(bayesCDR)]
        
        #FWR
        P = c(apply(mat[,c(7,8)],1,function(x){(x[1]/sum(x))}),0.5)
        N = c(apply(mat[,c(3,4)],1,function(x){(sum(x))}),0)
        X = c(mat[,3],0)
        bayesFWR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})
        bayesFWR = bayesFWR[-length(bayesFWR)]
    }
    
    if(test=="Imbalanced"){
        #CDR
        P = c(apply(mat[,c(5:8)],1,function(x){((x[1]+x[2])/sum(x))}),0.5)
        N = c(apply(mat[,c(1:4)],1,function(x){(sum(x))}),0)
        X = c(apply(mat[,c(1:2)],1,function(x){(sum(x))}),0)
        bayesCDR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})
        bayesCDR = bayesCDR[-length(bayesCDR)]
        
        #FWR
        P = c(apply(mat[,c(5:8)],1,function(x){((x[3]+x[4])/sum(x))}),0.5)
        N = c(apply(mat[,c(1:4)],1,function(x){(sum(x))}),0)
        X = c(apply(mat[,c(3:4)],1,function(x){(sum(x))}),0)
        bayesFWR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})
        bayesFWR = bayesFWR[-length(bayesFWR)]
    }
    
    if(test=="ImbalancedSilent"){
        #CDR
        P = c(apply(mat[,c(6,8)],1,function(x){((x[1])/sum(x))}),0.5)
        N = c(apply(mat[,c(2,4)],1,function(x){(sum(x))}),0)
        X = c(apply(mat[,c(2,4)],1,function(x){(x[1])}),0)
        bayesCDR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})
        bayesCDR = bayesCDR[-length(bayesCDR)]
        
        #FWR
        P = c(apply(mat[,c(6,8)],1,function(x){((x[2])/sum(x))}),0.5)
        N = c(apply(mat[,c(2,4)],1,function(x){(sum(x))}),0)
        X = c(apply(mat[,c(2,4)],1,function(x){(x[2])}),0)
        bayesFWR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})
        bayesFWR = bayesFWR[-length(bayesFWR)]
    }
    
    if(flagOneSeq==T){
        bayesCDR = bayesCDR[1]
        bayesFWR = bayesFWR[1]
    }
    return( list("CDR"=bayesCDR, "FWR"=bayesFWR) )
}


# Computes the 95% CI for a pdf
# calcBayesCI <- function(Pdf,low=0.025,up=0.975,max_sigma=20, length_sigma=4001){
#   if(length(Pdf)!=length_sigma) return(NA)
#   sigma_s=seq(-max_sigma,max_sigma,length.out=length_sigma)
#   cdf = cumsum(Pdf)
#   cdf = cdf/cdf[length(cdf)]
#   return( c(sigma_s[findInterval(low,cdf)-1] , sigma_s[findInterval(up,cdf)]) )
# }
calcBayesCI <- function(Pdf,low=0.025,up=0.975,max_sigma=20, length_sigma=4001){
    if(length(Pdf)!=length_sigma) return(NA)
    sigma_s=seq(-max_sigma,max_sigma,length.out=length_sigma)
    cdf = cumsum(Pdf)
    cdf = cdf/cdf[length(cdf)]
    intervalLow = findInterval(low,cdf)
    fractionLow = (low - cdf[intervalLow])/(cdf[intervalLow+1]-cdf[intervalLow])
    intervalUp = findInterval(up,cdf)
    fractionUp = (up - cdf[intervalUp])/(cdf[intervalUp]-cdf[intervalUp-1])
    sigmaLow = sigma_s[intervalLow]+fractionLow*(sigma_s[intervalLow+1]-sigma_s[intervalLow])
    sigmaUp = sigma_s[intervalUp]+fractionUp*(sigma_s[intervalUp+1]-sigma_s[intervalUp])
    return( c(sigmaLow,sigmaUp) )
}

# Computes a mean for a pdf
calcBayesMean <- function(Pdf,max_sigma=20,length_sigma=4001){
    if(length(Pdf)!=length_sigma) return(NA)
    sigma_s=seq(-max_sigma,max_sigma,length.out=length_sigma)
    norm = {length_sigma-1}/2/max_sigma
    return( (Pdf%*%sigma_s/norm)  )
}

# Returns the mean, and the 95% CI for a pdf
calcBayesOutputInfo <- function(Pdf,low=0.025,up=0.975,max_sigma=20, length_sigma=4001){
    if(is.na(Pdf))
        return(rep(NA,3))
    bCI = calcBayesCI(Pdf=Pdf,low=low,up=up,max_sigma=max_sigma,length_sigma=length_sigma)
    bMean = calcBayesMean(Pdf=Pdf,max_sigma=max_sigma,length_sigma=length_sigma)
    return(c(bMean, bCI))
}

# Computes the p-value of a pdf
computeSigmaP <- function(Pdf, length_sigma=4001, max_sigma=20){
    if(length(Pdf)>1){
        norm = {length_sigma-1}/2/max_sigma
        pVal = {sum(Pdf[1:{{length_sigma-1}/2}]) + Pdf[{{length_sigma+1}/2}]/2}/norm
        if(pVal>0.5){
            pVal = pVal-1
        }
        return(pVal)
    }else{
        return(NA)
    }
}

# Compute p-value of two distributions
compareTwoDistsFaster <-function(sigma_S=seq(-20,20,length.out=4001), N=10000, dens1=runif(4001,0,1), dens2=runif(4001,0,1)){
    #print(c(length(dens1),length(dens2)))
    if(length(dens1)>1 & length(dens2)>1 ){
        dens1<-dens1/sum(dens1)
        dens2<-dens2/sum(dens2)
        cum2 <- cumsum(dens2)-dens2/2
        tmp<- sum(sapply(1:length(dens1),function(i)return(dens1[i]*cum2[i])))
        #print(tmp)
        if(tmp>0.5)tmp<-tmp-1
        return( tmp )
    }
    else {
        return(NA)
    }
    #return (sum(sapply(1:N,function(i)(sample(sigma_S,1,prob=dens1)>sample(sigma_S,1,prob=dens2))))/N)
}

# get number of seqeunces contributing to the sigma (i.e. seqeunces with mutations)
numberOfSeqsWithMutations <- function(matMutations,test=1){
    if(test==4)test=2
    cdrSeqs <- 0
    fwrSeqs <- 0
    if(test==1){#focused
        cdrMutations <- apply(matMutations, 1, function(x){ sum(x[c(1,2,4)]) })
        fwrMutations <- apply(matMutations, 1, function(x){ sum(x[c(3,4,2)]) })
        if( any(which(cdrMutations>0)) ) cdrSeqs <- sum(cdrMutations>0)
        if( any(which(fwrMutations>0)) ) fwrSeqs <- sum(fwrMutations>0)
    }
    if(test==2){#local
        cdrMutations <- apply(matMutations, 1, function(x){ sum(x[c(1,2)]) })
        fwrMutations <- apply(matMutations, 1, function(x){ sum(x[c(3,4)]) })
        if( any(which(cdrMutations>0)) ) cdrSeqs <- sum(cdrMutations>0)
        if( any(which(fwrMutations>0)) ) fwrSeqs <- sum(fwrMutations>0)
    }
    return(c("CDR"=cdrSeqs, "FWR"=fwrSeqs))
}



shadeColor <- function(sigmaVal=NA,pVal=NA){
    if(is.na(sigmaVal) & is.na(pVal)) return(NA)
    if(is.na(sigmaVal) & !is.na(pVal)) sigmaVal=sign(pVal)
    if(is.na(pVal) || pVal==1 || pVal==0){
        returnColor = "#FFFFFF";
    }else{
        colVal=abs(pVal);
        
        if(sigmaVal<0){
            if(colVal>0.1)
                returnColor = "#CCFFCC";
            if(colVal<=0.1)
                returnColor = "#99FF99";
            if(colVal<=0.050)
                returnColor = "#66FF66";
            if(colVal<=0.010)
                returnColor = "#33FF33";
            if(colVal<=0.005)
                returnColor = "#00FF00";
            
        }else{
            if(colVal>0.1)
                returnColor = "#FFCCCC";
            if(colVal<=0.1)
                returnColor = "#FF9999";
            if(colVal<=0.05)
                returnColor = "#FF6666";
            if(colVal<=0.01)
                returnColor = "#FF3333";
            if(colVal<0.005)
                returnColor = "#FF0000";
        }
    }
    
    return(returnColor)
}



plotHelp <- function(xfrac=0.05,yfrac=0.05,log=FALSE){
    if(!log){
        x = par()$usr[1]-(par()$usr[2]-par()$usr[1])*xfrac
        y = par()$usr[4]+(par()$usr[4]-par()$usr[3])*yfrac
    }else {
        if(log==2){
            x = par()$usr[1]-(par()$usr[2]-par()$usr[1])*xfrac
            y = 10^((par()$usr[4])+((par()$usr[4])-(par()$usr[3]))*yfrac)
        }
        if(log==1){
            x = 10^((par()$usr[1])-((par()$usr[2])-(par()$usr[1]))*xfrac)
            y = par()$usr[4]+(par()$usr[4]-par()$usr[3])*yfrac
        }
        if(log==3){
            x = 10^((par()$usr[1])-((par()$usr[2])-(par()$usr[1]))*xfrac)
            y = 10^((par()$usr[4])+((par()$usr[4])-(par()$usr[3]))*yfrac)
        }
    }
    return(c("x"=x,"y"=y))
}

# SHMulation

# Based on targeting, introduce a single mutation & then update the targeting
oneMutation <- function(){
    # Pick a postion + mutation
    posMutation = sample(1:(seqGermlineLen*4),1,replace=F,prob=as.vector(seqTargeting))
    posNucNumb = ceiling(posMutation/4)                    # Nucleotide number
    posNucKind = 4 - ( (posNucNumb*4) - posMutation )   # Nuc the position mutates to
    
    #mutate the simulation sequence
    seqSimVec <-  s2c(seqSim)
    seqSimVec[posNucNumb] <- NUCLEOTIDES[posNucKind]
    seqSim <<-  c2s(seqSimVec)
    
    #update Mutability, Targeting & MutationsTypes
    updateMutabilityNTargeting(posNucNumb)
    
    #return(c(posNucNumb,NUCLEOTIDES[posNucKind]))
    return(posNucNumb)
}

updateMutabilityNTargeting <- function(position){
    min_i<-max((position-2),1)
    max_i<-min((position+2),nchar(seqSim))
    min_ii<-min(min_i,3)
    
    #mutability - update locally
    seqMutability[(min_i):(max_i)] <<- computeMutabilities(substr(seqSim,position-4,position+4))[(min_ii):(max_i-min_i+min_ii)]
    
    
    #targeting - compute locally
    seqTargeting[,min_i:max_i] <<- computeTargeting(substr(seqSim,min_i,max_i),seqMutability[min_i:max_i])
    seqTargeting[is.na(seqTargeting)] <<- 0
    #mutCodonPos = getCodonPos(position)
    mutCodonPos = seq(getCodonPos(min_i)[1],getCodonPos(max_i)[3])
    #cat(mutCodonPos,"\n")
    mutTypeCodon = getCodonPos(position)
    seqMutationTypes[,mutTypeCodon] <<- computeMutationTypesFast( substr(seqSim,mutTypeCodon[1],mutTypeCodon[3]) )
    # Stop = 0
    if(any(seqMutationTypes[,mutCodonPos]=="Stop",na.rm=T )){
        seqTargeting[,mutCodonPos][seqMutationTypes[,mutCodonPos]=="Stop"] <<- 0
    }
    
    
    #Selection
    selectedPos = (min_i*4-4)+(which(seqMutationTypes[,min_i:max_i]=="R"))
    # CDR
    selectedCDR = selectedPos[which(matCDR[selectedPos]==T)]
    seqTargeting[selectedCDR] <<-  seqTargeting[selectedCDR] *  exp(selCDR)
    seqTargeting[selectedCDR] <<- seqTargeting[selectedCDR]/baseLineCDR_K
    
    # FWR
    selectedFWR = selectedPos[which(matFWR[selectedPos]==T)]
    seqTargeting[selectedFWR] <<-  seqTargeting[selectedFWR] *  exp(selFWR)
    seqTargeting[selectedFWR] <<- seqTargeting[selectedFWR]/baseLineFWR_K
    
}



# Validate the mutation: if the mutation has not been sampled before validate it, else discard it.
validateMutation <- function(){
    if( !(mutatedPos%in%mutatedPositions) ){ # if it's a new mutation
        uniqueMutationsIntroduced <<- uniqueMutationsIntroduced + 1
        mutatedPositions[uniqueMutationsIntroduced] <<-  mutatedPos
    }else{
        if(substr(seqSim,mutatedPos,mutatedPos)==substr(seqGermline,mutatedPos,mutatedPos)){ # back to germline mutation
            mutatedPositions <<-  mutatedPositions[-which(mutatedPositions==mutatedPos)]
            uniqueMutationsIntroduced <<-  uniqueMutationsIntroduced - 1
        }
    }
}



# Places text (labels) at normalized coordinates
myaxis <- function(xfrac=0.05,yfrac=0.05,log=FALSE,w="text",cex=1,adj=1,thecol="black"){
    par(xpd=TRUE)
    if(!log)
        text(par()$usr[1]-(par()$usr[2]-par()$usr[1])*xfrac,par()$usr[4]+(par()$usr[4]-par()$usr[3])*yfrac,w,cex=cex,adj=adj,col=thecol)
    else {
        if(log==2)
            text(
                par()$usr[1]-(par()$usr[2]-par()$usr[1])*xfrac,
                10^((par()$usr[4])+((par()$usr[4])-(par()$usr[3]))*yfrac),
                w,cex=cex,adj=adj,col=thecol)
        if(log==1)
            text(
                10^((par()$usr[1])-((par()$usr[2])-(par()$usr[1]))*xfrac),
                par()$usr[4]+(par()$usr[4]-par()$usr[3])*yfrac,
                w,cex=cex,adj=adj,col=thecol)
        if(log==3)
            text(
                10^((par()$usr[1])-((par()$usr[2])-(par()$usr[1]))*xfrac),
                10^((par()$usr[4])+((par()$usr[4])-(par()$usr[3]))*yfrac),
                w,cex=cex,adj=adj,col=thecol)
    }
    par(xpd=FALSE)
}



# Count the mutations in a sequence
analyzeMutations <- function( inputMatrixIndex, model = 0 , multipleMutation=0, seqWithStops=0){
    
    paramGL = s2c(matInput[inputMatrixIndex,2])
    paramSeq = s2c(matInput[inputMatrixIndex,1])
    
    #if( any(paramSeq=="N") ){
    #  gapPos_Seq =  which(paramSeq=="N")
    #  gapPos_Seq_ToReplace = gapPos_Seq[paramGL[gapPos_Seq] != "N"]
    #  paramSeq[gapPos_Seq_ToReplace] =  paramGL[gapPos_Seq_ToReplace]
    #}
    mutations_val = paramGL != paramSeq
    
    if(any(mutations_val)){
        mutationPos = which(mutations_val)#{1:length(mutations_val)}[mutations_val]
        length_mutations =length(mutationPos)
        mutationInfo = rep(NA,length_mutations)
        
        pos<- mutationPos
        pos_array<-array(sapply(pos,getCodonPos))
        codonGL =  paramGL[pos_array]
        codonSeqWhole =  paramSeq[pos_array]
        codonSeq = sapply(pos,function(x){
            seqP = paramGL[getCodonPos(x)]
            muCodonPos = {x-1}%%3+1
            seqP[muCodonPos] = paramSeq[x]
            return(seqP)
        })
        GLcodons =  apply(matrix(codonGL,length_mutations,3,byrow=TRUE),1,c2s)
        SeqcodonsWhole =  apply(matrix(codonSeqWhole,length_mutations,3,byrow=TRUE),1,c2s)
        Seqcodons =   apply(codonSeq,2,c2s)
        
        mutationInfo = apply(rbind(GLcodons , Seqcodons),2,function(x){mutationType(c2s(x[1]),c2s(x[2]))})
        names(mutationInfo) = mutationPos
        
        mutationInfoWhole = apply(rbind(GLcodons , SeqcodonsWhole),2,function(x){mutationType(c2s(x[1]),c2s(x[2]))})
        names(mutationInfoWhole) = mutationPos
        
        mutationInfo <- mutationInfo[!is.na(mutationInfo)]
        mutationInfoWhole <- mutationInfoWhole[!is.na(mutationInfoWhole)]
        
        if(any(!is.na(mutationInfo))){
            
            #Filter based on Stop (at the codon level)
            if(seqWithStops==1){
                nucleotidesAtStopCodons = names(mutationInfoWhole[mutationInfoWhole!="Stop"])
                mutationInfo = mutationInfo[nucleotidesAtStopCodons]
                mutationInfoWhole = mutationInfo[nucleotidesAtStopCodons]
            }else{
                countStops = sum(mutationInfoWhole=="Stop")
                if(seqWithStops==2 & countStops==0) mutationInfo = NA
                if(seqWithStops==3 & countStops>0) mutationInfo = NA
            }
            
            if(any(!is.na(mutationInfo))){
                #Filter mutations based on multipleMutation
                if(multipleMutation==1 & !is.na(mutationInfo)){
                    mutationCodons = getCodonNumb(as.numeric(names(mutationInfoWhole)))
                    tableMutationCodons <- table(mutationCodons)
                    codonsWithMultipleMutations <- as.numeric(names(tableMutationCodons[tableMutationCodons>1]))
                    if(any(codonsWithMultipleMutations)){
                        #remove the nucleotide mutations in the codons with multiple mutations
                        mutationInfo <- mutationInfo[!(mutationCodons %in% codonsWithMultipleMutations)]
                        #replace those codons with Ns in the input sequence
                        paramSeq[unlist(lapply(codonsWithMultipleMutations, getCodonNucs))] = "N"
                        matInput[inputMatrixIndex,1] <<- c2s(paramSeq)
                    }
                }
                
                #Filter mutations based on the model
                if(any(mutationInfo)==T | is.na(any(mutationInfo))){
                    
                    if(model==1 & !is.na(mutationInfo)){
                        mutationInfo <- mutationInfo[mutationInfo=="S"]
                    }
                    if(any(mutationInfo)==T | is.na(any(mutationInfo))) return(mutationInfo)
                    else return(NA)
                }else{
                    return(NA)
                }
            }else{
                return(NA)
            }
            
            
        }else{
            return(NA)
        }
        
        
    }else{
        return (NA)
    }
}

analyzeMutationsFixed <- function( inputArray, model = 0 , multipleMutation=0, seqWithStops=0){
    
    paramGL = s2c(inputArray[2])
    paramSeq = s2c(inputArray[1])
    inputSeq <- inputArray[1]
    #if( any(paramSeq=="N") ){
    #  gapPos_Seq =  which(paramSeq=="N")
    #  gapPos_Seq_ToReplace = gapPos_Seq[paramGL[gapPos_Seq] != "N"]
    #  paramSeq[gapPos_Seq_ToReplace] =  paramGL[gapPos_Seq_ToReplace]
    #}
    mutations_val = paramGL != paramSeq
    
    if(any(mutations_val)){
        mutationPos = which(mutations_val)#{1:length(mutations_val)}[mutations_val]
        length_mutations =length(mutationPos)
        mutationInfo = rep(NA,length_mutations)
        
        pos<- mutationPos
        pos_array<-array(sapply(pos,getCodonPos))
        codonGL =  paramGL[pos_array]
        codonSeqWhole =  paramSeq[pos_array]
        codonSeq = sapply(pos,function(x){
            seqP = paramGL[getCodonPos(x)]
            muCodonPos = {x-1}%%3+1
            seqP[muCodonPos] = paramSeq[x]
            return(seqP)
        })
        GLcodons =  apply(matrix(codonGL,length_mutations,3,byrow=TRUE),1,c2s)
        SeqcodonsWhole =  apply(matrix(codonSeqWhole,length_mutations,3,byrow=TRUE),1,c2s)
        Seqcodons =   apply(codonSeq,2,c2s)
        
        mutationInfo = apply(rbind(GLcodons , Seqcodons),2,function(x){mutationType(c2s(x[1]),c2s(x[2]))})
        names(mutationInfo) = mutationPos
        
        mutationInfoWhole = apply(rbind(GLcodons , SeqcodonsWhole),2,function(x){mutationType(c2s(x[1]),c2s(x[2]))})
        names(mutationInfoWhole) = mutationPos
        
        mutationInfo <- mutationInfo[!is.na(mutationInfo)]
        mutationInfoWhole <- mutationInfoWhole[!is.na(mutationInfoWhole)]
        
        if(any(!is.na(mutationInfo))){
            
            #Filter based on Stop (at the codon level)
            if(seqWithStops==1){
                nucleotidesAtStopCodons = names(mutationInfoWhole[mutationInfoWhole!="Stop"])
                mutationInfo = mutationInfo[nucleotidesAtStopCodons]
                mutationInfoWhole = mutationInfo[nucleotidesAtStopCodons]
            }else{
                countStops = sum(mutationInfoWhole=="Stop")
                if(seqWithStops==2 & countStops==0) mutationInfo = NA
                if(seqWithStops==3 & countStops>0) mutationInfo = NA
            }
            
            if(any(!is.na(mutationInfo))){
                #Filter mutations based on multipleMutation
                if(multipleMutation==1 & !is.na(mutationInfo)){
                    mutationCodons = getCodonNumb(as.numeric(names(mutationInfoWhole)))
                    tableMutationCodons <- table(mutationCodons)
                    codonsWithMultipleMutations <- as.numeric(names(tableMutationCodons[tableMutationCodons>1]))
                    if(any(codonsWithMultipleMutations)){
                        #remove the nucleotide mutations in the codons with multiple mutations
                        mutationInfo <- mutationInfo[!(mutationCodons %in% codonsWithMultipleMutations)]
                        #replace those codons with Ns in the input sequence
                        paramSeq[unlist(lapply(codonsWithMultipleMutations, getCodonNucs))] = "N"
                        #matInput[inputMatrixIndex,1] <<- c2s(paramSeq)
                        inputSeq <- c2s(paramSeq)
                    }
                }
                
                #Filter mutations based on the model
                if(any(mutationInfo)==T | is.na(any(mutationInfo))){
                    
                    if(model==1 & !is.na(mutationInfo)){
                        mutationInfo <- mutationInfo[mutationInfo=="S"]
                    }
                    if(any(mutationInfo)==T | is.na(any(mutationInfo))) return(list(mutationInfo,inputSeq))
                    else return(list(NA,inputSeq))
                }else{
                    return(list(NA,inputSeq))
                }
            }else{
                return(list(NA,inputSeq))
            }
            
            
        }else{
            return(list(NA,inputSeq))
        }
        
        
    }else{
        return (list(NA,inputSeq))
    }
}

# triMutability Background Count
buildMutabilityModel <- function( inputMatrixIndex, model=0 , multipleMutation=0, seqWithStops=0, stopMutations=0){
    
    #rowOrigMatInput = matInput[inputMatrixIndex,]
    seqGL =  gsub("-", "", matInput[inputMatrixIndex,2])
    seqInput = gsub("-", "", matInput[inputMatrixIndex,1])
    #matInput[inputMatrixIndex,] <<- cbind(seqInput,seqGL)
    tempInput <- cbind(seqInput,seqGL)
    seqLength = nchar(seqGL)
    list_analyzeMutationsFixed<- analyzeMutationsFixed(tempInput, model, multipleMutation, seqWithStops)
    mutationCount <- list_analyzeMutationsFixed[[1]]
    seqInput <- list_analyzeMutationsFixed[[2]]
    BackgroundMatrix = mutabilityMatrix
    MutationMatrix = mutabilityMatrix
    MutationCountMatrix = mutabilityMatrix
    if(!is.na(mutationCount)){
        if((stopMutations==0 & model==0) | (stopMutations==1 & (sum(mutationCount=="Stop")<length(mutationCount))) | (model==1 & (sum(mutationCount=="S")>0)) ){
            
            fivermerStartPos = 1:(seqLength-4)
            fivemerLength <- length(fivermerStartPos)
            fivemerGL <- substr(rep(seqGL,length(fivermerStartPos)),(fivermerStartPos),(fivermerStartPos+4))
            fivemerSeq <- substr(rep(seqInput,length(fivermerStartPos)),(fivermerStartPos),(fivermerStartPos+4))
            
            #Background
            for(fivemerIndex in 1:fivemerLength){
                fivemer = fivemerGL[fivemerIndex]
                if(!any(grep("N",fivemer))){
                    fivemerCodonPos = fivemerCodon(fivemerIndex)
                    fivemerReadingFrameCodon = substr(fivemer,fivemerCodonPos[1],fivemerCodonPos[3])
                    fivemerReadingFrameCodonInputSeq = substr(fivemerSeq[fivemerIndex],fivemerCodonPos[1],fivemerCodonPos[3])
                    
                    # All mutations model
                    #if(!any(grep("N",fivemerReadingFrameCodon))){
                    if(model==0){
                        if(stopMutations==0){
                            if(!any(grep("N",fivemerReadingFrameCodonInputSeq)))
                                BackgroundMatrix[fivemer] <- (BackgroundMatrix[fivemer] + 1)
                        }else{
                            if( !any(grep("N",fivemerReadingFrameCodonInputSeq)) & translateCodonToAminoAcid(fivemerReadingFrameCodon)!="*" ){
                                positionWithinCodon = which(fivemerCodonPos==3)#positionsWithinCodon[(fivemerCodonPos[1]%%3)+1]
                                BackgroundMatrix[fivemer] <- (BackgroundMatrix[fivemer] + probNonStopMutations[fivemerReadingFrameCodon,positionWithinCodon])
                            }
                        }
                    }else{ # Only silent mutations
                        if( !any(grep("N",fivemerReadingFrameCodonInputSeq)) & translateCodonToAminoAcid(fivemerReadingFrameCodon)!="*" & translateCodonToAminoAcid(fivemerReadingFrameCodonInputSeq)==translateCodonToAminoAcid(fivemerReadingFrameCodon)){
                            positionWithinCodon = which(fivemerCodonPos==3)
                            BackgroundMatrix[fivemer] <- (BackgroundMatrix[fivemer] + probSMutations[fivemerReadingFrameCodon,positionWithinCodon])
                        }
                    }
                    #}
                }
            }
            
            #Mutations
            if(stopMutations==1) mutationCount = mutationCount[mutationCount!="Stop"]
            if(model==1) mutationCount = mutationCount[mutationCount=="S"]
            mutationPositions = as.numeric(names(mutationCount))
            mutationCount = mutationCount[mutationPositions>2 & mutationPositions<(seqLength-1)]
            mutationPositions =  mutationPositions[mutationPositions>2 & mutationPositions<(seqLength-1)]
            countMutations = 0
            for(mutationPosition in mutationPositions){
                fivemerIndex = mutationPosition-2
                fivemer = fivemerSeq[fivemerIndex]
                GLfivemer = fivemerGL[fivemerIndex]
                fivemerCodonPos = fivemerCodon(fivemerIndex)
                fivemerReadingFrameCodon = substr(fivemer,fivemerCodonPos[1],fivemerCodonPos[3])
                fivemerReadingFrameCodonGL = substr(GLfivemer,fivemerCodonPos[1],fivemerCodonPos[3])
                if(!any(grep("N",fivemer)) & !any(grep("N",GLfivemer))){
                    if(model==0){
                        countMutations = countMutations + 1
                        MutationMatrix[GLfivemer] <- (MutationMatrix[GLfivemer] + 1)
                        MutationCountMatrix[GLfivemer] <- (MutationCountMatrix[GLfivemer] + 1)
                    }else{
                        if( translateCodonToAminoAcid(fivemerReadingFrameCodonGL)!="*" ){
                            countMutations = countMutations + 1
                            positionWithinCodon = which(fivemerCodonPos==3)
                            glNuc =  substr(fivemerReadingFrameCodonGL,positionWithinCodon,positionWithinCodon)
                            inputNuc =  substr(fivemerReadingFrameCodon,positionWithinCodon,positionWithinCodon)
                            MutationMatrix[GLfivemer] <- (MutationMatrix[GLfivemer] + HS5FModel@substitution[glNuc,inputNuc])
                            MutationCountMatrix[GLfivemer] <- (MutationCountMatrix[GLfivemer] + 1)
                        }
                    }
                }
            }
            
            seqMutability = MutationMatrix/BackgroundMatrix
            seqMutability = seqMutability/sum(seqMutability,na.rm=TRUE)
            #cat(inputMatrixIndex,"\t",countMutations,"\n")
            return(list("seqMutability"  = seqMutability,"numbMutations" = countMutations,"seqMutabilityCount" = MutationCountMatrix, "BackgroundMatrix"=BackgroundMatrix))
            
        }
    }
    
}

#Returns the codon position containing the middle nucleotide
fivemerCodon <- function(fivemerIndex){
    codonPos = list(2:4,1:3,3:5)
    fivemerType = fivemerIndex%%3
    return(codonPos[[fivemerType+1]])
}

#returns probability values for one mutation in codons resulting in R, S or Stop
probMutations <- function(typeOfMutation){
    matMutationProb <- matrix(0,ncol=3,nrow=125,dimnames=list(words(alphabet = c(NUCLEOTIDES,"N"), length=3),c(1:3)))
    for(codon in rownames(matMutationProb)){
        if( !any(grep("N",codon)) ){
            for(muPos in 1:3){
                matCodon = matrix(rep(s2c(codon),3),nrow=3,ncol=3,byrow=T)
                glNuc = matCodon[1,muPos]
                matCodon[,muPos] = canMutateTo(glNuc)
                substitutionRate = HS5FModel@substitution[glNuc,matCodon[,muPos]]
                typeOfMutations = apply(rbind(rep(codon,3),apply(matCodon,1,c2s)),2,function(x){mutationType(c2s(x[1]),c2s(x[2]))})
                matMutationProb[codon,muPos] <- sum(substitutionRate[typeOfMutations==typeOfMutation])
            }
        }
    }
    
    return(matMutationProb)
}




#Mapping Trinucleotides to fivemers
mapTriToFivemer <- function(triMutability=triMutability_Literature_Human){
    rownames(triMutability) <- triMutability_Names
    Fivemer<-rep(NA,1024)
    names(Fivemer)<-words(alphabet=NUCLEOTIDES,length=5)
    Fivemer<-sapply(names(Fivemer),function(Word)return(sum( c(triMutability[substring(Word,3,5),1],triMutability[substring(Word,2,4),2],triMutability[substring(Word,1,3),3]),na.rm=TRUE)))
    Fivemer<-Fivemer/sum(Fivemer)
    return(Fivemer)
}

collapseFivemerToTri<-function(Fivemer,Weights=MutabilityWeights,position=1,NUC="A"){
    Indices<-substring(names(Fivemer),3,3)==NUC
    Factors<-substring(names(Fivemer[Indices]),(4-position),(6-position))
    tapply(which(Indices),Factors,function(i)weighted.mean(Fivemer[i],Weights[i],na.rm=TRUE))
}



CountFivemerToTri<-function(Fivemer,Weights=MutabilityWeights,position=1,NUC="A"){
    Indices<-substring(names(Fivemer),3,3)==NUC
    Factors<-substring(names(Fivemer[Indices]),(4-position),(6-position))
    tapply(which(Indices),Factors,function(i)sum(Weights[i],na.rm=TRUE))
}

#Uses the real counts of the mutated fivemers
CountFivemerToTri2<-function(Fivemer,Counts=MutabilityCounts,position=1,NUC="A"){
    Indices<-substring(names(Fivemer),3,3)==NUC
    Factors<-substring(names(Fivemer[Indices]),(4-position),(6-position))
    tapply(which(Indices),Factors,function(i)sum(Counts[i],na.rm=TRUE))
}

bootstrap<-function(x=c(33,12,21),M=10000,alpha=0.05){
    N<-sum(x)
    if(N){
        p<-x/N
        k<-length(x)-1
        tmp<-rmultinom(M, size = N, prob=p)
        tmp_p<-apply(tmp,2,function(y)y/N)
        (apply(tmp_p,1,function(y)quantile(y,c(alpha/2/k,1-alpha/2/k))))
    }
    else return(matrix(0,2,length(x)))
}




bootstrap2<-function(x=c(33,12,21),n=10,M=10000,alpha=0.05){
    
    N<-sum(x)
    k<-length(x)
    y<-rep(1:k,x)
    tmp<-sapply(1:M,function(i)sample(y,n))
    if(n>1)tmp_p<-sapply(1:M,function(j)sapply(1:k,function(i)sum(tmp[,j]==i)))/n
    if(n==1)tmp_p<-sapply(1:M,function(j)sapply(1:k,function(i)sum(tmp[j]==i)))/n
    (apply(tmp_p,1,function(z)quantile(z,c(alpha/2/(k-1),1-alpha/2/(k-1)))))
}



p_value<-function(x=c(33,12,21),M=100000,x_obs=c(2,5,3)){
    n=sum(x_obs)
    N<-sum(x)
    k<-length(x)
    y<-rep(1:k,x)
    tmp<-sapply(1:M,function(i)sample(y,n))
    if(n>1)tmp_p<-sapply(1:M,function(j)sapply(1:k,function(i)sum(tmp[,j]==i)))
    if(n==1)tmp_p<-sapply(1:M,function(j)sapply(1:k,function(i)sum(tmp[j]==i)))
    tmp<-rbind(sapply(1:3,function(i)sum(tmp_p[i,]>=x_obs[i])/M),
               sapply(1:3,function(i)sum(tmp_p[i,]<=x_obs[i])/M))
    sapply(1:3,function(i){if(tmp[1,i]>=tmp[2,i])return(-tmp[2,i])else return(tmp[1,i])})
}

#"D:\\Sequences\\IMGT Germlines\\Human_SNPless_IGHJ.FASTA"
# Remove SNPs from IMGT germline segment alleles
generateUnambiguousRepertoire <- function(repertoireInFile,repertoireOutFile){
    repertoireIn <- read.fasta(repertoireInFile, seqtype="DNA",as.string=T,set.attributes=F,forceDNAtolower=F)
    alleleNames <- sapply(names(repertoireIn),function(x)strsplit(x,"|",fixed=TRUE)[[1]][2])
    SNPs <- tapply(repertoireIn,sapply(alleleNames,function(x)strsplit(x,"*",fixed=TRUE)[[1]][1]),function(x){
        Indices<-NULL
        for(i in 1:length(x)){
            firstSeq = s2c(x[[1]])
            iSeq = s2c(x[[i]])
            Indices<-c(Indices,which(firstSeq[1:320]!=iSeq[1:320] & firstSeq[1:320]!="." & iSeq[1:320]!="."  ))
        }
        return(sort(unique(Indices)))
    })
    repertoireOut <- repertoireIn
    repertoireOut <- lapply(names(repertoireOut), function(repertoireName){
        alleleName <- strsplit(repertoireName,"|",fixed=TRUE)[[1]][2]
        geneSegmentName <- strsplit(alleleName,"*",fixed=TRUE)[[1]][1]
        alleleSeq <- s2c(repertoireOut[[repertoireName]])
        alleleSeq[as.numeric(unlist(SNPs[geneSegmentName]))] <- "N"
        alleleSeq <- c2s(alleleSeq)
        repertoireOut[[repertoireName]] <- alleleSeq
    })
    names(repertoireOut) <- names(repertoireIn)
    write.fasta(repertoireOut,names(repertoireOut),file.out=repertoireOutFile)
    
}






############
groupBayes2 = function(indexes, param_resultMat){
    
    BayesGDist_Focused_CDR = calculate_bayesG( x=param_resultMat[indexes,1], N=apply(param_resultMat[indexes,c(1,2,4)],1,sum,na.rm=T), p=apply(param_resultMat[indexes,5:8],1,function(x){x[1]/(x[1]+x[2]+x[4])}))
    BayesGDist_Focused_FWR = calculate_bayesG( x=param_resultMat[indexes,3], N=apply(param_resultMat[indexes,c(3,2,4)],1,sum,na.rm=T), p=apply(param_resultMat[indexes,5:8],1,function(x){x[3]/(x[3]+x[2]+x[4])}))
    #BayesGDist_Local_CDR = calculate_bayesG( x=param_resultMat[indexes,1], N=apply(param_resultMat[indexes,c(1,2)],1,sum,na.rm=T), p=apply(param_resultMat[indexes,5:8],1,function(x){x[1]/(x[1]+x[2])}))
    #BayesGDist_Local_FWR = calculate_bayesG( x=param_resultMat[indexes,3], N=apply(param_resultMat[indexes,c(3,4)],1,sum,na.rm=T), p=apply(param_resultMat[indexes,5:8],1,function(x){x[3]/(x[3]+x[4])}))
    #BayesGDist_Global_CDR = calculate_bayesG( x=param_resultMat[indexes,1], N=apply(param_resultMat[indexes,c(1,2,3,4)],1,sum,na.rm=T), p=apply(param_resultMat[indexes,5:8],1,function(x){x[1]/(x[1]+x[2]+x[3]+x[4])}))
    #BayesGDist_Global_FWR = calculate_bayesG( x=param_resultMat[indexes,3], N=apply(param_resultMat[indexes,c(1,2,3,4)],1,sum,na.rm=T), p=apply(param_resultMat[indexes,5:8],1,function(x){x[3]/(x[1]+x[2]+x[3]+x[4])}))
    return ( list("BayesGDist_Focused_CDR"=BayesGDist_Focused_CDR,
                  "BayesGDist_Focused_FWR"=BayesGDist_Focused_FWR) )
    #"BayesGDist_Local_CDR"=BayesGDist_Local_CDR,
    #"BayesGDist_Local_FWR" = BayesGDist_Local_FWR))
    #                "BayesGDist_Global_CDR" = BayesGDist_Global_CDR,
    #                "BayesGDist_Global_FWR" = BayesGDist_Global_FWR) )
    
    
}


calculate_bayesG <- function( x=array(), N=array(), p=array(), max_sigma=20, length_sigma=4001){
    G <- max(length(x),length(N),length(p))
    x=array(x,dim=G)
    N=array(N,dim=G)
    p=array(p,dim=G)
    
    indexOfZero = N>0 & p>0
    N = N[indexOfZero]
    x = x[indexOfZero]
    p = p[indexOfZero]
    G <- length(x)
    
    if(G){
        
        cons<-array( dim=c(length_sigma,G) )
        if(G==1) {
            return(calculate_bayes(x=x[G],N=N[G],p=p[G],max_sigma=max_sigma,length_sigma=length_sigma))
        }
        else {
            for(g in 1:G) cons[,g] <- calculate_bayes(x=x[g],N=N[g],p=p[g],max_sigma=max_sigma,length_sigma=length_sigma)
            listMatG <- convolutionPowersOfTwoByTwos(cons,length_sigma=length_sigma)
            y<-calculate_bayesGHelper(listMatG,length_sigma=length_sigma)
            return( y/sum(y)/(2*max_sigma/(length_sigma-1)) )
        }
    }else{
        return(NA)
    }
}


calculate_bayesGHelper <- function( listMatG,length_sigma=4001 ){
    matG <- listMatG[[1]]
    groups <- listMatG[[2]]
    i = 1
    resConv <- matG[,i]
    denom <- 2^groups[i]
    if(length(groups)>1){
        while( i<length(groups) ){
            i = i + 1
            resConv <- weighted_conv(resConv, matG[,i], w= {{2^groups[i]}/denom} ,length_sigma=length_sigma)
            #cat({{2^groups[i]}/denom},"\n")
            denom <- denom + 2^groups[i]
        }
    }
    return(resConv)
}

weighted_conv<-function(x,y,w=1,m=100,length_sigma=4001){
    lx<-length(x)
    ly<-length(y)
    if({lx<m}| {{lx*w}<m}| {{ly}<m}| {{ly*w}<m}){
        if(w<1){
            y1<-approx(1:ly,y,seq(1,ly,length.out=m))$y
            x1<-approx(1:lx,x,seq(1,lx,length.out=m/w))$y
            lx<-length(x1)
            ly<-length(y1)
        }
        else {
            y1<-approx(1:ly,y,seq(1,ly,length.out=m*w))$y
            x1<-approx(1:lx,x,seq(1,lx,length.out=m))$y
            lx<-length(x1)
            ly<-length(y1)
        }
    }
    else{
        x1<-x
        y1<-approx(1:ly,y,seq(1,ly,length.out=floor(lx*w)))$y
        ly<-length(y1)
    }
    tmp<-approx(x=1:(lx+ly-1),y=convolve(x1,rev(y1),type="open"),xout=seq(1,lx+ly-1,length.out=length_sigma))$y
    tmp[tmp<=0] = 0
    return(tmp/sum(tmp))
}

########################




mutabilityMatrixONE <- rep(0,4)
names(mutabilityMatrixONE) <- NUCLEOTIDES[1:4]

# triMutability Background Count
buildMutabilityModelONE <- function( inputMatrixIndex, model=0 , multipleMutation=0, seqWithStops=0, stopMutations=0){
    
    #rowOrigMatInput = matInput[inputMatrixIndex,]
    seqGL =  gsub("-", "", matInput[inputMatrixIndex,2])
    seqInput = gsub("-", "", matInput[inputMatrixIndex,1])
    matInput[inputMatrixIndex,] <<- c(seqInput,seqGL)
    seqLength = nchar(seqGL)
    mutationCount <- analyzeMutations(inputMatrixIndex, model, multipleMutation, seqWithStops)
    BackgroundMatrix = mutabilityMatrixONE
    MutationMatrix = mutabilityMatrixONE
    MutationCountMatrix = mutabilityMatrixONE
    if(!is.na(mutationCount)){
        if((stopMutations==0 & model==0) | (stopMutations==1 & (sum(mutationCount=="Stop")<length(mutationCount))) | (model==1 & (sum(mutationCount=="S")>0)) ){
            
            #         ONEmerStartPos = 1:(seqLength)
            #         ONEmerLength <- length(ONEmerStartPos)
            ONEmerGL <- s2c(seqGL)
            ONEmerSeq <- s2c(seqInput)
            
            #Background
            for(ONEmerIndex in 1:seqLength){
                ONEmer = ONEmerGL[ONEmerIndex]
                if(ONEmer!="N"){
                    ONEmerCodonPos = getCodonPos(ONEmerIndex)
                    ONEmerReadingFrameCodon = c2s(ONEmerGL[ONEmerCodonPos])
                    ONEmerReadingFrameCodonInputSeq = c2s(ONEmerSeq[ONEmerCodonPos] )
                    
                    # All mutations model
                    #if(!any(grep("N",ONEmerReadingFrameCodon))){
                    if(model==0){
                        if(stopMutations==0){
                            if(!any(grep("N",ONEmerReadingFrameCodonInputSeq)))
                                BackgroundMatrix[ONEmer] <- (BackgroundMatrix[ONEmer] + 1)
                        }else{
                            if( !any(grep("N",ONEmerReadingFrameCodonInputSeq)) & translateCodonToAminoAcid(ONEmerReadingFrameCodonInputSeq)!="*"){
                                positionWithinCodon = which(ONEmerCodonPos==ONEmerIndex)#positionsWithinCodon[(ONEmerCodonPos[1]%%3)+1]
                                BackgroundMatrix[ONEmer] <- (BackgroundMatrix[ONEmer] + probNonStopMutations[ONEmerReadingFrameCodon,positionWithinCodon])
                            }
                        }
                    }else{ # Only silent mutations
                        if( !any(grep("N",ONEmerReadingFrameCodonInputSeq)) & translateCodonToAminoAcid(ONEmerReadingFrameCodonInputSeq)!="*" & translateCodonToAminoAcid(ONEmerReadingFrameCodonInputSeq)==translateCodonToAminoAcid(ONEmerReadingFrameCodon) ){
                            positionWithinCodon = which(ONEmerCodonPos==ONEmerIndex)
                            BackgroundMatrix[ONEmer] <- (BackgroundMatrix[ONEmer] + probSMutations[ONEmerReadingFrameCodon,positionWithinCodon])
                        }
                    }
                }
            }
        }
        
        #Mutations
        if(stopMutations==1) mutationCount = mutationCount[mutationCount!="Stop"]
        if(model==1) mutationCount = mutationCount[mutationCount=="S"]
        mutationPositions = as.numeric(names(mutationCount))
        mutationCount = mutationCount[mutationPositions>2 & mutationPositions<(seqLength-1)]
        mutationPositions =  mutationPositions[mutationPositions>2 & mutationPositions<(seqLength-1)]
        countMutations = 0
        for(mutationPosition in mutationPositions){
            ONEmerIndex = mutationPosition
            ONEmer = ONEmerSeq[ONEmerIndex]
            GLONEmer = ONEmerGL[ONEmerIndex]
            ONEmerCodonPos = getCodonPos(ONEmerIndex)
            ONEmerReadingFrameCodon = c2s(ONEmerSeq[ONEmerCodonPos])
            ONEmerReadingFrameCodonGL =c2s(ONEmerGL[ONEmerCodonPos])
            if(!any(grep("N",ONEmer)) & !any(grep("N",GLONEmer))){
                if(model==0){
                    countMutations = countMutations + 1
                    MutationMatrix[GLONEmer] <- (MutationMatrix[GLONEmer] + 1)
                    MutationCountMatrix[GLONEmer] <- (MutationCountMatrix[GLONEmer] + 1)
                }else{
                    if( translateCodonToAminoAcid(ONEmerReadingFrameCodonGL)!="*" ){
                        countMutations = countMutations + 1
                        positionWithinCodon = which(ONEmerCodonPos==ONEmerIndex)
                        glNuc =  substr(ONEmerReadingFrameCodonGL,positionWithinCodon,positionWithinCodon)
                        inputNuc =  substr(ONEmerReadingFrameCodon,positionWithinCodon,positionWithinCodon)
                        MutationMatrix[GLONEmer] <- (MutationMatrix[GLONEmer] + HS5FModel@substitution[glNuc,inputNuc])
                        MutationCountMatrix[GLONEmer] <- (MutationCountMatrix[GLONEmer] + 1)
                    }
                }
            }
        }
        
        seqMutability = MutationMatrix/BackgroundMatrix
        seqMutability = seqMutability/sum(seqMutability,na.rm=TRUE)
        #cat(inputMatrixIndex,"\t",countMutations,"\n")
        return(list("seqMutability"  = seqMutability,"numbMutations" = countMutations,"seqMutabilityCount" = MutationCountMatrix, "BackgroundMatrix"=BackgroundMatrix))
        #         tmp<-list("seqMutability"  = seqMutability,"numbMutations" = countMutations,"seqMutabilityCount" = MutationCountMatrix)
    }
}

################
# $Id: trim.R 989 2006-10-29 15:28:26Z ggorjan $

trim <- function(s, recode.factor=TRUE, ...)
    UseMethod("trim", s)

trim.default <- function(s, recode.factor=TRUE, ...)
    s

trim.character <- function(s, recode.factor=TRUE, ...)
{
    s <- sub(pattern="^ +", replacement="", x=s)
    s <- sub(pattern=" +$", replacement="", x=s)
    s
}

trim.factor <- function(s, recode.factor=TRUE, ...)
{
    levels(s) <- trim(levels(s))
    if(recode.factor) {
        dots <- list(x=s, ...)
        if(is.null(dots$sort)) dots$sort <- sort
        s <- do.call(what=reorder.factor, args=dots)
    }
    s
}

trim.list <- function(s, recode.factor=TRUE, ...)
    lapply(s, trim, recode.factor=recode.factor, ...)

trim.data.frame <- function(s, recode.factor=TRUE, ...)
{
    s[] <- trim.list(s, recode.factor=recode.factor, ...)
    s
}
#######################################
# Compute the expected for each sequence-germline pair by codon
getExpectedIndividualByCodon <- function(matInput){
    if( any(grep("multicore",search())) ){
        facGL <- factor(matInput[,2])
        facLevels = levels(facGL)
        LisGLs_MutabilityU = mclapply(1:length(facLevels),  function(x){
            computeMutabilities(facLevels[x])
        })
        facIndex = match(facGL,facLevels)
        
        LisGLs_Mutability = mclapply(1:nrow(matInput),  function(x){
            cInput = rep(NA,nchar(matInput[x,1]))
            cInput[s2c(matInput[x,1])!="N"] = 1
            LisGLs_MutabilityU[[facIndex[x]]] * cInput
        })
        
        LisGLs_Targeting =  mclapply(1:dim(matInput)[1],  function(x){
            computeTargeting(matInput[x,2],LisGLs_Mutability[[x]])
        })
        
        LisGLs_MutationTypes  = mclapply(1:length(matInput[,2]),function(x){
            #print(x)
            computeMutationTypes(matInput[x,2])
        })
        
        LisGLs_R_Exp = mclapply(1:nrow(matInput),  function(x){
            Exp_R <-  rollapply(as.zoo(1:VLENGTH),width=3,by=3,
                                function(codonNucs){
                                    RPos = which(LisGLs_MutationTypes[[x]][,codonNucs]=="R")
                                    sum( LisGLs_Targeting[[x]][,codonNucs][RPos], na.rm=T )
                                }
            )
        })
        
        LisGLs_S_Exp = mclapply(1:nrow(matInput),  function(x){
            Exp_S <-  rollapply(as.zoo(1:VLENGTH),width=3,by=3,
                                function(codonNucs){
                                    SPos = which(LisGLs_MutationTypes[[x]][,codonNucs]=="S")
                                    sum( LisGLs_Targeting[[x]][,codonNucs][SPos], na.rm=T )
                                }
            )
        })
        
        Exp_R = matrix(unlist(LisGLs_R_Exp),nrow=nrow(matInput),ncol=VLENGTH/3,T)
        Exp_S = matrix(unlist(LisGLs_S_Exp),nrow=nrow(matInput),ncol=VLENGTH/3,T)
        return( list( "Expected_R"=Exp_R, "Expected_S"=Exp_S) )
    }else{
        facGL <- factor(matInput[,2])
        facLevels = levels(facGL)
        LisGLs_MutabilityU = lapply(1:length(facLevels),  function(x){
            computeMutabilities(facLevels[x])
        })
        facIndex = match(facGL,facLevels)
        
        LisGLs_Mutability = lapply(1:nrow(matInput),  function(x){
            cInput = rep(NA,nchar(matInput[x,1]))
            cInput[s2c(matInput[x,1])!="N"] = 1
            LisGLs_MutabilityU[[facIndex[x]]] * cInput
        })
        
        LisGLs_Targeting =  lapply(1:dim(matInput)[1],  function(x){
            computeTargeting(matInput[x,2],LisGLs_Mutability[[x]])
        })
        
        LisGLs_MutationTypes  = lapply(1:length(matInput[,2]),function(x){
            #print(x)
            computeMutationTypes(matInput[x,2])
        })
        
        LisGLs_R_Exp = lapply(1:nrow(matInput),  function(x){
            Exp_R <-  rollapply(as.zoo(1:VLENGTH),width=3,by=3,
                                function(codonNucs){
                                    RPos = which(LisGLs_MutationTypes[[x]][,codonNucs]=="R")
                                    sum( LisGLs_Targeting[[x]][,codonNucs][RPos], na.rm=T )
                                }
            )
        })
        
        LisGLs_S_Exp = lapply(1:nrow(matInput),  function(x){
            Exp_S <-  rollapply(as.zoo(1:VLENGTH),width=3,by=3,
                                function(codonNucs){
                                    SPos = which(LisGLs_MutationTypes[[x]][,codonNucs]=="S")
                                    sum( LisGLs_Targeting[[x]][,codonNucs][SPos], na.rm=T )
                                }
            )
        })
        
        Exp_R = matrix(unlist(LisGLs_R_Exp),nrow=nrow(matInput),ncol=VLENGTH/3,T)
        Exp_S = matrix(unlist(LisGLs_S_Exp),nrow=nrow(matInput),ncol=VLENGTH/3,T)
        return( list( "Expected_R"=Exp_R, "Expected_S"=Exp_S) )
    }
}

# getObservedMutationsByCodon <- function(listMutations){
#   numbSeqs <- length(listMutations)
#   obsMu_R <- matrix(0,nrow=numbSeqs,ncol=VLENGTH/3,dimnames=list(c(1:numbSeqs),c(1:(VLENGTH/3))))
#   obsMu_S <- obsMu_R
#   temp <- mclapply(1:length(listMutations), function(i){
#     arrMutations = listMutations[[i]]
#     RPos = as.numeric(names(arrMutations)[arrMutations=="R"])
#     RPos <- sapply(RPos,getCodonNumb)
#     if(any(RPos)){
#       tabR <- table(RPos)
#       obsMu_R[i,as.numeric(names(tabR))] <<- tabR
#     }
#
#     SPos = as.numeric(names(arrMutations)[arrMutations=="S"])
#     SPos <- sapply(SPos,getCodonNumb)
#     if(any(SPos)){
#       tabS <- table(SPos)
#       obsMu_S[i,names(tabS)] <<- tabS
#     }
#   }
#   )
#   return( list( "Observed_R"=obsMu_R, "Observed_S"=obsMu_S) )
# }

getObservedMutationsByCodon <- function(listMutations){
    numbSeqs <- length(listMutations)
    obsMu_R <- matrix(0,nrow=numbSeqs,ncol=VLENGTH/3,dimnames=list(c(1:numbSeqs),c(1:(VLENGTH/3))))
    obsMu_S <- obsMu_R
    temp <- lapply(1:length(listMutations), function(i){
        arrMutations = listMutations[[i]]
        RPos = as.numeric(names(arrMutations)[arrMutations=="R"])
        RPos <- sapply(RPos,getCodonNumb)
        if(any(RPos)){
            tabR <- table(RPos)
            obsMu_R[i,as.numeric(names(tabR))] <- tabR
        }
        
        SPos = as.numeric(names(arrMutations)[arrMutations=="S"])
        SPos <- sapply(SPos,getCodonNumb)
        if(any(SPos)){
            tabS <- table(SPos)
            obsMu_S[i,names(tabS)] <- tabS
        }
    }
    )
    return( list( "Observed_R"=obsMu_R, "Observed_S"=obsMu_S) )
}


