# Region definition extension for including FWR2-4 and CDR2-3

#' @include Shazam.R
NULL

#### Extend region definition to CDR3 and FWR4 ####

#' Generate a data frame from graph object and clone data frame 
#'
#' @details
#' \code{makeGraphDf} adds columns and rows to the clones database: 
#' 
#' Additional \strong{columns} are added for parent_sequence and parent 
#' (which is the parent sequence id).
#' 
#' Additional \strong{rows} are added for inferred sequences and the germline of the clone graph.
#' 
#' \link{makeGraphDf} also renames sequence_id content according to the following 
#' (assume clone number is 34):  
#' 
#' 34_Germline, 34_Inferred1, 34_1, 34_2, 34_3, 34_Inferred2, 34_4, etc.
#' 
#' The original sequence id is kept under a new column named \code{orig_sequence_id}, 
#' and the original parent sequence id is kept under a new column named \code{orig_parent}.
#' 
#' Note: sequence field name of \code{curCloneGraph} argument must be "sequence".
#' 
#' @param   curCloneGraph   an \link{igraph} object of the specific clone.
#' @param   curCloneObj     \code{ChangeoClone} object of the specific clone. 
#' @param   objSeqId        sequence id field name of \code{curCloneObj}
#' @param   objSeq           sequence field name of \code{curCloneObj}
#' 
#' @return A \code{ChangeoClone} object with additional columns (parent_sequence and parent)
#'         and additional rows (for germline and inferred sequences)
#'         
#' @examples 
#' \dontrun{
#' library("igraph")
#' library("dplyr")
#' # Load and subset example data:
#' data(ExampleDb, package="alakazam")
#' clone_3170_db <- subset(ExampleDb, clone_id == 3170)
#' clone_3170_obj <- makeChangeoClone(clone_3170_db, 
#'                       seq="sequence_alignment",
#'                        germ="germline_alignment")
#' dnapars_exec <- "~/apps/phylip-3.69/dnapars"
#' clone_3170_graph <- buildPhylipLineage(clone_3170_obj, 
#'                        dnapars_exec, rm_temp = TRUE)  
#' clone_3170_GraphDf <- makeGraphDf(clone_3170_graph, clone_3170_obj)
#' }
#' @export
makeGraphDf <- function(curCloneGraph, curCloneObj,objSeqId="sequence_id",objSeq="sequence") {
    # extracting the cur_clone_num from the inputs to function:
    cur_clone_num <- curCloneObj@clone
    # generating a data frame from the clone igraph object 
    # (- for getting the inferred sequences from the graph,
    # and the parent information):
    curCloneGraph_df <- summarizeSubtrees(curCloneGraph, fields="sequence")
    # merging the db from clone object and from graph:
    cur_clone_merged_df <- merge(x=curCloneObj@data, y=curCloneGraph_df, 
                                 by.x=objSeqId, by.y="name", all=T)
    # Renaming sequence_id column to orig_sequence_id, and renaming parent 
    # to orig_parent:
    cur_clone_merged_df$orig_sequence_id <- cur_clone_merged_df[,objSeqId]
    cur_clone_merged_df$orig_parent <- cur_clone_merged_df$parent
    
    # uniquifying some values in merged data frame, and filling some 
    # missing values:
    
    #1. Replacing inferred sequences names with a unique name 
    # (using the clone number). 
    # Doing so for both sequence_id and parent and graph vertices
    cur_clone_merged_df$parent <- gsub(pattern="Inferred", 
                                       x=cur_clone_merged_df$parent, 
                                       replacement=paste("Inferred_", 
                                                         cur_clone_num, "_", 
                                                         sep=""))
    cur_clone_merged_df[,objSeqId] <- gsub(pattern="Inferred", 
                                            x=cur_clone_merged_df[,objSeqId], 
                                            replacement=paste("Inferred_", 
                                                              cur_clone_num, 
                                                              "_", sep=""))
    V(curCloneGraph)$label <- gsub(pattern="Inferred", 
                                     x=V(curCloneGraph)$label, 
                                     replacement=paste("Inferred_", cur_clone_num, 
                                                       "_", sep=""))

    #2. Replacing Germline sequence name with a unique name 
    # (using the clone number). 
    # Doing so for both sequence_id and parent and graph vertices:
    cur_clone_merged_df$parent <- gsub(pattern="Germline", 
                                       x=cur_clone_merged_df$parent, 
                                       replacement=paste("Germline_", 
                                                         cur_clone_num, sep=""))
    cur_clone_merged_df$sequence_id <- gsub(pattern="Germline", 
                                            x=cur_clone_merged_df[,objSeqId], 
                                            replacement=paste("Germline_", 
                                                              cur_clone_num, 
                                                              sep=""))
    V(curCloneGraph)$label <- gsub(pattern="Germline", 
                                     x=V(curCloneGraph)$label, 
                                     replacement= paste("Germline_", 
                                                        cur_clone_num, sep=""))
    
    #3. Now need to fill in missing values for germline sequence and 
    #   inffered sequences:
    cur_clone_merged_df$clone <- cur_clone_num
    cur_clone_merged_df$v_call <- curCloneObj@v_gene
    cur_clone_merged_df$j_call <- curCloneObj@j_gene
    cur_clone_merged_df$junction_length <- curCloneObj@junc_len
    cur_clone_merged_df$germline_imgt <- curCloneObj@germline
    
    #4. setting a new sequence_id column with following format:  
    #<CLONE_NUM>_<SEQUENCE_SERIAL_NUM>  
    #Except for Germline and Inferred names which will remain as is.   
    cur_clone_merged_df[,objSeqId] <- paste(cur_clone_num, "_", 
                                             c(1:length(cur_clone_merged_df$orig_sequence_id)), 
                                             sep="")
    cur_clone_merged_df[,objSeqId] <- ifelse(grepl("Germline|Inferred", 
                                                        cur_clone_merged_df$orig_sequence_id), 
                                                     paste(cur_clone_num, "_", 
                                                           cur_clone_merged_df$orig_sequence_id, 
                                                           sep=""), 
                                                  cur_clone_merged_df[,objSeqId])
    
    #5. Doing the same for parent column: 
    # setting a new parent column with following format:  
    # <CLONE_NUM>_<SEQUENCE_SERIAL_NUM>  
    #Except for Germline and Inferred names which will remain as is.
    cur_clone_merged_df$parent <- cur_clone_merged_df[match(cur_clone_merged_df$orig_parent,
                                                            cur_clone_merged_df$orig_sequence_id), 
                                                      objSeqId]
    
    #6. There are 2 sequence columns: one from curCloneGraph_df (names "sequence") 
    #   and one from curCloneObj@data (called per argument objSeq).
    # so first checking if objSeq=="sequence", and taking care accordingly:
    # Removing the sequence column that came from curCloneObj@data as it does not
    # include sequences of germline and inferred.
    # Setting the "sequence" column to be named "sequence"
    if (objSeq=="sequence") {
        cur_clone_merged_df <- cur_clone_merged_df %>% select(-!!rlang::sym("sequence.x"))
        cur_clone_merged_df <- rename(cur_clone_merged_df, sequence=!!rlang::sym("sequence.y"))
    } else {
        cur_clone_merged_df1<-cur_clone_merged_df1[!(names(cur_clone_merged_df1) %in% c(objSeq))]
    }
    
    #7. Adding the parent sequence as a new column:
    cur_clone_merged_df$parent_sequence <- cur_clone_merged_df[match(cur_clone_merged_df$parent, 
                                                                     cur_clone_merged_df[,objSeqId]), 
                                                               objSeq]         
    # filling the parent sequece of the Germline sequence to be its own sequence 
    # (=GERMLINE_IMGT):
    #cur_clone_merged_df <- mutate(cur_clone_merged_df, 
    #                              parent_sequence=ifelse(is.na(parent_sequence), 
    #                                                     GERMLINE_IMGT, 
    #                                                     parent_sequence))
    
    # now checking if the germline sequence is equal to its (only) child sequence. 
    # For example if "250_7" sequence parent is the "250_Germline" sequence, 
    # then mergin them to one line called "250_7_Germline". 
    germ_seq_line <- filter(cur_clone_merged_df,!!rlang::sym("orig_sequence_id")=="Germline")
    germ_seq <- germ_seq_line[,objSeq]
    germ_son_seq_line <- filter(cur_clone_merged_df,!!rlang::sym("orig_parent")=="Germline")
    germ_son_seq <- germ_son_seq_line[,objSeq]
    if (seqDist(germ_seq, germ_son_seq) == 0) {
        # removing from db the line of the germline:
        cur_clone_merged_df <- filter(cur_clone_merged_df, 
                                      !!rlang::sym("orig_sequence_id")!="Germline")
        # renaming the sequence id of the germline son - to include "Germline"
        # in its name:
        cur_clone_merged_df[,objSeqId]<-ifelse(cur_clone_merged_df[,"orig_parent"] == "Germline", 
                                                 paste(cur_clone_merged_df[,objSeqId], "_", 
                                                       "Germline", sep=""),
                                                 cur_clone_merged_df[,objSeqId])
                                                 
        # Change the parent SEQUENCE to be NA (as it is the Germline)
        cur_clone_merged_df <- mutate(cur_clone_merged_df, 
                                      parent=ifelse(!!rlang::sym("orig_parent") == "Germline", "NA", 
                                                    !!rlang::sym("parent")))
    }
    return(cur_clone_merged_df)
}

#' Defining a Region that will include also CDR3 and FWR4 based on junction length and sequence
#' 
#' @details  
#' This function gets as input a junction length and an imgt aligned sequence
#' and outputs a \link{RegionDefinition} object that includes following regions:   
#' 
#' \strong{For \code{regionDefinition="IMGT_VDJ_BY_REGIONS"}:}
#' 
#'- \strong{fwr1}: Bases 1 to 78 (based on \link{IMGT_V_BY_REGIONS} definitions)  
#'
#'- \strong{cdr1}: Bases 79 to 114 (based on \link{IMGT_V_BY_REGIONS} definitions) 
#'
#'- \strong{fwr2}: Bases 115 to 165 (based on \link{IMGT_V_BY_REGIONS} definitions) 
#' 
#'- \strong{cdr2}: Bases 166 to 195 (based on \link{IMGT_V_BY_REGIONS} definitions) 
#' 
#'- \strong{fwr3}: Bases 196 to 312 (based on \link{IMGT_V_BY_REGIONS} definitions)
#'
#'- \strong{cdr3}: Bases 313 to (313 + \code{juncLength} - 6) - since junction sequnece 
#'                 includes (on the left) the last codon from fwr3, and (on the right)  
#'                 the first codon from fwr4.  
#'        
#'- \strong{fwr4}: Bases (313 + \code{juncLength} - 6 + 1) to sequence_length. 
#' 
#' \strong{For \code{regionDefinition}="IMGT_VDJ":}
#' 
#'- \strong{fwr}: Bases	1 to 78 (based on \link{IMGT_V_BY_REGIONS} definitions)
#'
#'                Bases	115 to 165 (based on \link{IMGT_V_BY_REGIONS} definitions)
#'  
#'                Bases	196 to 312 (based on \link{IMGT_V_BY_REGIONS} definitions)
#'  
#'                Bases	(313 + \code{juncLength} - 6 + 1) to sequence_length.
#'            
#'- \strong{cdr}: Bases	79 to 114 (based on \link{IMGT_V_BY_REGIONS} definitions)
#'  
#'                Bases	166 to 195 (based on \link{IMGT_V_BY_REGIONS} definitions) 
#' 
#'                Bases	313 to (313 + \code{juncLength} - 6) - since junction sequnece 
#'                includes (on the left) the last codon from fwr3, and (on the right) 
#'                the first codon from fwr4.  
#'
#' Note: In case the \code{regionDefinition} argument is not one of the extended
#'       regions (\code{IMGT_VDJ_BY_REGIONS} or \code{IMGT_VDJ}) - then this
#'       function will return the \code{regionDefinition} as is.
#'
#' @param  juncLength         The junction length of the sequence
#' @param  sequenceImgt       The imgt aligned sequence
#' @param  regionDefinition   The \link{RegionDefinition} type to calculate
#'                            the regionDefinition for. Can be one of 2: 
#'                            \code{"IMGT_VDJ_BY_REGIONS"} or \code{"IMGT_VDJ"}. 
#'                            Only these 2 regions include all
#'                            CDR1/2/3 and FWR1/2/3/4 regions.
#'                            
#' @return a \link{RegionDefinition} object that includes CDR1/2/3 and 
#'         FWR1/2/3/4 for the specific \code{sequenceImgt}, 
#'         \code{juncLength} and \code{regionDefinition}.
#'         
#' @examples 
#' # Load and subset example data
#' data(ExampleDb, package="alakazam")  
#' juncLength <-ExampleDb[['junction_length']][1]
#' sequenceImgt<-ExampleDb[['sequence_alignment']][1]
#' seq_1_reg_def<-makeRegion(juncLength = juncLength, 
#'                           sequenceImgt = sequenceImgt, 
#'                           regionDefinition = IMGT_VDJ_BY_REGIONS)
#' @export
makeRegion <- function(juncLength, sequenceImgt,
                       regionDefinition=NULL) {

    if (!is.null(regionDefinition)) {
        
        if (!is(regionDefinition, "RegionDefinition")) {
            stop(deparse(substitute(regionDefinition)), " is not a valid RegionDefinition object")
        }
        
        if (regionDefinition@name %in% c("IMGT_VDJ_BY_REGIONS","IMGT_VDJ")) { 
            # all slots except for boundaries and seqLength are already defined in regionDefinition
            # First need to extract sequence length from sequence:
            seqLength <- nchar(sequenceImgt)
            
            # juncLength doesn't include alignment gaps, which are in sequenceImgt
            # and need to be added to correctly identify the boundaries
            # Also, sequence_alignment can have `.` that represent indels
            junction_length_helper <- !strsplit(sequenceImgt[[1]], "")[[1]] %in% c("-",".")
            junction_length_helper[1:310-1] <- 0
            if (juncLength>0) {
                junction_end <- which(cumsum(junction_length_helper[1:length(junction_length_helper)])==juncLength[[1]])[1]
                if (is.na(junction_end)) {
                    warning("junction ends past `sequenceImgt`")
                    junction_end <- nchar(sequenceImgt)
                    num_gaps <- sum(!junction_length_helper[310:junction_end])
                    cdr3_end <- junction_end
                } else {
                    num_gaps <- sum(!junction_length_helper[310:junction_end])
                    juncLength <- juncLength + num_gaps
                    cdr3_end <- 313 + as.integer(juncLength) - 6 - 1
                }

            } else {
                cdr3_end <- 0
            }
            # now for the boundaries slot:
            boundaries <- factor(shazam::IMGT_V_BY_REGIONS@boundaries, 
                                 levels=c(levels(shazam::IMGT_V_BY_REGIONS@boundaries), "cdr3", "fwr4"))
            if (cdr3_end > 312) {
                boundaries[313:cdr3_end] <- factor("cdr3")
                if (cdr3_end < nchar(sequenceImgt)) {
                    boundaries[(cdr3_end+1):seqLength] <- factor("fwr4")      
                } 
            } else {
                # I you are here, the junction is too short, <= 6nt
                warning("CDR3 end < CDR3 start. Couldn't identify CDR3 and FWR4. Aligned junction length is: ", juncLength) 
            }

            if (regionDefinition@name == "IMGT_VDJ") {
                boundaries <- gsub(pattern="fwr.", replacement = "fwr", x=boundaries, perl=TRUE)
                boundaries <- gsub(pattern="cdr.", replacement = "cdr", x=boundaries, perl=TRUE)
                boundaries <- factor(boundaries, levels=c("fwr", "cdr"))
            } 
            new("RegionDefinition", 
                name=regionDefinition@name, 
                description=regionDefinition@description, 
                boundaries=boundaries, seqLength=unname(seqLength), 
                regions=regionDefinition@regions, 
                labels=regionDefinition@labels, 
                citation=regionDefinition@citation)
        } else {
            regionDefinition  
        }
    } else {
        makeNullRegionDefinition(nchar(sequenceImgt))
    }
}


# Calculating an extended (=that includes cdr1/2/3 and fwr1/2/3/4) region definition 
# for a specific clone in database.
# Inputs: 
# - clone_num:        the clone number for which to calculate the region definition.
# - db:               a ChangeoClone database that includes clone numbers. 
# - seq_col:          the name of the db column containing the sequence that is imgt aligned.
# - juncLengthColumn:       the name of the db column containing the junction length.
# - clone_col:        the name of the db column containing the clone number.
# - regionDefinition: the region definition type to be output for this clone. 
# Output:
# A regionDefinition object for the specific clone 
# Note: regionDefinition needs to be calculated specifically for the clone if it
#       is of type IMGT_VDJ or IMGT_VDJ_BY_REGIONS, as it includes also cdr3 and fwr4
#       which are specific to clone.
# Note: The region definition is same for all sequences in clone - so doing it 
#       based on first sequence in clone.  
getCloneRegion <- function(clone_num, db, seq_col="sequence", 
                           juncLengthColumn="junction_length", 
                           clone_col="clone", 
                           regionDefinition=NULL) {
    # Check region definition
    if (!is.null(regionDefinition) & !is(regionDefinition, "RegionDefinition")) {
        stop(deparse(substitute(regionDefinition)), " is not a valid RegionDefinition object")
    }
    # subseting the db to lines for specific clone
    clone_db <- db[db[[clone_col]] == clone_num,]
    if ( length(unique(clone_db[[juncLengthColumn]])) >1 ) {
        stop("Expecting clones where all sequences have the same junction lenth. Different lengths found for clone ", clone_num)
    }
    # getting one of the sequences of the specific clone: 
    seq <- clone_db[[seq_col]][1]
    junc_len <- clone_db[[juncLengthColumn]][1]
    reg <- makeRegion(juncLength=junc_len, sequenceImgt=seq, 
                      regionDefinition=regionDefinition)
    return(reg)
}


# Status: experimental, not exported function
# data(oneseq_db, package="alakazam")
# germline_db <- list(
# "IGHV3-11*05"="CAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTCAAGCCTGGAGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTC............AGTGACTACTACATGAGCTGGATCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTTCATACATTAGTAGTAGT......AGTAGTTACACAAACTACGCAGACTCTGTGAAG...GGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTGTATTACTGTGCGAGAGA",
# "IGHD3-10*01"="GTATTACTATGGTTCGGGGAGTTATTATAAC",
# "IGHJ5*02"="ACAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
# )
# pja <- shazam:::plotJunctionAlignment(oneseq_db, germline_db,regionDefinition=IMGT_VDJ_BY_REGIONS)
# pja$p
plotJunctionAlignment <- function(db_row, germline_db, 
                                  sequence_alignment="sequence_alignment", 
                                  v_call="v_call", d_call="d_call", j_call="j_call",
                                  v_germline_start="v_germline_start", v_germline_end="v_germline_end",
                                  d_germline_start="d_germline_start", d_germline_end="d_germline_end",
                                  j_germline_start="j_germline_start", j_germline_end="j_germline_end",
                                  np1_length="np1_length", # np2_length="np2_length",
                                  junction="junction", junction_length="junction_length",
                                  germline_alignment="germline_alignment",
                                  regionDefinition=NULL
                                  
) {
    
    
    # Check for valid columns
    check <- checkColumns(db_row, c(sequence_alignment, 
                                    v_call, j_call, # d_call,
                                    v_germline_start, v_germline_end,
                                    # d_germline_start, # d_germline_end,
                                    j_germline_start, j_germline_end,
                                    np1_length, # np2_length, 
                                    junction, junction_length)
    )
    if (check != TRUE) { stop(check) }
    
    if (!is.null(d_call)) {
        d_columns <- c(d_call, d_germline_start, d_germline_end)
        found <- d_columns %in% colnames(db_row)
        if (any(!found)) {
            stop( "Column(s) ", paste(d_columns[!found], sep="", collpase=",") ," not found.")
        }
    }
    
    if (!germline_alignment %in% colnames(db_row)) {
        warning("The column germline_alignment doesn't exist. Assigned NA value.")
        db_row[[germline_alignment]] <- NA
    }
    
    # Check region definition
    if (!is.null(regionDefinition)) {
        if (!is(regionDefinition, "RegionDefinition")) {
            stop(deparse(substitute(regionDefinition)), " is not a valid RegionDefinition object")
        }
    } 
    
    junction_seq <- db_row[[junction]]
    v_allele <- getAllele(db_row[[v_call]], first=T)
    d_allele <- getAllele(db_row[[d_call]], first=T)
    j_allele <- getAllele(db_row[[j_call]], first=T)
    
    
    #if (!is.na(d_allele)) {
    v_germline <- germline_db[[v_allele]]
    if (!is.na(d_allele)) {
        d_germline <- germline_db[[d_allele]]
    }
    j_germline <- germline_db[[j_allele]]
    
    seq_aln <- db_row[[sequence_alignment]]
    germ_aln <- db_row[[germline_alignment]]
    
    getPositionDf <- function(seq, label) {
        sequence_chars <- strsplit(seq,"")[[1]]
        sequence_chars
        position_df <- data.frame("nucleotide"=sequence_chars, 
                                  "pos"=1:nchar(seq),
                                  "label"=label,
                                  stringsAsFactors = F)
        position_df
    }   
    sequence_aln_df <- getPositionDf(seq_aln,sequence_alignment)
    sequence_aln_df[['aligned']] <- T
    
    if (!is.na(germ_aln)) {
        germ_aln_df <- getPositionDf(germ_aln,germline_alignment)
        germ_aln_df[['aligned']] <- T   
    } else {
        germ_aln_df <- data.frame()
    }
    
    v_germline_df <- getPositionDf(v_germline,"v_germline")
    v_germ_start <- as.numeric(db_row[[v_germline_start]])
    v_germ_end <- as.numeric(db_row[[v_germline_end]])
    v_germline_df[['aligned']] <- v_germline_df$pos >= v_germ_start & v_germline_df$pos <= v_germ_end
    v_germ_start_off <- 1 - v_germ_start
    v_germline_df$pos <- v_germline_df$pos + v_germ_start_off
    v_length <- v_germ_end - v_germ_start + 1
    
    np1_len <- db_row[[np1_length]]
    
    if (!is.na(d_allele)) {
        d_germline_df <- getPositionDf(d_germline,"d_germline")
        d_germ_start <- as.numeric(db_row[[d_germline_start]])
        d_germ_end <- as.numeric(db_row[[d_germline_end]])
        d_germline_df[['aligned']] <- d_germline_df$pos >= d_germ_start & d_germline_df$pos <= d_germ_end
        d_germ_start_off <- 1-d_germ_start
        d_germline_df$pos <- v_length + np1_len + d_germline_df$pos + d_germ_start_off
        # d_length <- d_germ_end - d_germ_start + 1
    } else {
        d_germline_df <- data.frame()
    }
    # np2_len <- db_row[[np2_length]]
    
    j_germline_df <- getPositionDf(j_germline,"j_germline")
    j_germ_start <- as.numeric(db_row[[j_germline_start]])
    j_germ_end <- as.numeric(db_row[[j_germline_end]])
    j_germline_df[['aligned']] <- j_germline_df$pos >= j_germ_start & j_germline_df$pos <= j_germ_end
    j_length <- j_germ_end - j_germ_start + 1
    # Use end as reference
    j_germline_df$pos <- nchar(seq_aln) - j_length + 1 + j_germline_df$pos - j_germ_start
    
    df <- bind_rows(sequence_aln_df, v_germline_df, d_germline_df, j_germline_df, germ_aln_df)
    
    if (is.na(junction_seq) | db_row[[junction_length]]<1 ) {
        message("No junction available for this sequence.")   
        junction_df <- NULL
    } else {
        # junction_position <- #stri_locate(seq_aln,fixed=junction_seq)
        junction_start <- 310
        j_len <- db_row[[junction_length]]
        # junction_end <- junction_start + j_len - 1
        
        # Get aligned junction, with gaps
        junction_alignment_helper <- !stringi::stri_split_boundaries(seq_aln, type="char")[[1]] %in% c("-",".")
        junction_alignment_helper[1:junction_start-1] <- 0
        junction_end <- which(cumsum(junction_alignment_helper[1:length(junction_alignment_helper)])>j_len)[1] - 1
        if (is.na(junction_end)) {
            warning("The junction ends past sequence_alignment. Using the last position as junction_end.")
            junction_end <-  length(junction_alignment_helper)
        }
        
        junction_alignment <- stringi::stri_sub(seq_aln,310,junction_end)
        
        
        # junction_end <- junction_start + nchar(db_row[[junction]]) -1
        junction_df <- data.frame("nucleotide"=strsplit(junction_alignment,"")[[1]], stringsAsFactors = F)
        junction_df[['pos']] <-  c(junction_start:junction_end)
        junction_df[['label']] <- "junction"
        junction_df[['aligned']] <- T
    }        
    
    missing_junction_nt <- db_row[[junction_length]] - sum(junction_df$aligned)
    if ( missing_junction_nt > 0 ) {
        junction_nt <- gsub("[\\.\\-]","",db_row[['junction']])
        missing_nt <- stringi::stri_sub(junction_nt,
                                        nchar(junction_nt)-missing_junction_nt+1,
                                        nchar(junction_nt))
        junction_df_missing <- getPositionDf(missing_nt,"junction")
        junction_df_missing[['aligned']] <- F
        junction_df_missing[['pos']] <- junction_df_missing[['pos']] +  max(junction_df[['pos']])
        
        junction_df <- bind_rows(junction_df, junction_df_missing)
    }
    
    # addRegionDefinition boundaries
    region_definition <- makeRegion(db_row[[junction_length]], db_row[[sequence_alignment]], regionDefinition )
    
    if (!is.null(regionDefinition)) {
        rdf <- data.frame(
            "nucleotide"=as.character(region_definition@boundaries),
            "pos"=1:length(region_definition@boundaries),
            "label"="region_definition",
            "aligned"=T,
            stringsAsFactors = F)
    } else {
        rdf <- data.frame()
    }
    
    df <- bind_rows(df, junction_df, rdf)
    ordered_labels <- c("region_definition",sequence_alignment, "junction", "v_germline", "d_germline","j_germline", germline_alignment)
    df[['label']] <- factor(df[['label']], levels=rev(ordered_labels), ordered = T)
    
    color_palette <- c(
        "A" = "lightskyblue2",
        "T" = "khaki2",
        "G" = "palegreen3",
        "C" = "coral2",
        "." = "white", 
        "N" = "grey80",
        "-" = "black",
        "cdr1" = "#FFDCD2",
        "cdr2" = "#FFDCD2",
        "cdr3" = "#FFB189",
        "fwr1" = "#8CCEDB",
        "fwr2" = "#8CCEDB",
        "fwr3" = "#8CCEDB",
        "fwr4" = "#8CCEDB"         
    )
    
    fig_theme <- function(font_size=7) {
        theme_bw() +
            theme(text=element_text(size=font_size),
                 axis.title=element_text(size=font_size),
                 axis.text=element_text(size=font_size),
                 axis.text.x=element_text(size=font_size),
                 axis.text.y=element_text(size=font_size),
                 #axis.ticks=element_blank(),
                 #axis.ticks=theme_segment(colour = "black"),
                 #panel.background=element_rect(fill = NA, colour = "black", size = 0.25), 
                 panel.border=element_blank(),
                 #panel.grid.major=element_line(colour = "grey", size = 0.05),
                 #panel.grid.minor=element_line(colour = "grey", size = 0.05),
                 panel.grid.major.y=element_blank(),
                 panel.grid.minor.y=element_blank(),
                 panel.grid.major.x=element_blank(),
                 panel.grid.minor.x=element_blank(),
                 panel.spacing=unit(0.25, "lines"),
                 plot.title=element_text(size=font_size, 
                                        face="bold",
                                        lineheight=0.8),
                 legend.position="top",
                 legend.text=element_text(size=font_size),
                 # legend.title = element_text(size=font_size), 
                 legend.title=element_blank(),
                 legend.spacing=unit(0.25, "lines"),
                 legend.box="horizontal",
                 legend.box.spacing=unit(0.25, "lines"),
                 legend.key.height=unit(1,"line"),
                 legend.key.width=unit(1,"line"),
                 strip.text = element_text(size = font_size, face="plain"),
                 strip.background = element_blank(),
                 plot.margin= unit(c(0, 0, 0, 0), "lines"))   
    }
    
    p <- ggplot(data=df %>% 
                    dplyr::filter(.data$label == "region_definition") %>%
                    dplyr::mutate(label=factor(.data$label, levels = ordered_labels, ordered = T)), 
                aes(x=.data$pos, y=.data$label, fill=.data$nucleotide, alpha=.data$aligned)) +
        geom_tile( height=0.3) + 
        geom_tile(data=df %>% 
                      dplyr::filter(.data$label != "region_definition") %>%
                      dplyr::mutate(label=factor(.data$label, levels = ordered_labels, ordered = T)), 
                  color="grey50") + 
        fig_theme() +  
        scale_fill_manual(values=color_palette) +
        scale_alpha_manual(values=c("TRUE"=1, "FALSE"=0.2), guide="none") +
        scale_y_discrete(breaks=ordered_labels, labels=ordered_labels, limits=rev(ordered_labels), expand=c(0, 0)) +
        ylab("") + xlab("IMGT position") +
        guides(fill=guide_legend(nrow=1)) +
        scale_x_continuous(expand=c(0, 0)) 
    
   
    list(p=p, data=df)
}
