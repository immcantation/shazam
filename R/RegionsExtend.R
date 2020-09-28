# Region definition extention for including fwr/2/3/4 and cdr/2/3

#' @include Shazam.R
NULL


#### functions foe extending region definition also to cdr3 and fwr4 ####

# Generates a list of clones existing in db.
# Inputs:
# -db: a ChangeoClone database that includes a clone columns.
# - clone_col: the name of the clone columns in db.
# output: 
# A vector of all clone numbers in db (each clone number 
#  will appear only once in the list).  
# the transpose is needed to ge tthe proper length of the coerced value
makeClonesList <-  function(db, clone_col="clone_id") {
    clones_list <- as.list(t(unique(db[, c(clone_col)])))
    return(clones_list)
}


# Creating a ChangeoClone object for one specific clone:  
# Inputs: 
# - db: a ChangeoClone database
# - cur_clone_num: clone number for which to generate the ChangeoClone object.
# - id: name of the column containing sequence identifiers.
# - seq: name of the column containing observed DNA sequences. 
#        All sequences in this column must be multiple aligned.
# - germ: name of the column containing germline DNA sequences. 
#         All entries in this column should be identical for any given clone, 
#         and they must be multiple aligned with the data in the seq column.
# - v_call: name of the column containing V-segment allele assignments. 
#          All entries in this column should be identical to the gene level.
# - j_call: name of the column containing J-segment allele assignments. 
#          All entries in this column should be identical to the gene level.
# - junc_len: name of the column containing the length of the junction as a numeric value. 
#             All entries in this column should be identical for any given clone.
# - clone: name of the column containing the identifier for the clone. 
#          All entries in this column should be identical.
# - mask_char: character to use for masking and padding.
# - max_mask: maximum number of characters to mask at the leading and trailing sequence ends.
#             If NULL then the upper masking bound will be automatically determined from the
#             maximum number of observed leading or trailing Ns amongst all sequences. 
#             If set to 0 (default) then masking will not be performed
# - pad_end: if TRUE pad the end of each sequence with mask_char to make every sequence the same length.
# Output:
# A ChangeoClone object of the specific clone_num
#   - ChangeoClone object of the specific clone_num 
#   - Clone size (after collapsing relevant lines)
makeChangeoCloneCurClone <- function(db, cur_clone_num, id="sequence_id",  
                                    seq="SEQUENCE_IMGT", germ="GERMLINE_IMGT",  
                                    v_call="v_call", j_call="j_call",  
                                    junc_len="JUNCTION_LENGTH", clone="clone_id", 
                                    mask_char="N", max_mask=0,pad_end=FALSE,text_fields=NULL, 
                                    num_fields=NULL, seq_fields=NULL,
                                    add_count=TRUE, verbose=FALSE) {
    #clone_header <- as.character(as.name(clone))
    #cur_clone_db <- subset(db, clone_header == cur_clone_num)
    # subseting the db to lines for specific clone
    cur_clone_db <- db[db[[clone]] == cur_clone_num,]
    curCloneObj <- makeChangeoClone(data=cur_clone_db, id=id, seq=seq, germ=germ, 
                                      v_call=v_call, j_call=j_call, junc_len=junc_len, 
                                      clone=clone, mask_char=mask_char, 
                                      max_mask=max_mask, pad_end=pad_end, 
                                      text_fields=text_fields, 
                                      num_fields=num_fields, seq_fields=seq_fields,
                                      add_count=add_count, verbose=verbose)
    curCloneObj@data$v_call <- curCloneObj@v_gene
    curCloneObj@data$j_call <- curCloneObj@j_gene
    curCloneObj@data$clone_id <- curCloneObj@clone
    curCloneObj@data$germline_alignment <- curCloneObj@germline
    curCloneObj@data$junction_length <- curCloneObj@junc_len
    # marking the curCloneObj data size:
    # If only one line, it means that even though the original clone size was of more than 1 line 
    # - then all the lines were the same, and they all collapsed to one line while building a clone object.
    # special handling will be needed here.
    #cur_clone_size <- dim(curCloneObj@data)[1]
    return(curCloneObj)
}

# 
# Making an igraph object out of ChangeoClone object:
# Inputs:
# - curCloneObj: a chaneOclone object of one clone
# - dnapars_exec: absolute path to the PHYLIP dnapars executable
# Output:
# An Igrpah object of the specific clone ChangeoClone object
makeGraphCurClone <- function(curCloneObj,dnapars_exec,seq_id_col="sequence_id") {
    # The below line is in order to avoid warning message
    # (message: "binding factor and character vector, coercing into character vector")
    curCloneObj@data[[seq_id_col]] <- as.character(curCloneObj@data[[seq_id_col]]) 
    curCloneGraph <- buildPhylipLineage(curCloneObj,dnapars_exec,rm_temp = TRUE)  
    # in case that the cur_clone_num does not include at least 2 unique sequences - the curCloneGraph will be "NULL".
    #note: "N" instead of A/C/G/T does not make a sequence unique.
    return(curCloneGraph)
}

#' Generating a data frame from graph object, and mering it with clone data frame 
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
#' @examples 
#' \dontrun{
#' library("igraph")
#' library("dplyr")
#' # Load and subset example data:
#' data(ExampleDb)
#' clone_3170_db <- subset(ExampleDb, clone_id == 3170)
#' clone_3170_obj <- makeChangeoClone(clone_3170_db, seq="sequence_alignment", germ="germline_alignment")
#' dnapars_exec <- "~/apps/phylip-3.69/dnapars"
#' clone_3170_graph <- buildPhylipLineage(clone_3170_obj, dnapars_exec, rm_temp = TRUE)  
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
        cur_clone_merged_df <- subset(cur_clone_merged_df, select=c(-sequence.x))
        cur_clone_merged_df <- rename(cur_clone_merged_df, sequence=sequence.y)
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
    germ_seq_line <- filter(cur_clone_merged_df,orig_sequence_id=="Germline")
    germ_seq <- germ_seq_line[,objSeq]
    germ_son_seq_line <- filter(cur_clone_merged_df,orig_parent=="Germline")
    germ_son_seq <- germ_son_seq_line[,objSeq]
    if (seqDist(germ_seq, germ_son_seq) == 0) {
        # removing from db the line of the germline:
        cur_clone_merged_df <- filter(cur_clone_merged_df, 
                                      orig_sequence_id!="Germline")
        # renaming the sequence id of the germline son - to include "Germline"
        # in its name:
        cur_clone_merged_df[,objSeqId]<-ifelse(cur_clone_merged_df[,"orig_parent"] == "Germline", 
                                                 paste(cur_clone_merged_df[,objSeqId], "_", 
                                                       "Germline", sep=""),
                                                 cur_clone_merged_df[,objSeqId])
                                                 
        # Change the parent SEQUENCE to be NA (as it is the Germline)
        cur_clone_merged_df <- mutate(cur_clone_merged_df, 
                                      parent=ifelse(orig_parent == "Germline", "NA", 
                                                    parent))
    }
    return(cur_clone_merged_df)
}
#' Plotting a tree-plot of a specific clone: 
#' 
#' @details 
#' This function will take as input any \code{ChangeoClone} db, and a specific
#' clone number. It will build a \code{ChangeoClone} object and an \link{igraph}
#' object for the specific clone, and will plot a tree plot for it.
#' 
#' Notes: 
#' 
#' 1. This function will give an error in case that all the sequences in the 
#' db of the specific clone are the same (as no lineage tree can be built on that).
#' 
#' 2. In case one of the sequences is equal to the germline sequence - then the 
#' graph may show that sequence ID as the graph root (and not the "Germline" sequence ID).
#'  
#' @param   db          a \code{ChangeoClone} database
#' @param   curClone   clone number for which to plot the tree plot.
#' @param   id          name of the column containing sequence identifiers.
#' @param   seq         name of the column containing observed DNA sequences. 
#'                       All sequences in this column must be multiple aligned.
#' @param   germ        name of the column containing germline DNA sequences. 
#'                       All entries in this column should be identical for any given clone, 
#'                       and they must be multiple aligned with the data in the seq column.
#' @param   v_call        name of the column containing V-segment allele assignments. 
#'                       All entries in this column should be identical to the gene level.
#' @param   j_call        name of the column containing J-segment allele assignments. 
#'                       All entries in this column should be identical to the gene level.
#' @param   junc_len     name of the column containing the length of the junction as a numeric value. 
#'                       All entries in this column should be identical for any given clone.
#' @param   clone        name of the column containing the identifier for the clone. 
#'                       All entries in this column should be identical.
#' @param   mask_char    character to use for masking and padding.
#' @param   max_mask     maximum number of characters to mask at the leading and trailing sequence ends.
#'                       If NULL then the upper masking bound will be automatically determined from the
#'                       maximum number of observed leading or trailing Ns amongst all sequences. 
#'                       If set to 0 (default) then masking will not be performed
#' @param   pad_end      if TRUE pad the end of each sequence with mask_char to make every sequence the same length.
#' @param   dnapars_exec absolute path to the PHYLIP dnapars executable
#' 
#' @return  a plot of the specific clone linegae tree.
#' @examples 
#' library ("igraph")
#' library("alakazam")
#' data(ExampleDb)
#' dnapars_exec <- "~/apps/phylip-3.69/dnapars"
#' plotCloneLineageTree(db=ExampleDb, curClone=3139, seq="sequence_alignment",
#'                      id="sequence_id", germ="germline_alignment", dnapars_exec = dnapars_exec)
#' @export   
plotCloneLineageTree <- function(db, curClone=NULL, id="sequence_id", 
                                 seq="SEQUENCE_IMGT", germ="GERMLINE_IMGT", 
                                 v_call="v_call", j_call="j_call", 
                                 junc_len="JUNCTION_LENGTH", clone="clone_id", 
                                 mask_char="N", max_mask=0, pad_end=FALSE, 
                                 dnapars_exec) {
    curCloneObj <- makeChangeoCloneCurClone(db=db, cur_clone_num=curClone, id=id, 
                                         seq=seq,germ=germ, v_call=v_call, 
                                         j_call=j_call, junc_len=junc_len, 
                                         clone=clone, mask_char=mask_char, 
                                         max_mask=max_mask, pad_end=pad_end)
    curCloneGraph <- makeGraphCurClone(curCloneObj, dnapars_exec)
    V(curCloneGraph)$color <- ifelse(grepl("Germline", V(curCloneGraph)$label), 
                                       "lightgreen", 
                                       ifelse(grepl("Inferred", 
                                                    V(curCloneGraph)$label), 
                                              "lightpink", "lightblue"))
    #making the graph vertices show the simple sequence_id
    par(mar=c(0, 0, 0, 0) + 1)
    plot(curCloneGraph, layout=layout_as_tree, edge.arrow.mode=0, 
         vertex.frame.color="black", vertex.label.cex=0.8,
         vertex.label.color="black", vertex.size=18, 
         vertex.label=V(curCloneGraph)$seq_num, 
         main=paste("Lineage tree for clone ", curClone))
    legend("topright", c("Germline", "Inferred", "Sample"), 
           fill=c("lightgreen", "lightpink", "lightblue"), cex=0.75)
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
#' @return a \link{RegionDefinition} object that includes CDR1/2/3 and 
#'         FWR1/2/3/4 for the specific \code{sequenceImgt}, 
#'         \code{juncLength} and \code{regionDefinition}.
#' @examples 
#' Load and subset example data
#' data(ExampleDb)  
#' juncLength<-ExampleDb[1,"junction_length"]
#' sequenceImgt<-ExampleDb[1,"sequence_alignment"]
#' seq_1_reg_def<-makeRegion(juncLength = juncLength, 
#'                           sequenceImgt = sequenceImgt, 
#'                           regionDefinition = IMGT_VDJ_BY_REGIONS)
#' @export

makeRegion <- function(juncLength, sequenceImgt,  
                       regionDefinition=IMGT_VDJ_BY_REGIONS) {
    if (!is(regionDefinition, "RegionDefinition")) {
        stop(deparse(substitute(regionDefinition)), " is not a valid RegionDefinition object")
    }
    # all slots except for boundaries and seqLength are already defined in regionDefinition
    # First need to extract sequence length from sequence:
    seqLength <- nchar(sequenceImgt)
    # now for the boundaries slot:
    boundaries <- factor(IMGT_V_BY_REGIONS@boundaries, 
                         levels=c(levels(IMGT_V_BY_REGIONS@boundaries), "cdr3", "fwr4"))
    boundaries[313:(313 + as.integer(juncLength) - 6 - 1)] <- factor("cdr3")
    boundaries[(313 + as.integer(juncLength) - 6):seqLength] <- factor("fwr4")
    if (regionDefinition@name == "IMGT_VDJ") {
        boundaries <- gsub(pattern="fwr.", replacement = "fwr", x=boundaries, perl=TRUE)
        boundaries <- gsub(pattern="cdr.", replacement = "cdr", x=boundaries, perl=TRUE)
        boundaries <- factor(boundaries, levels=c("fwr", "cdr"))
    } 
    region_out <- new("RegionDefinition", 
                      name=regionDefinition@name, 
                      description=regionDefinition@description, 
                      boundaries=boundaries, seqLength=unname(seqLength), 
                      regions=regionDefinition@regions, 
                      labels=regionDefinition@labels, 
                      citation=regionDefinition@citation)
    # taking care of non-extended region definitions:
    if ((regionDefinition@name != "IMGT_VDJ_BY_REGIONS") & 
        (regionDefinition@name != "IMGT_VDJ")) {
        region_out <- regionDefinition
    }
    return(region_out)
}


# Calculating an extended (=that includes cdr1/2/3 and fwr1/2/3/4) region definition 
# for a specific clone in database.
# Inputs: 
# - clone_num:        the clone number for which to calculate the region definition.
# - db:               a ChangeoClone database that includes clone numbers. 
# - seq_col:          the name of the db column containing the sequence that is imgt aligned.
# - juncLenCol:       the name of the db column containing the junction length.
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
                           juncLenCol="junction_length", 
                           clone_col="clone", 
                           regionDefinition=IMGT_VDJ_BY_REGIONS) {
    # subseting the db to lines for specific clone
    clone_db <- db[db[[clone_col]] == clone_num,]
    # getting one of the sequences of the specific clone: 
    seq <- clone_db[1, seq_col]
    junc_len <- clone_db[1, juncLenCol]
    reg <- makeRegion(juncLength=junc_len, sequenceImgt=seq, 
                      regionDefinition=regionDefinition)
    return(reg)
}

# calculating consensus sequence for one clone only.  
# Note: works same as collapseClones function, but on one clone only, and thus
# can work with extended region definitions such as IMGT_VDJ and IMGT_VDJ_BY_REGIONS.
# Inputs:
# - clone_num:        the clone number for which to calculate the consensus sequence.
# - db:               data.frame containing sequence data. 
# - juncLenCol:       column name of junction length
# - cloneColumn:      name of the column containing clonal identifiers
# - sequenceColumn:   name of the column containing input sequences. 
#                     The length of each input sequence should match that of its 
#                     corresponding germline sequence.
# - regionDefinition: RegionDefinition object defining the regions and boundaries 
#                     of the Ig sequences. 
# - germlineColumn:   name of the column containing germline sequences.
#                     The length of each germline sequence should match that of 
#                     its corresponding input sequence.
# - muFreqColumn:     name of the column containing mutation frequency. 
#                     Applicable to the "mostMutated" and "leastMutated" methods. 
#                     If not supplied, mutation frequency is computed by calling 
#                     observedMutations. Default is NULL. See Cautions for note on usage.
# - method:           method for calculating input consensus sequence. 
#                     One of "thresholdedFreq", "mostCommon", "catchAll", "mostMutated", 
#                     or "leastMutated". See "Methods" for details.
# - minimumFrequency: frequency threshold for calculating input consensus sequence. 
#                     Applicable to and required for the "thresholdedFreq" method. 
#                     A canonical choice is 0.6. Default is NULL.
# - includeAmbiguous: whether to use ambiguous characters to represent positions at 
#                     which there are multiple characters with frequencies that are 
#                     at least minimumFrequency or that are maximal (i.e. ties). 
#                     Applicable to and required for the "thresholdedFreq" and "mostCommon" 
#                     methods. Default is FALSE. See "Choosing ambiguous characters" 
#                     for rules on choosing ambiguous characters.
# - breakTiesStochastic: In case of ties, whether to randomly pick a sequence from 
#                        sequences that fulfill the criteria as consensus. 
#                        Applicable to and required for all methods except for "catchAll". 
#                        Default is FALSE. See "Methods" for details.
# - breakTiesByColumns: A list of the form list(c(col_1, col_2, ...), c(fun_1, fun_2, ...)), 
#                       where col_i is a character name of a column in db, and fun_i 
#                       is a function to be applied on that column. Currently, 
#                       only max and min are supported. Note that the two c()'s in 
#                       list() are essential (i.e. if there is only 1 column, the 
#                       list should be of the form list(c(col_1), c(func_1)). 
#                       Applicable to and optional for the "mostMutated" and "leastMutated" 
#                       methods. If supplied, fun_i's are applied on col_i's to 
#                       help break ties. Default is NULL. See "Methods" for details.
# - expandedDb:         logical, indicating whether or not to return the expanded db, 
#                       containing all the sequences (as opposed to returning just 
#                       one sequence per clone).
# - nproc:              Number of cores to distribute the operation over. 
#                       If the cluster has already been set earlier, then pass 
#                       the cluster. This will ensure that it is not reset.
# Output:               
#                       consensus sequence for the specific clone and reion definition.
# 
collapseOneClone <- function(clone_num, db, juncLenCol="junction_length", 
                             cloneColumn = "clone_id", sequenceColumn = "sequence_alignment", 
                             regionDefinition = IMGT_VDJ_BY_REGIONS,
                             germlineColumn = "germline_alignment_d_mask", 
                             muFreqColumn = NULL, 
                             method=c("mostCommon","thresholdedFreq","catchAll","mostMutated","leastMutated"),
                             minimumFrequency = NULL,includeAmbiguous = FALSE, 
                             breakTiesStochastic = FALSE,
                             breakTiesByColumns = NULL, expandedDb = FALSE, nproc = 1) {
    clone_db <- db[db[[cloneColumn]] == clone_num,]
    clone_reg_def <- getCloneRegion(clone_num=clone_num, db=clone_db, 
                                    seq_col=sequenceColumn, 
                                    juncLenCol=juncLenCol, 
                                    clone_col=cloneColumn, 
                                    regionDefinition=regionDefinition)
    collapsed_clone <- collapseClonesL(db=clone_db, cloneColumn = cloneColumn, 
                                       sequenceColumn=sequenceColumn,
                                       germlineColumn=germlineColumn, 
                                       muFreqColumn = muFreqColumn, 
                                       regionDefinition = clone_reg_def, 
                                       method = method, 
                                       minimumFrequency = minimumFrequency, 
                                       includeAmbiguous = includeAmbiguous,
                                       breakTiesStochastic = breakTiesStochastic, 
                                       breakTiesByColumns=breakTiesByColumns, 
                                       expandedDb=expandedDb, nproc=nproc)
    return(collapsed_clone)
}


# This function calculates baseline for one clone in db.  
# It will first calculate the specific clone region definition 
# (which is specific to the clone).  
# Inputs:
# - clone_num:          the clone number for which to calculate the consensus sequence.
# - db:                 data.frame containing sequence data and annotations
# - sequenceColumn:     character name of the column in db containing input sequences.
# - cloneColumn:        the name of the db column containing the clone number.
# - juncLenCol:         column name of junction length 
# - germlineColumn:     character name of the column in db containing germline sequences.
# - testStatistic:      character indicating the statistical framework used to 
#                       test for selection. One of c("local", "focused", "imbalanced").
# - regionDefinition:   RegionDefinition object defining the regions and 
#                       boundaries of the Ig sequences.
# - targetingModel:     TargetingModel object. Default is HH_S5F.
# - mutationDefinition: MutationDefinition object defining replacement and silent 
#                       mutation criteria. If NULL then replacement and silent 
#                       are determined by exact amino acid identity. 
#                       Note, if the input data.frame already contains observed 
#                       and expected mutation frequency columns then mutations 
#                       will not be recalculated and this argument will be ignored.
# - calcStats:          logical indicating whether or not to calculate the summary 
#                       statistics data.frame stored in the stats slot of a 
#                       Baseline object.
# - nproc:              number of cores to distribute the operation over. 
#                       If nproc=0 then the cluster has already been set and will 
#                       not be reset.
# Output:
#                       BASELINe posterior probability density function (PDF) for 
#                       the specific clone and region definition.



calcBaselineOneClone <- function(clone_num, db, sequenceColumn = "sequence_alignment", 
                                 cloneColumn="clone_id",
                                 juncLenCol="junction_length",
                                 germlineColumn = "germline_alignment",
                                 testStatistic = c("local", "focused", "imbalanced"), 
                                 regionDefinition = IMGT_VDJ_BY_REGIONS,
                                 targetingModel = HH_S5F, 
                                 mutationDefinition = NULL,
                                 calcStats = FALSE, 
                                 nproc = 1) {
    clone_db <- db[db[[cloneColumn]] == clone_num,]
    reg_def <- getCloneRegion(clone_num=clone_num, db=clone_db, 
                              seq_col=sequenceColumn,
                              juncLenCol=juncLenCol,
                              clone_col=cloneColumn, regionDefinition=regionDefinition)
    clone_baseline <- calcBaselineL(db=clone_db,sequenceColumn = sequenceColumn,
                                    germlineColumn = germlineColumn, 
                                    testStatistic = testStatistic,
                                    regionDefinition = reg_def,
                                    targetingModel = targetingModel,
                                    mutationDefinition = mutationDefinition,
                                    calcStats = calcStats, nproc = nproc,
                                    cloneColumn=cloneColumn,
                                    juncLengthColumn=juncLenCol)
    return(clone_baseline)
}

