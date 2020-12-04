# Generates a list of clones existing in db.
# Inputs:
# -db: a ChangeoClone database that includes a clone columns.
# - clone_col: the name of the clone columns in db.
# output:
# A vector of all clone numbers in db (each clone number
#  will appear only once in the list).
# the transpose is needed to get the proper length of the coerced value
makeClonesList <-  function(db, clone_col="clone_id") {
    clones_list <- as.list(t(unique(db[, c(clone_col)])))
    return(clones_list)
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
                                 seq="sequence_alignment", germ="germline_alignment",
                                 v_call="v_call", j_call="j_call",
                                 junc_len="junction_length", clone="clone_id",
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
                                     seq="sequence_alignment", germ="germline_alignment",
                                     v_call="v_call", j_call="j_call",
                                     junc_len="junction_length", clone="clone_id",
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