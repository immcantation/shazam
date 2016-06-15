##### Functions for shmulate

### getPathLengths and getFounder from Topology.R ###

#' Calculate path lengths from a lineage tree root
#'
#' \code{getPathLengths} calculates the unweighted (number of steps) and weighted (distance) 
#' path lengths from the root of a lineage tree in an igraph object.
#'
#' @param    graph     igraph object with node annotations.
#' @param    root      name of the root (germline) node.
#' @param    field     annotation field to use for exclusion of nodes from step count.
#' @param    exclude   annotation values specifying which nodes to exclude from step count. 
#'                     if \code{NULL} consider all nodes. This does not affect the weighted
#'                     (distance) path length calculation.
#' @return   A data.frame with columns: name, steps, and distance.
#' 
#' @seealso  \link[igraph]{shortest.paths} and \link{getFounder}.
#' @examples
#' # Define simple graph
#' library(igraph)
#' graph <- graph(c(1,2,2,3,2,4,3,5,3,6), directed=TRUE)
#' V(graph)$name <- c("Germline", "Inferred", "Seq1", "Seq2", "Seq3", "Seq4")
#' V(graph)$isotype <- c(NA, NA, "IgM", "IgG", "IgA", "IgA")
#' E(graph)$weight <- c(10, 3, 6, 4, 1)
#' E(graph)$label <- E(graph)$weight
#' V(graph)$label <- V(graph)$name
#'
#' Plot graph with a tree layout
#' ly <- layout.reingold.tilford(graph, root="Germline", circular=F, flip.y=T)
#' plot(graph, layout=ly)
#' 
#' # Consider all nodes
#' getPathLengths(graph, root="Germline")
#' 
#' # Exclude nodes without an isotype annotation from step count
#' getPathLengths(graph, root="Germline", field="isotype", exclude=NA)
#' 
#' @export
getPathLengths <- function(graph, root="Germline", field=NULL, exclude=NULL) {
    # TODO:  should steps be level?
    # Define path length data.frame
    path_df <- data.frame(name=V(graph)$name, stringsAsFactors=FALSE)
    
    # Get indices of excluded vertices
    skip_idx <- which(path_df$name == root)
    if (!is.null(field)) {
        g <- get.vertex.attribute(graph, name=field)
        skip_idx <- union(skip_idx, which(g %in% exclude))
    }
    
    # Get paths
    step_list <- get.shortest.paths(graph, root, mode="out", weights=NA, output="vpath")
    step_list <- step_list$vpath
    
    # Get path lengths
    for (i in 1:length(step_list)) {
        v <- step_list[[i]]
        path_df[i, "steps"] <- sum(!(v %in% skip_idx)) 
        path_df[i, "distance"] <- sum(E(graph, path=v)$weight)
    }
    
    return(path_df)
}


#' Retrieve the first non-root node of a lineage tree
#' 
#' \code{getFounder} returns the set of founding node(s) for an igraph object. The founding 
#' node(s) are defined as the set of nodes within the minimum weighted or unweighted path length
#' from the root (germline) of the lineage tree, allowing for exclusion of specific groups of
#' nodes.
#'
#' @param    graph    igraph object with vertex annotations.
#' @param    path     string defining whether to use unweighted (steps) or weighted (distance) 
#'                    measures for determining the founder node set. Must be one of 
#'                    \code{c("steps", "distance")}. 
#' @param    root     name of the root (germline) node.
#' @param    field    annotation field to use for both unweighted path length exclusion and
#'                    consideration as a founder node. if \code{NULL} do not exclude any nodes.
#' @param    exclude  vector of annotation values in the given field to exclude from potential 
#'                    founder set. If \code{NULL} do not exclude any nodes. Has no effect if 
#'                    \code{field=NULL}.
#' @return   A data.frame of the founder node(s).
#' 
#' @seealso  \code{\link{getPathLengths}}.
#' @examples
#' # Define simple graph
#' library(igraph)
#' graph <- graph(c(1,2,2,3,2,4,3,5,3,6), directed=TRUE)
#' V(graph)$name <- c("Germline", "Inferred", "Seq1", "Seq2", "Seq3", "Seq4")
#' V(graph)$isotype <- c(NA, NA, "IgM", "IgG", "IgA", "IgA")
#' V(graph)$label <- c("Germline", "Inferred", "IgM", "IgG", "IgA", "IgA")
#' E(graph)$weight <- c(10, 3, 6, 4, 1)
#' E(graph)$label <- E(graph)$weight
#'
#' Plot graph with a tree layout
#' ly <- layout.reingold.tilford(graph, root="Germline", circular=F, flip.y=T)
#' plot(graph, layout=ly)
#' 
#' # Use unweighted path length and do not exclude any nodes
#' getFounder(graph, path="steps", root="Germline")
#'
#' # Exclude nodes without an isotype annotation and use weighted path length
#' getFounder(graph, path="distance", root="Germline", field="isotype", exclude=NA)
#' 
#' @export
getFounder <- function(graph, path="distance", root="Germline", field=NULL, exclude=NULL) {
    # Get distance from root
    path_df <- getPathLengths(graph, root=root, field=field, exclude=exclude)
    
    # Get indices of excluded vertices
    skip_idx <- which(path_df$name == root)
    if (!is.null(field)) {
        g <- get.vertex.attribute(graph, name=field)
        skip_idx <- union(skip_idx, which(g %in% exclude))
    }
    
    # Get founder nodes
    if (path == "distance") { 
        path_len <- setNames(path_df$distance, 1:nrow(path_df))
    } else if (path == "steps") {
        path_len <- setNames(path_df$steps, 1:nrow(path_df))
    } else {
        stop("Invalid value for 'path' parameter. Must be one of c('distance', 'steps').\n")
    }
    path_len <- path_len[-skip_idx]
    root_idx <- as.numeric(names(path_len)[which(path_len == min(path_len))])
    root_df <- get.data.frame(graph, what="vertices")[root_idx, ]
    root_df$steps <- path_df$steps[root_idx]
    root_df$distance <- path_df$distance[root_idx]
    
    return(root_df)
}