# SHMulate

#' @include shazam.R
NULL

#' Creates an igraph object from edge and node data.frames
#'
#' @param   edge.df   data.frame of edges with columns [from, to, weight]
#' @param   node.df   data.frame of nodes with columns [taxa, seq]
#' @return  an igraph object with added vertex annotations
#' @export
getGraph <- function(edge_df, node_df) {
  # Create igraph object
  g_obj <- graph.data.frame(edge_df, directed=T)
  V(g_obj)$number <- match(V(g_obj)$name, node_df$taxa)
  sub_df <- node_df[V(g_obj)$number, ]
  # Add annotations from node_df to graph
  for(col in names(sub_df)) {
    g_obj <- set.vertex.attribute(g_obj, name=col, value=sub_df[,col])
  }
  return(g_obj)
}
