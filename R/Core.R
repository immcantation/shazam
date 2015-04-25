#### Plotting functions ####

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


#### Transformation functions ####

# Converts a matrix to a vector
#
# \code{clearConsole} clears the console.
# 
# @examples
# # Generate a sampel mutations_array 
# sample_matrix <- matrix(sample(20,4),nrow=2, dimnames=list( c("CDR","FWR"), c("R","S") ))
# collapseMatrixToVector(sample_matrix)
#
collapseMatrixToVector <- function(mat, byrow = FALSE){    
    # Get the row and column names
    rnames <- rownames(mat)
    cnames <- colnames(mat)
    if (is.null(rnames)) { rnames <- paste0("Row", 1:nrow(mat)) }
    if (is.null(cnames)) { cnames <- paste0("Column", 1:ncol(mat)) }
    # Combine the row and columns names
    cobminedNames <- outer(rnames, cnames, paste, sep = "_")
    
    # Collapse the matrix to a vector
    if (byrow) {
        collapsed_mat <- c(t(mat))
        names(collapsed_mat) <- c(t(cobminedNames))
    }
    else{
        collapsed_mat <- c(mat)
        names(collapsed_mat) <- c(cobminedNames)
    }
    return(collapsed_mat)
}

#### Logging functions ####

# Make a progress bar when using %dopar% with a parallel backend
# Adapated from http://stackoverflow.com/questions/5423760/how-do-you-create-a-progress-bar-when-using-the-foreach-function-in-r
#
# @return  NULL
doparProgressBar <- function(n){
    pb <- txtProgressBar(min=1, max=n-1,style=3)
    count <- 0
    function(...) {
        count <<- count + length(list(...)) - 1
        setTxtProgressBar(pb,count)
        #Sys.sleep(0.01)
        flush.console()
        rbind(...)
        
    }
}


#### OS interaction functions ####

#' Determines the OS plotform being used
#'
#' @return  The OS platform
#' 
#' @examples
#' getPlatform()
#' 
#' @export
getPlatform <-function(){
    return(.Platform$OS.typ)
}


#' Determines the (maximum) numbers of cores/cpus availalbe
#'
#' @return  The number of cores/cpus available. Returns 1 if undeterminable.
#'
#' @examples
#' getnproc()
#' 
#' @export
getnproc <-function(){
    platform <- getPlatform()
    nproc <- 1
    if(platform=="windows") nproc <- Sys.getenv('NUMBER_OF_PROCESSORS')
    if(platform=="unix") nproc <- system("nproc", intern=TRUE)
    return(as.numeric(nproc))
}


#' Clears the console
#'
#' \code{clearConsole} clears the console.
#' 
#' @examples
#' clearConsole()
#'
#' @export
clearConsole <- function(){
    platform <- getPlatform()
    nproc <- 1
    if(platform=="windows") cat("\014")
    if(platform=="unix") system('clear')
}