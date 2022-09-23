# Deprecated and defunct functions

#' @include MutationProfiling.R
NULL

#### Deprecated ####

#' slideWindowTunePlot - plotSlideWindowTune backward compatability
#'
#' Wrapper function for \link{plotSlideWindowTune}
#' 
#' @param    tuneList            a list of logical matrices returned by \link{slideWindowTune}.
#' @param    plotFiltered        whether to plot the number of filtered (\code{TRUE} or \code{filtered}), 
#'                               or remaining (FALSE or remaining) sequences for each mutation threshold. 
#'                               Use \code{NULL} or \code{per_mutation} to plot the number of sequences 
#'                               at each mutation value. Default is \code{TRUE}.
#' @param    percentage          whether to plot on the y-axis the percentage of filtered sequences
#'                               (as opposed to the absolute number). Default is \code{FALSE}.                             
#' @param    jitter.x            whether to jitter x-axis values. Default is \code{FALSE}.                               
#' @param    jitter.x.amt        amount of jittering to be applied on x-axis values if 
#'                               \code{jitter.x=TRUE}. Default is 0.1.
#' @param    jitter.y            whether to jitter y-axis values. Default is \code{FALSE}.
#' @param    jitter.y.amt        amount of jittering to be applied on y-axis values if 
#'                               \code{jitter.y=TRUE}. Default is 0.1.                               
#' @param    pchs                point types to pass on to \link{plot}.
#' @param    ltys                line types to pass on to \link{plot}.
#' @param    cols                colors to pass on to \link{plot}.                             
#' @param    plotLegend          whether to plot legend. Default is \code{TRUE}.
#' @param    legendPos           position of legend to pass on to \link{legend}. Can be either a
#'                               numeric vector specifying x-y coordinates, or one of 
#'                               \code{"topright"}, \code{"center"}, etc. Default is \code{"topright"}.
#' @param    legendHoriz         whether to make legend horizontal. Default is \code{FALSE}.
#' @param    legendCex           numeric values by which legend should be magnified relative to 1.
#' @param    title               plot main title. Default is NULL (no title)
#' @param    returnRaw           Return a data.frame with sequence counts (TRUE) or a
#'                               plot. Default is \code{FALSE}.
#' 
#' @details  For each \code{windowSize}, if \code{plotFiltered=TRUE}, the x-axis 
#'           represents a mutation threshold range, and the y-axis the number of
#'           sequences that have at least that number of mutations. If 
#'           \code{plotFiltered=TRUE}, the y-axis represents the number of sequences
#'           that have less mutations than the mutation threshold range. For the same
#'           window size, a sequence can be included in the counts for different
#'           mutation thresholds. For example, sequence "CCACCAAAA" with germline
#'           "AAAAAAAAA" has 4 mutations. This sequence has at least 2 mutations 
#'           and at least 3 mutations, in a window of size 4. the sequence will
#'           be included in the sequence count for mutation thresholds 2 and 3.
#'           If \code{plotFiltered=TRUE}, the sequences are counted only once for
#'           each window size, at their largest mutation threshold. The above 
#'           example sequence would be included in the sequence count for 
#'           mutation threshold 3. 
#'           
#'           When plotting, a user-defined \code{amount} of jittering can be applied on values plotted
#'           on either axis or both axes via adjusting \code{jitter.x}, \code{jitter.y}, 
#'           \code{jitter.x.amt} and \code{jitter.y.amt}. This may be help with visually distinguishing
#'           lines for different window sizes in case they are very close or identical to each other. 
#'           If plotting percentages (\code{percentage=TRUE}) and using jittering on the y-axis values 
#'           (\code{jitter.y=TRUE}), it is strongly recommended that \code{jitter.y.amt} be set very
#'           small (e.g. 0.01). 
#'           
#'           \code{NA} for a combination of \code{mutThresh} and \code{windowSize} where 
#'           \code{mutThresh} is greater than \code{windowSize} will not be plotted. 
#' 
#' @seealso  See \link{slideWindowTune} for how to get \code{tuneList}. See \link{jitter} for 
#'           use of \code{amount} of jittering.
#' 
#' @examples
#' # Use an entry in the example data for input and germline sequence
#' data(ExampleDb, package="alakazam")
#' 
#' # Try out thresholds of 2-4 mutations in window sizes of 3-5 nucleotides 
#' # on a subset of ExampleDb
#' tuneList <- slideWindowTune(db = ExampleDb[1:10, ], 
#'                            mutThreshRange = 2:4, windowSizeRange = 3:5,
#'                            verbose = FALSE)
#'
#' # Visualize
#' # Plot numbers of sequences filtered without jittering y-axis values
#' slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3, 
#'                     plotFiltered=TRUE, jitter.y=FALSE)
#'                     
#' # Notice that some of the lines overlap
#' # Jittering could help
#' slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3,
#'                     plotFiltered=TRUE, jitter.y=TRUE)
#'                     
#' # Plot numbers of sequences remaining instead of filtered
#' slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3, 
#'                     plotFiltered=FALSE, jitter.y=TRUE, 
#'                     legendPos="bottomright")
#'                     
#' # Plot percentages of sequences filtered with a tiny amount of jittering
#' slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3,
#'                     plotFiltered=TRUE, percentage=TRUE, 
#'                     jitter.y=TRUE, jitter.y.amt=0.01)
#' @export
slideWindowTunePlot <- function(tuneList, 
                                plotFiltered = c(TRUE,FALSE,NULL,'filtered','remaining','per_mutation'), 
                                percentage = FALSE,
                                jitter.x = FALSE, jitter.x.amt = 0.1,
                                jitter.y = FALSE, jitter.y.amt = 0.1,
                                pchs = 1, ltys = 2, cols = 1,
                                plotLegend = TRUE, legendPos = "topright", 
                                legendHoriz = FALSE, legendCex = 1, title=NULL,
                                returnRaw=FALSE){
    .Deprecated("plotSlideWindowTune",
                msg="slideWindowTunePlot() is deprecated, please see plotSlideWindowTune() for future use")

    # input validation
    plotFiltered_choices <- c(TRUE,FALSE,NULL,'filtered','remaining','per_mutation')
    plotFiltered <- plotFiltered[1]
    if (!is.null(plotFiltered)) {
        if (!plotFiltered %in% plotFiltered_choices) {
            stop("`plotFiltered` must be one of: ", paste(plotFiltered_choices,collapse=", "))
        }   
    }
    
    # logic for converting T/F/NULL to new values
    if (plotFiltered %in% c(TRUE,'filtered')) {
        plotFilteredMapped <- 'filtered'
    } else if (plotFiltered %in% c(FALSE,'remaining')) {
        plotFilteredMapped <- 'remaining'
    } else {
        plotFilteredMapped <- 'per_mutation'
    }
    
    plotSlideWindowTune(tuneList, plotFiltered = plotFilteredMapped, 
                        percentage = percentage,
                        jitter.x = jitter.x, jitter.x.amt = jitter.x.amt,
                        jitter.y = jitter.y, jitter.y.amt = jitter.y.amt,
                        pchs = pchs, ltys = ltys, cols = cols,
                        plotLegend = plotLegend, legendPos = legendPos, 
                        legendHoriz = legendHoriz, legendCex = legendCex, title=title,
                        returnRaw=returnRaw)
}