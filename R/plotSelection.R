#' Plot the BASELINe selection values
#'
#'
#' @param    db              data.frame returned by calcGroupedBaseline
#' @param    groups          column names by which you want the anlaysis grouped by
#' @param    colors          values are colors to assign to those group values
#' @param    main_title      string specifying the plot title.
#' @param    legend_title    string specifying the legend title.
#' @return   a \code{ggplot} object
#' @export
plotSelection <- function(db, groups=NULL, colors=NULL,
                          main_title="Selection Strength",
                          legend_title=NULL) {
  # Define plot elements
  plotCDR <-  ggplot(db, aes_string(x=groups, y="CDR_Sigma")) +
              theme_bw() +
              geom_errorbar(aes_string(ymin="CDR_CI95_Lower", ymax="CDR_CI95_Upper",colors=groups), width=.1) +
              geom_point(aes_string(color=groups), size=3)
  plot(plotCDR)
#
#     ggplot(db, aes(CDR_Sigma)) +
#     ggtitle(main_title) +
#     #getBaseTheme() +
#     xlab('Density') +
#     ylab('Selection Strength') +
#     #geom_ribbon(aes(ymin=CDR_CI95_Lower, ymax=CDR_CI95_Upper, fill=groups), alpha=0.25) +
#     geom_line(aes(color=group)) +
#     guides(color=guide_legend(title=legend_title),
#            fill=guide_legend(title=legend_title))
#   if (!is.null(colors)) {
#     p1 <- p1 + scale_color_manual(name=legend_title, values=colors) +
#       scale_fill_manual(name=legend_title, values=colors)
#   }
  return(plotCDR)
}
