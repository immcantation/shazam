#' Plot the BASELINe selection values
#'
#'
#' @param    db              data.frame returned by calcGroupedBaseline
#' @param    groups          column names by which you want the anlaysis grouped by
#' @param    colors          values are colors to assign to those group values
#' @param    outputName         Filename for the output files
#' @param    outputPath         Path to save output
#' @export
plotSelection <- function(db, groups=NULL, colors=NULL,outputName=NULL, outputPath=NULL) {

  figurePath = file.path(outputPath, paste0(outputName, ".pdf"))
  cat(figurePath)
  plotCDR <-  ggplot(db, aes_string(x=groups, y="CDR_Sigma")) +
              theme_bw() +
              theme( axis.title.x = element_text(face="bold", colour="#990000", size=16),
                     axis.title.y = element_text(face="bold", colour="#990000", size=16),
                     axis.text.x  = element_text(angle=0, size=12),
                     axis.text.y  = element_text(angle=0, size=12),
                     title = element_text(angle=0, vjust=0.5, size=20))+
              theme(legend.title = element_text(size = 8))+
              geom_errorbar(aes_string(ymin="CDR_CI95_Lower", ymax="CDR_CI95_Upper",colors=groups), width=.1) +
              geom_point(aes_string(color=groups), size=5)+
              ylab("CDR")+
              ggtitle("Selection Strength")+
              scale_y_continuous(limits=c(-2,1), breaks=seq(-1, 1, 0.5), oob=rescale_none)


  plotFWR <-  ggplot(db, aes_string(x=groups, y="FWR_Sigma")) +
    theme_bw() +
    theme( axis.title.x = element_text(face="bold", colour="#990000", size=16),
           axis.title.y = element_text(face="bold", colour="#990000", size=16),
           axis.text.x  = element_text(angle=0,  size=12),
           axis.text.y  = element_text(angle=0,  size=12),
           title = element_text(angle=0, vjust=0.5, size=20))+
    theme(legend.title = element_text(size = 8))+
    geom_errorbar(aes_string(ymin="FWR_CI95_Lower", ymax="FWR_CI95_Upper",colors=groups), width=.1) +
    geom_point(aes_string(color=groups), size=5)+
    ylab("FWR")+
    scale_y_continuous(limits=c(-3,1), breaks=seq(-3, 1, 0.5), oob=rescale_none)


  pdf(figurePath,8,10)
  multiplot(plotCDR, plotFWR, cols=1)
  dev.off()
}
