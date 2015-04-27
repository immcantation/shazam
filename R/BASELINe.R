#' @include RegionDefinitions.R
#' @include shm.R
NULL

#### Classes ####

#' S4 class defining a Baseline (selection) object
#' 
#' \code{Baseline} defines a common data structure for defining the BASELINe Selection
#' results
#' 
#' @slot    description         General information regarding the sequences, selection 
#'                              analysis and/or object.
#' @slot    db                  \code{data.frame} containing annotation information about 
#'                              the sequences/selection results.
#' @slot    regionDefinition    \code{\link{RegionDefinition}} object defining the regions
#'                              and boundaries of the Ig sequences. Note, only the part of
#'                              sequences defined in \code{regionDefinition} are analyzed.
#'                              Any mutations outside the definition will be ignored. E.g.
#'                              If the default \code{\link{IMGT_V_NO_CDR3}} definition is
#'                              used, then mutations in positions greater than 312 will not
#'                              be counted.                   
#' @slot    testStatistic       The statistical framework used to test for selection.
#'                              E.g.
#'                              \code{local} = CDR_R / (CDR_R + CDR_S)
#'                              \code{focused} = CDR_R / (CDR_R + CDR_S + FWR_S).
#'                              For \code{focused} the \code{regionDefinition} must only
#'                              contain two regions. If more than two regions are defined
#'                              the \code{local} test statistic will be used.
#'                              See Uduman et al. (2011) for further information.
#' @slot    regions             \code{character} vector of the regions the BASELINe 
#'                              selection was carried out on. E.g. "CDR & "FWR" or
#'                              "CDR1", "CDR2", "CDR3" etc.
#' @slot    numbOfSeqs          \code{matrix} of r x c dimensions, where:
#'                              r = number of rows = number of groups/sequencs and
#'                              c = number of columns = number of regions
#'                              The matrix contains the number of non-NA sequences or the
#'                              number of BASELINe posterior probability distribution 
#'                              functions (PDF) that were convoluted (grouped) together. 
#'                              This information is essential to further regroup PDFs and 
#'                              weight them accordingly.
#'                              
#' @slot    pdfs                A \code{list} (one item for each defined region e.g. CDR &
#'                              FWR) of \code{matrices} of r x c deminetions, where:
#'                              r = number of sequences or groups and
#'                              c = number of columns = length of the PDF (default 4001)
#'                              The matrix contains the BASELINe posterior probability 
#'                              distribution functions (PDF) for each sequence or goup of
#'                              sequences. 
#' @slot    stats               A melted \code{data.frame} of BASELINe statistics, 
#'                              including: Selection strength (Sigma), 95\% confidence 
#'                              intervals and P values.
#'                                
#'                          
#' @name Baseline
#' @export
setClass("Baseline", 
         slots = c( description="character",
                    db="data.frame",
                    regionDefinition="RegionDefinition",
                    testStatistic="character",
                    regions="character",
                    numbOfSeqs="matrix",
                    pdfs="list",
                    stats="data.frame"
         )
)

#### Baseline building functions #####

#' Creates a Baseline object
#' 
#' \code{createBaseline} creates a \code{RegionDefinition}.
#'
#' @param   description         General information regarding the sequences, selection 
#'                              analysis and/or object.
#' @param   db                  \code{data.frame} containing annotation information about 
#'                              the sequences/selection results.
#' @param   regionDefinition    \code{\link{RegionDefinition}} object defining the regions
#'                              and boundaries of the Ig sequences. Note, only the part of
#'                              sequences defined in \code{regionDefinition} are analyzed.
#'                              Any mutations outside the definition will be ignored. E.g.
#'                              If the default \code{\link{IMGT_V_NO_CDR3}} definition is
#'                              used, then mutations in positions greater than 312 will not
#'                              be counted.                   
#' @param   testStatistic       The statistical framework used to test for selection.
#'                              E.g.
#'                              \code{local} = CDR_R / (CDR_R + CDR_S)
#'                              \code{focused} = CDR_R / (CDR_R + CDR_S + FWR_S).
#'                              For \code{focused} the \code{regionDefinition} must only
#'                              contain two regions. If more than two regions are defined
#'                              the \code{local} test statistic will be used.
#'                              See Uduman et al. (2011) for further information.
#' @param   numbOfSeqs          \code{matrix} of r x c dimensions, where"
#'                              r = number of rows = number of groups/sequencs and
#'                              c = number of columns = number of regions
#'                              The matrix contains the number of non-NA sequences or the
#'                              number of BASELINe posterior probability distribution 
#'                              functions (PDF) that were convoluted (grouped) together. 
#'                              This information is essential to further regroup PDFs and 
#'                              weight them accordingly.
#' @param   pdfs                 A \code{list} (one item for each defined region e.g. CDR &
#'                              FWR) of \code{matrices} of r x c deminetions, where:
#'                              r = number of sequences or groups and
#'                              c = number of columns = length of the PDF (default 4001)
#'                              The matrix contains the BASELINe posterior probability 
#'                              distribution functions (PDF) for each sequence or goup of
#'                              sequences.
#' @param   stats               A melted \code{data.frame} of BASELINe statistics, 
#'                              including: Selection strength (Sigma), 95\% confidence 
#'                              intervals and P values.#'                                
#' 
#' @return   A \code{Baseline} object.
#' 
#' @seealso  See \code{\link{Baseline}} for the return object.
#' 
#' @examples
#' library(shm)
#' 
#' 
#' @export
createBaseline <- function( description="",
                            db="",
                            regionDefinition="",
                            testStatistic="",
                            regions="",
                            numbOfSeqs="",
                            pdfs="",
                            stats="") {
    
    if ( stats=="" ) {
        stats <- data.frame( GROUP=character(),
                             REGION=character(),
                             BASELINE_SIGMA=character(),
                             BASELINE_CI_LOWER=character(),
                             BASELINE_CI_UPPER=character(),
                             BASELINE_CI_PVALUE=character(),
                             stringsAsFactors=FALSE) 
    }
    # Define RegionDefinition object
    baseline <- new( "Baseline",
                     description=description,
                     db=db,
                     regionDefinition=regionDefinition,
                     testStatistic=testStatistic,
                     regions=regionDefinition@regions,
                     numbOfSeqs=numbOfSeqs,
                     pdfs=pdfs,
                     stats=stats)
    
    return(baseline)
}


# Edit the Baseline object
# 
# \code{editBaseline} edits a \code{Baseline}.
#
# param   baseline     The \code{Baseline} S4 object to be edited.
# param   field_name   Name of the field in the \code{Baseline} S4 object to be edited.
# param   value        The value to set the \code{field_name}.
# 
# return   A \code{Baseline} object.
# 
# seealso  See \code{\link{Baseline}} for the return object.
editBaseline <- function ( baseline,
                           field_name,
                           value ) {
    if ( !match(field_name, slotNames(baseline)) ) { stop("field_name not part of BASELINe object!")}
    slot(baseline, field_name) = value
    return(baseline)
}


#### Baseline selection calculating functions ####

#' Calculate the BASELINe PDFs
#' 
#' \code{calcBaselinePdfs} calculates the BASELINe posterior probability distribution 
#' functions (PDFs) for sequences in the given db
#'
#' @param   db                  \code{data.frame} containing sequence data and annotation.
#' @param   sequenceColumn      Name of the column containing sample/input sequences.
#' @param   germlineColumn      Name of the column containing germline sequences.
#' @param   testStatistic       The statistical framework used to test for selection.
#'                              E.g.
#'                              \code{local} = CDR_R / (CDR_R + CDR_S)
#'                              \code{focused} = CDR_R / (CDR_R + CDR_S + FWR_S).
#'                              For \code{focused} the \code{regionDefinition} must only
#'                              contain two regions. If more than two regions are defined
#'                              the \code{local} test statistic will be used.
#'                              See Uduman et al. (2011) for further information.
#' @param   regionDefinition    \code{\link{RegionDefinition}} object defining the regions
#'                              and boundaries of the Ig sequences. Note, only the part of
#'                              sequences defined in \code{regionDefinition} are analyzed.
#'                              Any mutations outside the definition will be ignored. E.g.
#'                              If the default \code{\link{IMGT_V_NO_CDR3}} definition is
#'                              used, then mutations in positions greater than 312 will not
#'                              be counted.
#' @param   nproc               number of cores to distribute the operation over. If 
#'                              \code{nproc} = 0 then the \code{cluster} has already been
#'                              set and will not be reset.
#' 
#' @return  A \code{Baseline} object, containing the modified \code{db} and the BASELINe 
#'          posterior probability distribution functions (PDF) for each of the sequences
#'           
#' @details Calculates the BASELINe posterior probability distribution function (PDF) for 
#'          sequences in the provided \code{db}. 
#'          
#'          If the \code{db} does not contain the 
#'          required columns to calculate the PDFs (namely OBSERVED & EXPECTED mutations)
#'          then the function will:
#'          1. Collapse the sequences by the CLONE column (if present)
#'          2. Calculate the numbers of observed mutations
#'          3. Calculate the expected frequencies of mutations
#'          and modify the provided \code{db} (this will be returned as part of the 
#'          \code{Baseline} object).
#'          
#'          
#' @seealso To calculate BASELINe statistics, such as the mean selection strength
#'          and the 95\% confidence interval, see .
#'          To group the sequence PDFs and get a combined PDF see 
#'          \code{\link{groupBaseline}}.
#' 
#' 
#' @references
#' \enumerate{
#'   \item  Gur Yaari; Mohamed Uduman; Steven H. Kleinstein. Quantifying selection 
#'          in high-throughput Immunoglobulin sequencing data sets. Nucleic Acids Res.
#'           2012 May 27. 
#'  \item   Mohamed Uduman; Gur Yaari; Uri Hershberg; Mark J. Shlomchik; Steven H. 
#'          Kleinstein. Detecting selection in immunoglobulin sequences. Nucleic Acids 
#'          Res. 2011 Jul;39(Web Server issue):W499-504.  
#'  \item   Hershberg U, Uduman M, Shlomchik MJ, Kleinstein SH. Improved methods for
#'          detecting selection by mutation analysis of Ig V region sequences. Int Immunol.
#'          2008 May;20(5):683-94. doi: 10.1093/intimm/dxn026. Epub 2008 Apr 7. 
#' }
#' 
#' @examples
#' # Load example data
#' library("shm")
#' dbPath <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(dbPath)
#'                      
#' @export
calcBaselinePdfs <- function( db,
                              sequenceColumn="SEQUENCE_IMGT",
                              germlineColumn="GERMLINE_IMGT_D_MASK",
                              testStatistic=c("local","focused"),
                              regionDefinition=IMGT_V_NO_CDR3,
                              nproc=1 ) {
    
    # Evaluate argument choices
    testStatistic <- match.arg(testStatistic, c("local", "focused"))
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, getnproc())
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    if(nproc>1){        
        cluster <- makeCluster(nproc, type = "SOCK")
        clusterExport( cluster, list( 'db',
                                      'sequenceColumn', 'germlineColumn', 
                                      'regionDefinition',
                                      'break2chunks', 'PowersOfTwo', 
                                      'convolutionPowersOfTwo', 
                                      'convolutionPowersOfTwoByTwos', 
                                      'weighted_conv', 
                                      'calculate_bayesGHelper', 
                                      'groupPosteriors', 'fastConv',
                                      'calcBaselinePdfs_Helper'), 
                       envir=environment() )
        clusterEvalQ(cluster, library(shm))
        clusterEvalQ(cluster, library(seqinr))
        registerDoSNOW(cluster)
    } else if( nproc==1 ) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
    
    # If db does not contain the required columns to calculate the PDFs (namely OBSERVED 
    # & EXPECTED mutations), then the function will:
    #          1. Collapse the sequences by the CLONE column (if present)
    #          2. Calculate the numbers of observed mutations
    #          3. Calculate the expected frequencies of mutations    
    # After that BASELINe prob. densities can be calcualted per seqeunce.    
    observedColumns <- paste0("OBSERVED_", regionDefinition@labels)
    expectedColumns <- paste0("EXPECTED_", regionDefinition@labels)
    
    if( !all( c(observedColumns,expectedColumns) %in% colnames(db) ) ) {
        
        # If the germlineColumn & sequenceColumn are not found in the db error and quit
        if( !all( c(sequenceColumn, germlineColumn) %in% colnames(db) ) ) {
            stop( paste0("Both ", sequenceColumn, " & ", germlineColumn, 
                         " columns need to be present in the db") )
        }
        
        # Collapse the sequences by the CLONE column (if present)
        if( "CLONE" %in% colnames(db) ) {                       
            db <- getClonalConsensus( db, 
                                      cloneColumn="CLONE", 
                                      sequenceColumn=sequenceColumn,
                                      germlineColumn=germlineColumn,
                                      collapseByClone=TRUE, nproc=cluster)            
            sequenceColumn="CLONAL_CONSENSUS_SEQUENCE"
        }
        
        # Calculate the numbers of observed mutations
        db <- getObservedMutations( db,
                                    sequenceColumn=sequenceColumn,
                                    germlineColumn="GERMLINE_IMGT_D_MASK",
                                    regionDefinition=IMGT_V_NO_CDR3,
                                    nproc=0 )
        
        # Calculate the expected frequencies of mutations
        db <- getExpectedMutationFrequencies( db,
                                              sequenceColumn="CLONAL_CONSENSUS_SEQUENCE",
                                              germlineColumn="GERMLINE_IMGT_D_MASK",
                                              regionDefinition=IMGT_V_NO_CDR3,
                                              nproc=0 )
    }
    
    # Calculate PDFs for each sequence
    
    # Print status to console
    cat("Calculating BASELINe probability density functions...\n")
    
    # Number of sequences (used in foreach)
    totalNumbOfSequences <- nrow(db)
    # The column indexes of the OBSERVED_ and EXPECTED_
    cols_observed <- grep( paste0("OBSERVED_"),  colnames(db) ) 
    cols_expected <- grep( paste0("EXPECTED_"),  colnames(db) ) 
    
    # Exporting additional environment variables and functions needed to run foreach 
    if( nproc!=1 ) {
        clusterExport( 
            cluster, list('cols_observed', 'cols_expected'), 
            envir=environment() 
        )
        registerDoSNOW(cluster)
    }
    
    list_pdfs <- list()
    regions <- regionDefinition@regions
    # For every region (e.g. CDR, FWR etc.)
    for (region in regions) {
        
        # Foreach returns a list of PDFs
        list_region_pdfs <- 
            foreach( i=icount(totalNumbOfSequences) ) %dopar% {                
                calcBaselinePdfs_Helper( 
                    observed = db[i,cols_observed],
                    expected = db[i,cols_expected],
                    region = region,
                    testStatistic = testStatistic,
                    regionDefinition = regionDefinition
                )
            }
        # Convert the list of the region's PDFs into a matrix                
        list_pdfs[[region]] <- 
            do.call( rbind, 
                     lapply( 
                         list_region_pdfs, 
                         function(x) { 
                             length(x) <- 4001 
                             return(x)
                         }
                     )
            ) 
    }
    
    # Initialize numbOfSeqs
    # This holds the number of non NA sequences
    numbOfSeqs <- matrix( NA, 
                          ncol=length(regions), 
                          nrow=1,
                          dimnames=list( 1, regions )
    )
    
    # Create a Baseline object with the above results to return
    baseline <- createBaseline( description="",
                                db=db,
                                regionDefinition=regionDefinition,
                                testStatistic=testStatistic,
                                regions=regionDefinition@regions,
                                numbOfSeqs=numbOfSeqs,
                                pdfs=list_pdfs )
    
    # Stop SNOW cluster
    if(nproc > 1) { stopCluster(cluster) }
    
    return(baseline)
    
}



# Helper function for calcBaselinePdfs
#
# \code{calcBaselinePdfs_Helper} calculates 
#
# @param   observed
# @param   expected
# @param   region
# @param   testStatistic
# @param   regionDefinition
# 
# @return  A modified \code{\link{Baseline}} object with the BASELINe probability 
#          density function calculated for the regions defined in the \code{regionDefinition}.
#           
# @details
# \code{getBaseline} calculates 
# 
# @seealso  
# See also
# @examples
# # Load example data
# library("shm")
# dbPath <- system.file("extdata", "Influenza_IB.tab", package="shm")
# db <- readChangeoDb(dbPath)
#    
# export
calcBaselinePdfs_Helper  <- function( observed,
                                      expected,
                                      region,
                                      testStatistic="local",
                                      regionDefinition=IMGT_V_NO_CDR3 ) {
    
    # Evaluate argument choices
    testStatistic <- match.arg(testStatistic, c("local", "focused"))
    
    #If there are more than two regions (e.g. CDR and FWR then you cannot perform the focused test)
    if (testStatistic=="focused" & length(regionDefinition@regions)!=2) {
        testStatistic="local"    
    }    
    
    # local test statistic
    if (testStatistic == "local") { 
        obsX_Index <- grep( paste0("OBSERVED_", region,"_R"),  names(observed) )
        obsN_Index <- grep( paste0("OBSERVED_", region),  names(observed) )
        
        expX_Index <- grep( paste0("EXPECTED_", region,"_R"),  names(expected) )
        expN_Index <- grep( paste0("EXPECTED_", region),  names(expected) )       
    }
    
    # focused test statistic
    if (testStatistic == "focused") { 
        obsX_Index <- grep( paste0("OBSERVED_", region,"_R"),  names(observed) )
        obsN_Index <- 
            grep( 
                paste0( 
                    "OBSERVED_", region, "|", 
                    "OBSERVED_", regionDefinition@regions[regionDefinition@regions!=region], "_S"
                ),
                names(observed) 
            )
        
        expX_Index <- grep( paste0("EXPECTED_", region,"_R"),  names(expected) )
        expN_Index <- 
            grep( 
                paste0( 
                    "EXPECTED_", region, "|", 
                    "EXPECTED_",  regionDefinition@regions[regionDefinition@regions!=region], "_S"
                ),
                names(expected) 
            )        
    }     
    
    obsX <- as.numeric( observed[obsX_Index] )
    obsN <- as.numeric( sum( observed[obsN_Index], na.rm=T ) )
    
    expP <-
        as.numeric( 
            expected[expX_Index] / 
                sum( expected[expN_Index], na.rm=T )
        )
    
    return( calcBaselineBinomialPdf( x=obsX, n=obsN, p=expP) )
}

# Calculate the BASELINe probability function in a
# binomial framework.
calcBaselineBinomialPdf <- function ( x=3, 
                                      n=10, 
                                      p=0.33,
                                      CONST_i=CONST_I,
                                      max_sigma=20,
                                      length_sigma=4001 ) {
    if(n!=0){
        sigma_s<-seq(-max_sigma,max_sigma,length.out=length_sigma)
        sigma_1<-log({CONST_i/{1-CONST_i}}/{p/{1-p}})
        index<-min(n,60)
        y<- dbeta(CONST_i,x+BAYESIAN_FITTED[index],n+BAYESIAN_FITTED[index]-x)*(1-p)*p*exp(sigma_1)/({1-p}^2+2*p*{1-p}*exp(sigma_1)+{p^2}*exp(2*sigma_1))
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


#' Group BASELINe PDFs
#' 
#' \code{groupBaseline} convolutes the BASELINe posterior probability distribution 
#' functions (PDFs) of sequences to get a combined (grouped) pdf
#'
#' @param   baseline            \code{Baseline} object, containing the \code{db} and the 
#'                              BASELINe posterior probability distribution functions 
#'                              (PDF) for each of the sequences. This would be returned by
#'                              \code{\link{calcBaselinePdfs}}.
#' @param   groupBy             The columns in the \code{db} slot of the \code{Baseline}
#'                              object to group the sequence PDFs by.
#' @param   nproc               number of cores to distribute the operation over. If 
#'                              \code{nproc} = 0 then the \code{cluster} has already been
#'                              set and will not be reset.
#' 
#' @return  A \code{Baseline} object, containing the modified \code{db} and the BASELINe 
#'          posterior probability distribution functions (PDF) for each of the groups.
#'           
#' @details Calculates the BASELINe posterior probability distribution function (PDF) for 
#'          sequences in the provided \code{db}, convoluted (grouped)  according to the 
#'          annotations in the \code{groupBy} argument.
#'          
#'          
#' @seealso To calculate BASELINe statistics, such as the mean selection strength
#'          and the 95\% confidence interval, see .
#' 
#' 
#' @references
#' \enumerate{
#'   \item  Gur Yaari; Mohamed Uduman; Steven H. Kleinstein. Quantifying selection 
#'          in high-throughput Immunoglobulin sequencing data sets. Nucleic Acids Res.
#'           2012 May 27. 
#' }
#' 
#' @examples
#' # Load example data
#' library("shm")
#' dbPath <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(dbPath)
#'                      
#' @export
groupBaseline <- function( baseline,
                           groupBy,
                           nproc=1 ) {
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, getnproc())

    # Convert the db (data.frame) to a data.table & set keys
    # This is an efficient way to get the groups of CLONES, instead of doing dplyr
    dt <- data.table(baseline@db)
    # Get the group indexes
    groupByFormatted <- paste(groupBy, collapse=",", sep=",")
    dt <- dt[ , list( yidx = list(.I) ) , by=groupByFormatted ]
    groups <- dt[,yidx] 
    df <- as.data.frame(dt)    
   

    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    if(nproc>1){        
        cluster <- makeCluster(nproc, type = "SOCK")
        clusterExport( cluster, list( 'baseline', 'groups',
                                      'break2chunks', 'PowersOfTwo', 
                                      'convolutionPowersOfTwo', 
                                      'convolutionPowersOfTwoByTwos', 
                                      'weighted_conv', 
                                      'calculate_bayesGHelper', 
                                      'groupPosteriors', 'fastConv'), 
                       envir=environment() )
        clusterEvalQ(cluster, library(shm))
        registerDoSNOW(cluster)
    } else if( nproc==1 ) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
 
    
    # Print status to console
    cat("Grouping BASELINe probability density functions...\n")
    
    # Number of total groups
    numbOfTotalGroups <- length(groups)
    list_pdfs <- list()
    regions <- baseline@regions
    
    # Initialize numbOfSeqs
    # This holds the number of non NA sequences
    numbOfSeqs <- matrix( NA, 
                          ncol=length(baseline@regions), 
                          nrow=numbOfTotalGroups,
                          dimnames=list( 1:numbOfTotalGroups, regions )
    )    
    

    # For every region (e.g. CDR, FWR etc.)
    for (region in regions) {
        
        list_region_pdfs  <-
            foreach( i=icount(numbOfTotalGroups)) %dopar% {
                matrix_GroupPdfs <- (baseline@pdfs[[region]])[groups[[i]],]
                
                list_GroupPdfs <- 
                    lapply( 1:nrow(matrix_GroupPdfs), 
                            function(rowIndex) {
                                rowVals <- matrix_GroupPdfs[rowIndex,]
                                if( !all(is.na(rowVals)) ) { matrix_GroupPdfs[rowIndex,] }
                            })
                
                list_GroupPdfs <- Filter(Negate(is.null), list_GroupPdfs)
                numbOfNonNASeqs <- length(list_GroupPdfs)
                
                # If all the sequences in the group are NAs, return a PDF of NAs
                if( length(list_GroupPdfs) == 0 ) { list_GroupPdfs = list((rep(NA,4001))) }
                
                return( c( groupPosteriors(list_GroupPdfs), numbOfNonNASeqs ) )
            }
        
        # Convert the list of the region's PDFs into a matrix                
        matrix_region_pdfs <- 
            do.call( rbind, 
                     lapply( 
                         list_region_pdfs, 
                         function(x) { 
                             length(x) <- 4002 
                             return(x)
                         }
                     )
            )
        
        
        list_pdfs[[region]] <- matrix_region_pdfs[,1:4001]
        #numbOfSeqs[,region] <- matrix_region_pdfs[,4002]
        #colnames(numbOfSeqs) <- paste0("NUMB_SEQUENCES_", colnames(numbOfSeqs))
    }
 
    
    # Create the db, which will now contain the group information
    db <- cbind( df[,groupBy], numbOfSeqs)
    if(!class(db)=="data.frame") { 
        db <- as.data.frame(db) 
        colnames(db)[1] <- groupBy
    }
    
    # Create a Baseline object with the above results to return
    baseline <- createBaseline( description="",
                                db=db,
                                regionDefinition=baseline@regionDefinition,
                                testStatistic=baseline@testStatistic,
                                regions=regions,
                                numbOfSeqs=numbOfSeqs,
                                pdfs=list_pdfs )
    
    # Stop SNOW cluster
    if(nproc > 1) { stopCluster(cluster) }
    
    return(baseline)
    
}


#' Calculate BASELINe statistics
#'
#' \code{calcBaselineStats} calculates BASELINe statistics such as the Selection Strength
#' (Sigma), the 95% confidence intervals & P-values.
#'
#' @param   baseline    \code{Baseline} object, containing the \code{db} and the 
#'                      BASELINe posterior probability distribution functions 
#'                      (PDF) for each of the sequences. This would be returned by
#'                      \code{\link{calcBaselinePdfs}}.
#' @param   nproc       number of cores to distribute the operation over. If 
#'                      \code{nproc} = 0 then the \code{cluster} has already been
#'                      set and will not be reset.
#' 
#' @return  A modified \code{Baseline} object with the BASELINe selection strength ,
#'          95\% confidence intervals and P-value. This information is updated in the 
#'          \code{stats} slot of the \code{baseline} object passed as an argument. 
#'           
#' @details \code{getBASELINeStats} calculates BASELINe statistics such as the Selection 
#'          Strength (Sigma), the 95% confidence intervals & P-values.
#' 
#' @seealso See \code{link{calcBaselinePdfs}} and \code{link{groupBaseline}}.
#' 
#' @examples
#' # Load example data
#' library("shm")
#' dbPath <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(dbPath)
#'                      
#' @export
calcBaselineStats <- function ( baseline,
                                nproc=1 ) {
    
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, getnproc())
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.  
    if(nproc>1){        
        cluster <- makeCluster(nproc, type = "SOCK")
        clusterExport( cluster, list( 'baseline',
                                      'calcBaselineSigma',
                                      'calcBaselineCI',
                                      'calcBaselinePvalue' ), 
                       envir=environment() )
        clusterEvalQ( cluster, library("shm") )
        registerDoSNOW(cluster)
    } else if( nproc==1 ) {
        # If needed to run on a single core/cpu then, regsiter DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    }
    
    # Printing status to console
    cat("Calculating BASELINe statistics...\n")
    
    # Calculate stats for each sequence/group
    numbOfTotalSeqs <- nrow(baseline@db)
    regions <- baseline@regions
    db <- baseline@db
    if ( "SEQUENCE_ID" %in% colnames(db) ) { db <- subset(db, select="SEQUENCE_ID")}
    list_stats <-
        foreach( i=icount(numbOfTotalSeqs) ) %dopar% {
            df_baseline_seq <- NULL
            db_seq <- db[i,]
            for (region in regions) {
                baseline_pdf <- baseline@pdfs[[region]][i,]
                baseline_ci <- calcBaselineCI(baseline_pdf)
                df_baseline_seq_region <- 
                    data.frame( db[i,],
                                REGION=region,
                                BASELINE_SIGMA=calcBaselineSigma(baseline_pdf),
                                BASELINE_CI_LOWER=baseline_ci[1],
                                BASELINE_CI_UPPER=baseline_ci[2],
                                BASELINE_CI_PVALUE=calcBaselinePvalue(baseline_pdf)
                    )
                df_baseline_seq <- rbind(df_baseline_seq, df_baseline_seq_region)
            }
            return(df_baseline_seq)
        }
    
    # Stop SNOW cluster
    if(nproc > 1) { stopCluster(cluster) }
    
    # Convert list of BASELINe stats into a data.frame
    stats <- do.call(rbind, list_stats)
    
    # Append stats to baseline object
    editBaseline(baselien, field_name = "stats", stats )
    
    return(baseline)   
}


# Given a BASELIne PDF calculate mean sigma
calcBaselineSigma <- function ( baseline_pdf,
                                     max_sigma=20,
                                     length_sigma=4001 ) {
    
    if ( length(baseline_pdf)!=length_sigma) { return(NA) }
    
    sigma_s <- seq(-max_sigma, max_sigma, length.out=length_sigma)
    norm = {length_sigma-1}/2/max_sigma
    return( (baseline_pdf%*%sigma_s/norm)  )
}


# Given a BASELIne PDF calculate Confidence Interval
calcBaselineCI <- function ( baseline_pdf,
                                  low=0.025,
                                  up=0.975,
                                  max_sigma=20,
                                  length_sigma=4001 ){
    
    if ( length(baseline_pdf)!=length_sigma ) { return( c(NA,NA) ) }
    
    sigma_s <- seq(-max_sigma, max_sigma, length.out=length_sigma)
    cdf <- cumsum(baseline_pdf)
    cdf <- cdf/cdf[length(cdf)]
    intervalLow <- findInterval(low,cdf)
    fractionLow <- (low - cdf[intervalLow])/(cdf[intervalLow+1]-cdf[intervalLow])
    intervalUp <- findInterval(up,cdf)
    fractionUp <- (up - cdf[intervalUp])/(cdf[intervalUp]-cdf[intervalUp-1])
    sigmaLow <- sigma_s[intervalLow]+fractionLow*(sigma_s[intervalLow+1]-sigma_s[intervalLow])
    sigmaUp <- sigma_s[intervalUp]+fractionUp*(sigma_s[intervalUp+1]-sigma_s[intervalUp])
    return( c(sigmaLow,sigmaUp) )
}

# Given a BASELIne PDF calculate P value
calcBaselinePvalue <- function ( baseline_pdf, 
                                      length_sigma=4001, 
                                      max_sigma=20 ){
    if ( length(baseline_pdf)>1 ) {
        norm = {length_sigma-1}/2/max_sigma
        pVal = {sum(baseline_pdf[1:{{length_sigma-1}/2}]) + baseline_pdf[{{length_sigma+1}/2}]/2}/norm
        if(pVal>0.5){
            pVal = pVal-1
        }
        return(pVal)
    }else{
        return(NA)
    }
}


# # Given a BASELIne PDF calculate Mean, Confidence Interval (lower & upper) and P value
# calcBaselineStats <- function ( baseline_pdf,
#                                      low=0.025,
#                                      up=0.975,
#                                      max_sigma=20, 
#                                      length_sigma=4001 ){
#     
#     # if NA (i.e. length of baseline_pdf is 1)
#     if ( length(baseline_pdf)==1  ) { return(rep(NA,4)) }
#     
#     
#     baselineSigma <- calculateBaselineSigma( baseline_pdf=baseline_pdf, 
#                                              max_sigma=max_sigma, 
#                                              length_sigma=length_sigma )
#     
#     
#     baselineCI <- calculateBaselineCI( baseline_pdf=baseline_pdf,
#                                        low=low,
#                                        up=up,
#                                        max_sigma=max_sigma,
#                                        length_sigma=length_sigma )
#     
#     baselinePvalue <- calculateBaselinePvalue( baseline_pdf=baseline_pdf,
#                                                max_sigma=max_sigma,
#                                                length_sigma=length_sigma )
#     
#     return( c( "Sigma"=baselineSigma, 
#                "CI_Lower"=baselineCI[1], 
#                "CI_Upper"=baselineCI[2],
#                "Pvalue"=baselinePvalue 
#     ) 
#     )
# }


#### Plotting functions ####

#' Plots the results of BASELINe analysis
#' 
#' \code{plotBaseline} plots the results of selection analysis using the BASELINe method.
#'
#' @param    baseline     \code{Baseline} object containing selection probability 
#'                        density functions.
#' @param    idColumn     name of the column in the \code{data} slot of \code{baseline} 
#'                        containing primary identifiers.
#' @param    groupColumn  name of the column in the \code{data} slot of \code{baseline} 
#'                        containing secondary grouping identifiers. If \code{NULL}, 
#'                        organize the plot only on values in \code{idColumn}.
#' @param    regions      character vector naming the regions to plot, correspoding to the
#'                        regions for which the \code{baseline} data was calculated.
#' @param    style        type of plot to draw. One of:
#'                        \itemize{
#'                          \item \code{"point"}:  plots the mean and confidence interval of
#'                                                 selection scores for each value in 
#'                                                 \code{idColumn}, grouped and colored
#'                                                 by values in \code{groupColumn}.
#'                          \item \code{"curve"}:  plots a set of curves for each probability 
#'                                                 density function in \code{baseline}, 
#'                                                 colored by values in \code{idColumn}.
#'                        }
#' @param    size         numeric scaling factor for lines, points and text in the plot.
#' @param    silent       if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                        object; if \code{FALSE} draw the plot.
#' @param    ...          additional arguments to pass to ggplot2::theme.
#' 
#' @return   A ggplot object defining the plot.
#'
#' @seealso  Takes as input a \code{\link{Baseline}} object.
#'           See \code{\link{getBaseline}} for generating selection probability density 
#'           functions.
#' 
#' @examples
#' library(alakazam)
#' db_file <- system.file("extdata", "Influenza_IB.tab", package="shm")
#' db <- readChangeoDb(db_file)
#' 
#' # Plot mean and confidence interval
#' plotBaseline(baseline, "BARCODE", style="point")
#' 
#' # Plot density function
#' plotBaseline(baseline, "BARCODE", style="curve")
#' 
#' @export
plotBaseline <- function(baseline, idColumn, groupColumn=NULL, regions=c("CDR", "FWR"), 
                         style=c("point", "curve"), size=1, silent=FALSE, ...) {
    # TODO:  add group/subgroup color specification.
    # TODO:  add subgroup plotting
    # TODO:  curve is not right
    
    #idColumn="BARCODE"
    #groupColumn="BARCODE"
    #regions=c("CDR", "FWR")
    #style="point"
    #size=1
    #silent=FALSE
    
    # Check input
    style <- match.arg(style)
    
    # Set base plot settings
    base_theme <- theme_bw() +
        theme(panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        theme(strip.background=element_rect(fill="white"))
    #theme(axis.title.x=element_blank(),
    #      axis.text.x=element_blank(), 
    #      axis.ticks.x=element_blank()) +
    #theme(legend.position="top")
    
    # Define universal plot settings
    base_theme <- theme_bw()
    
    # Plot BASELINe summary
    if (style == "point") {
        # Get summary statistics
        #baseline_df <- getBASELINeStats(baseline[[1]], baseline[[2]], nproc=1)
        baseline_df <- as.data.frame(baseline[[1]], row.names=1:nrow(baseline[[1]]), stringsAsFactors=FALSE)
        baseline_df <- transform(baseline_df, SIGMA=as.numeric(SIGMA), CI_LOWER=as.numeric(CI_LOWER), 
                                 CI_UPPER=as.numeric(CI_UPPER), P_VALUE=as.numeric(P_VALUE))
        
        # Subset to regions of interest
        sub_df <- subset(baseline_df, REGION %in% regions)
        
        # Build plot
        p1 <- ggplot(sub_df, aes_string(x=idColumn, y="SIGMA", ymax=max("SIGMA"))) +
            base_theme + 
            #theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
            #ggtitle(region) +
            xlab("") +
            ylab(expression(Sigma)) +
            #scale_shape_manual(name="Isotype", values=isotype_shapes) +
            #scale_color_manual(name="Sample Type", values=CLASS_COLORS) +
            #scale_fill_manual(name="Sample Type", values=CLASS_COLORS) +
            geom_hline(yintercept=0, size=1*size, linestyle=2, color="grey") +
            geom_point(size=3*size, position=position_dodge(0.4)) +
            geom_errorbar(aes(ymin=CI_LOWER, ymax=CI_UPPER), 
                          width=0.2, size=0.5*size, alpha=0.8, position=position_dodge(0.4))
        if (!is.null(groupColumn)) {
            p1 <- p1 + aes_string(color=groupColumn)
        }
        p1 <- p1 + facet_grid(REGION ~ .)
    } else if (style == "curve") {
        # Mess with input to get curves
        id_names <- unique(baseline[[1]][, idColumn])
        dens_names <- names(baseline[[2]][[1]]@probDensity)
        
        # Subset to regions of interest
        dens_names <- dens_names[dens_names %in% regions]
        dens_list <- list()
        for (n in dens_names) {
            # Extract density into matrix
            tmp_mat <- t(sapply(baseline[[2]], function(x) x@probDensity[[n]]))
            colnames(tmp_mat) <- 1:ncol(tmp_mat)
            rownames(tmp_mat) <- id_names
            dens_list[[n]] <- tmp_mat
        }
        
        # Melt density matrices
        melt_list <- list()
        for (n in dens_names) {
            tmp_melt <- reshape2::melt(dens_list[[n]], varnames=c(idColumn, "x"), value.name="SIGMA")
            melt_list[[n]] <- tmp_melt
        }
        dens_df <- ldply(melt_list, .id="REGION")
        
        # Build plot
        p1 <- ggplot(dens_df, aes_string(x="x", y="SIGMA")) +
            base_theme + 
            xlab("x") +
            ylab(expression(Sigma)) +
            geom_line(aes_string(color=idColumn), size=1*size)
        if (!is.null(groupColumn)) {
            p1 <- p1 + aes_string(linetype=groupColumn)
        }
        p1 <- p1 + facet_grid(REGION ~ .)
    }
    
    # Add additional theme elements
    p1 <- p1 + do.call(theme, list(...))
    
    # Plot
    if (!silent) { 
        plot(p1)
    }
    
    invisible(p1)
}


#### Original BASELINe functions ####

##Covolution
break2chunks<-function(G=1000){
    base<-2^round(log(sqrt(G),2),0)
    return(c(rep(base,floor(G/base)-1),base+G-(floor(G/base)*base)))
}

PowersOfTwo <- function(G=100){
    exponents <- array()
    i = 0
    while(G > 0){
        i=i+1
        exponents[i] <- floor( log2(G) )
        G <- G-2^exponents[i]
    }
    return(exponents)
}

convolutionPowersOfTwo <- function( cons, length_sigma=4001 ){
    G = ncol(cons)
    if(G>1){
        for(gen in log(G,2):1){
            ll<-seq(from=2,to=2^gen,by=2)
            sapply(ll,function(l){cons[,l/2]<<-weighted_conv(cons[,l],cons[,l-1],length_sigma=length_sigma)})
        }
    }
    return( cons[,1] )
}

convolutionPowersOfTwoByTwos <- function( cons, length_sigma=4001,G=1 ){
    if(length(ncol(cons))) G<-ncol(cons)
    groups <- PowersOfTwo(G)
    matG <- matrix(NA, ncol=length(groups), nrow=length(cons)/G )
    startIndex = 1
    for( i in 1:length(groups) ){
        stopIndex <- 2^groups[i] + startIndex - 1
        if(stopIndex!=startIndex){
            matG[,i] <- convolutionPowersOfTwo( cons[,startIndex:stopIndex], length_sigma=length_sigma )
            startIndex = stopIndex + 1
        }
        else {
            if(G>1) matG[,i] <- cons[,startIndex:stopIndex]
            else matG[,i] <- cons
            #startIndex = stopIndex + 1
        }
    }
    return( list( matG, groups ) )
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

# Given a list of PDFs, returns a convoluted PDF
groupPosteriors <- function( listPosteriors, max_sigma=20, length_sigma=4001 ,Threshold=2 ){
    listPosteriors = listPosteriors[ !is.na(listPosteriors) ]
    Length_Postrior<-length(listPosteriors)
    if(Length_Postrior>1 & Length_Postrior<=Threshold){
        cons = matrix(unlist(listPosteriors),length(listPosteriors[[1]]),length(listPosteriors))
        listMatG <- convolutionPowersOfTwoByTwos(cons,length_sigma=length_sigma)
        y<-calculate_bayesGHelper(listMatG,length_sigma=length_sigma)
        return( y/sum(y)/(2*max_sigma/(length_sigma-1)) )
    }else if(Length_Postrior==1) return(listPosteriors[[1]])
    else  if(Length_Postrior==0) return(NA)
    else {
        cons = matrix(unlist(listPosteriors),length(listPosteriors[[1]]),length(listPosteriors))
        y = fastConv(cons,max_sigma=max_sigma, length_sigma=length_sigma )
        return( y/sum(y)/(2*max_sigma/(length_sigma-1)) )
    }
}

fastConv<-function(cons, max_sigma=20, length_sigma=4001){
    chunks<-break2chunks(G=ncol(cons))
    if(ncol(cons)==3) chunks<-2:1
    index_chunks_end <- cumsum(chunks)
    index_chunks_start <- c(1,index_chunks_end[-length(index_chunks_end)]+1)
    index_chunks <- cbind(index_chunks_start,index_chunks_end)
    
    case <- sum(chunks!=chunks[1])
    if(case==1) End <- max(1,((length(index_chunks)/2)-1))
    else End <- max(1,((length(index_chunks)/2)))
    
    firsts <- sapply(1:End,function(i){
        indexes<-index_chunks[i,1]:index_chunks[i,2]
        convolutionPowersOfTwoByTwos(cons[ ,indexes])[[1]]
    })
    if(case==0){
        result<-calculate_bayesGHelper( convolutionPowersOfTwoByTwos(firsts) )
    }else if(case==1){
        last<-list(calculate_bayesGHelper(
            convolutionPowersOfTwoByTwos( cons[ ,index_chunks[length(index_chunks)/2,1]:index_chunks[length(index_chunks)/2,2]] )
        ),0)
        result_first<-calculate_bayesGHelper(convolutionPowersOfTwoByTwos(firsts))
        result<-calculate_bayesGHelper(
            list(
                cbind(
                    result_first,last[[1]]),
                c(log(index_chunks_end[length(index_chunks)/2-1],2),log(index_chunks[length(index_chunks)/2,2]-index_chunks[length(index_chunks)/2,1]+1,2))
            )
        )
    }
    return(as.vector(result))
}
