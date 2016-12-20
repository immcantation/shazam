Version 0.1.4.999: December 20, 2016
-------------------------------------------------------------------------------
Targeting Models:

+ Added `numMutationsOnly` argument to `createSubstitutionMatrix()`, enabling
  parameter tuning for `minNumMutations`.  
+ Added functions `minNumMutationsTune()` and `minNumSeqMutationsTune()` to 
  tune for parameters `minNumMutations` and `minNumSeqMutations` in functions 
  `createSubstitutionMatrix()` and `createMutabilityMatrix()` respectively. 
  Also added function `plotTune()` which helps visualize parameter tuning using
  the abovementioned two new functions. 
+ Added human kappa and lambda light chain, silent, 5-mer, functional targeting
  model (`HKL_S5F`).
+ Renamed `HS5FModel` as `HH_S5F`, `MRS5NFModel` as `MK_RS5NF`, and `U5NModel` 
  as `U5N`.
+ Added human heavy chain, silent, 1-mer, functional substitution model (`HH_S1F`),
  human kappa and lambda light chain, silent, 1-mer, functional substitution model 
  (`HKL_S1F`), and mouse kappa light chain, replacement and silent, 1-mer, 
  non-functional substitution model (`MK_RS1NF`).
+ Added `makeDegenerate5merSub` and `makeDegenerate5merMut` which make degenerate
  5-mer substitution and mutability models respectively based on the 1-mer models. 
  Also added `makeAverage1merSub` and `makeAverage1merMut` which make 1-mer 
  substitution and mutability models respectively by averaging over the 5-mer models. 

Mutation Profiling:

+ Added `returnRaw` argument to `calcObservedMutations()`, which if true returns 
  the positions of point mutations and their corresponding mutation types, as 
  opposed to counts of mutations (hence "raw"). 
+ Added new functions `slideWindowSeq()` and `slideWindowDb()` which implement 
  a sliding window approach towards filtering a single sequence or sequences in
  a data.frame which contain(s) equal to or more than a given number of mutations 
  in a given number of consecutive nucleotides.
+ Added new function `slideWindowTune()` which allows for parameter tuning for
  using `slideWindowSeq()` and `slideWindowDb()`.
+ Added new function `slideWindowTunePlot()` which visualizes parameter tuning 
  by `slideWindowTune()`.
  
Distance Calculation:

+ Added `findThreshold` function to infer clonal distance threshold from 
  nearest neighbor distances returned by `distToNearest`.
+ Removed `HS1FDistance` (replaced by `HH_S1F`) and `M1NDistance` (replaced by 
  `MK_RS1NF`). 
+ Updated `calcTargetingDistance()` to enable calculation of a symmetric distance
  matrix given a 1-mer substitution matrix normalized by row, such as `HH_S1F`.  

Region Definition:

+ Deleted `IMGT_V_NO_CDR3` and `IMGT_V_BY_REGIONS_NO_CDR3`. Updated `IMGT_V` 
  and `IMGT_V_BY_REGIONS` so that neither includes CDR3 now.


Version 0.1.4:  August 5, 2016
-------------------------------------------------------------------------------

Selection Analysis:

+ Fixed a bug in calcBaseline wherein the germline column was incorrected 
  hardcoded, leading to erroneous mutation counts for some clonal consensus 
  sequences.

Targeting Models:

+ Added `numSeqMutationsOnly` argument to `createMutabilityMatrix()`, enabling
  parameter tuning for `minNumSeqMutations`.
  

Version 0.1.3:  July 31, 2016
-------------------------------------------------------------------------------

General:

+ Added ape and igraph dependency
+ Removed the `InfluenzaDb` data object, in favor of the updated `ExampleDb`
  provided in alakazam 0.2.4.
+ Added conversion of sequence to uppercase for several functions to support
  data that was not generated via Change-O.

Distance Calculation:

+ Added the `cross` argument to `distToNearest()` which allows restriction of 
  distances to only distances across samples (ie, excludes within-sample 
  distances).
+ Added `mst` flag to `distToNearest()`, which will return all distances to 
  neighboring nodes in a minimum spanning tree.
+ Updated single nucleotide distance models to use the new C++ distance
  methods in alakazam 0.2.4 for better performance.
+ Fixed a bug leading to failed distance calculations for the `aa` model 
  of `distToNearest()`.
+ Fixed a bug wherein gap characters where being translated into Ns (Asn) 
  rather than Xs within the `aa` model of `distToNearest()`.

Mutation Profiling:

+ Added the `MutationDefinition` `VOLUME_MUTATIONS`.
+ Added the functions `shmulateSeq()` and `shmulateTree()` to simulate
  mutations on sequences and lineage trees, respectively, using a 5-mer
  targeting model.
+ Renamed `collapseByClone`, `calcDbExpectedMutations` and 
  `calcDbObservedMutations` to `collapseClones`, `expectedMutations`,
  and `observedMutations`, respectively.
  
Selection Analysis:

+ Fixed a bug wherein passing a `Baseline` object through `groupBaseline()`
  multiple times resulted in incorrect normalization.
+ Added `title` options to `plotBaselineSummary()` and `plotBaselineDensity()`.
+ Added more control over colors and group ordering to `plotBaselineSummary()` 
  and `plotBaselineDensity()`.
+ Added the `testBaseline()` function to test the significance of 
  differences between two selection distributions.
+ Improved selection analysis vignette. 


Version 0.1.2:  February 20, 2016
-------------------------------------------------------------------------------

General:

+ Renamed package from shm to shazam.
+ Internal changes to conform to CRAN policies.
+ Compressed and moved example database to the data object `InfluenzaDb`.
+ Fixed several bugs where functions would not work properly when passed 
  a `dplyr::tbl_df` object instead of a `data.frame`.
+ Changed R dependency to R >= 3.1.2.
+ Added stringi dependency.

Distance Calculation:

+ Fixed a bug wherein `distToNearest()` did not return the nearest neighbor 
  with a non-zero distance.

Targeting Models:

+ Performance improvements to `createSubstitutionMatrix()`,  
  `createMutabilityMatrix()`, and `plotMutability()`.
+ Modified color scheme in `plotMutability()`.
+ Fixed errors in the targeting models vignette.

Mutation Profiling:

+ Added the `MutationDefinition` objects `MUTATIONS_CHARGE`, 
  `MUTATIONS_HYDROPATHY`, `MUTATIONS_POLARITY` providing alternate approaches
  to defining replacement and silent annotations to mutations when calling
  `calcDBObservedMutations()` and `calcDBExpectedMutations()`.
+ Fixed a few bugs where column names, region definitions or mutation models
  were not being recognized properly when non-default values were used.
+ Made the behavior of `regionDefinition=NULL` consistent for all mutation
  profiling functions. Now the entire sequence is used as the region and
  calculations are made accordingly.
+ `calcDBObservedMutations()` returns R and S mutations also 
  when `regionDefinition=NULL`. Older versions reported the sum of R and S 
  mutations. The function will add the columns `OBSERVED_SEQ_R` and
  `OBSERVED_SEQ_S` when `frequency=FALSE`, and `MU_FREQ_SEQ_R` and 
  `MU_FREQ_SEQ_R` when `frequency=TRUE`.
    

Version 0.1.1:  December 18, 2015
-------------------------------------------------------------------------------

General:

+ Swapped dependency on doSNOW for doParallel.
+ Swapped dependency on plyr for dplyr.
+ Swapped dependency on reshape2 for tidyr.
+ Documentation clean up.

Distance Calculation:

+ Changed underlying method of calcTargetingDistance to be negative log10 of
  the probability that is then centered at one by dividing by the mean 
  distance.
+ Added `symmetry` parameter to distToNearest to change behavior of how 
  asymmetric distances (A->B != B->A) are combined to get distance 
  between A and B. 
+ Updated error handling in distToNearest to issue warning when unrecognized
  character is in the sequence and return an NA.
+ Fixed bug in 'aa' model in distToNearest that was calculating distance
  incorrectly when normalizing by length.
+ Changed behavior to return nearest nonzero distance neighbor.

Mutation Profiling:

+ Renamed calcDBClonalConsensus to collapseByClone
  Also, renamed argument collapseByClone to expandedDb.
+ Fixed a (major) bug in calcExpectedMutations. Previously, the targeting calculation
  was incorrect and resulted in incorrect expected mutation frequencies. Note, that this
  also resulted in incorrect BASELINe Selection (Sigma) values.
+ Changed denominator in calcObservedMutations to be based on informative 
  (unambiguous) positions only.
+ Added nonTerminalOnly parameter to calcDBClonalConsensus indicating whether
  to consider mutations at leaves or not (defaults to false).
  
Selection Analysis:

+ Updated groupBaseline. Now when regrouping a Baseline object (i.e. grouping previously
  grouped PDFs) weighted convolution is performed. 
+ Added "imbalance" test statistic to the Baseline selection calculation.
+ Extended the Baseline Object to include binomK, binomN and binomP
  Similar to numbOfSeqs, each of these are a matrix. They contain binomial inputs for 
  each sequence and region. 
  
Targeting Models:

+ Added `minNumMutations` parameter to createSubstitutionMatrix. This is the 
  minimum number of observed 5-mers required for the substituion model. 
  The substitution rate of 5-mers with fewer number of observed mutations
  will be inferred from other 5-mers. 
+ Added `minNumSeqMutations` parameter to createMutabilityMatrix. This is the 
  minimum number of mutations required in sequences containing the 5-mers of 
  interest. The mutability of 5-mers with fewer number of observed mutations 
  in the sequences will be inferred. 
+ Added `returnModel` parameter to createSubstitutionMatrix. This gives user 
  the option to return 1-mer or 5-mer model.
+ Added `returnSource` parameter to createMutabilityMatrix. If TRUE, the 
  code will return a data frame indicating whether each 5-mer mutability is 
  observed or inferred. 
+ In createSubstitutionMatrix and createMutabilityMatrix, fixed a bug when
  multipleMutation is set to "ignore".
+ Changed inference procedure for the 5-mer substitution model.
+ Added inference procedure for 5-mers without enough observed mutations
  in the mutability model.
+ Fixed a bug in background 5-mer count for the RS model.
+ Fixed a bug in IMGT gap handling in createMutabilityMatrix.
+ Fixed a bug that occurs when sequences are in lower cases.


Version 0.1.0:  June 18, 2015
-------------------------------------------------------------------------------

Initial public release.

General:

+ Restructured the S4 class documentation.
+ Fixed bug wherein example `Influenza.tab` file did not load on Mac OS X.
+ Added citations for `citation("shazam")` command.
+ Added dependency on data.table >= 1.9.4 to fix bug that occured with 
  earlier versions of data.table.

Distance Calculation:

+ Added a human 1-mer substitution matrix, `HS1FDistance`, based on the
  Yaari et al, 2013 data.
+ Set the `hs1f` as the default distance model for `distToNearest()`.
+ Added conversion of sequences to uppercase in `distToNearest()`.
+ Fixed a bug wherein unrecongized (including lowercase) characters would
  lead to silenting returning a distance of 0 to the neared neighbor. 
  Unrecognized characters will now raise an error.

Mutation Profiling:

+ Fixed bug in `calcDBClonalConsensus()` so that the function now works 
  correctly when called with the argument `collapseByClone=FALSE`.
+ Added the `frequency` argument to `calcObservedMutations()` and
  `calcDBObservedMutations()`, which enables return of mutation frequencies
  rather the default of mutation counts.
  
Targeting Models:

+ Removed `M3NModel` and all options for using said model.
+ Fixed bug in `createSubstitutionMatrix()` and `createMutabilityMatrix()` 
  where IMGT gaps were not being handled.


Version 0.1.0.beta-2015-05-30:  May 30, 2015
-------------------------------------------------------------------------------

General:

+ Added more error checking.

Targeting Models:

+ Updated the targeting model workflow to include a clonal consensus step.


Version 0.1.0.beta-2015-05-11:  May 11, 2015
-------------------------------------------------------------------------------

Targeting Models:

+ Added the `U5NModel`, which is a uniform 5-mer model.
+ Improvements to `plotMutability()` output.


Version 0.1.0.beta-2015-05-05:  May 05, 2015
-------------------------------------------------------------------------------

Prerelease for review.
