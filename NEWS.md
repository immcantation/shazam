Version 1.1.2.999: February 03, 2023
-------------------------------------------------------------------------------

General:

+ Updated dependencies to ggplot2 >= 3.4.0.

Version 1.1.2: September 26, 2022
-------------------------------------------------------------------------------

Mutation Profiling:

+ Bug fix in parallelization set up for functions `slideWindowTune` 
  and `slideWindowDb`.
  
+ `plotSlideWindowTune` (`slideWindowTunePlot`). Updated the possible
   values of the parameter `plotFiltered`, for easier usage. The new values 
   (and their equivalent values in `slideWindowTunePlot`) are `filtered` (`TRUE`), 
   `remaining` (`FALSE`), and `per_mutation` (`NULL`). 
  
Deprecated:

  + Deprecated `slideWindowTunePlot` in favor of `plotSlideWindowTune`, for naming
  consistency.

Version 1.1.1: May 23, 2022
-------------------------------------------------------------------------------

General:

+ Removed dependency: kedd. The CRAN kedd package (by Arsalane Chouaib Guidoum) 
  has been scheduled for archival on 2022-05-25. We have adapted the functions
  used by shazam and removed the dependency.
  
New feature:

+ Added the function `convertNumbering` to convert between numbering systems (IMGT, Kabat).

Mutation Profiling:

+ `shmulateTree` has new argument `nproc` to specify the number of cores. Default
   values `mutThresh` and `windowSize` have been set to `mutThresh=6` and
   `windowSize=10`.
   
+ Added the option `plotFiltered=NULL` to `slideWindowTunePlot`.
   
+ Fixed a bug in `listObservedMutations` not returning a list when `db` had
  one sequence with one mutation.
  
+ Fixed bars shifted in `plotMutability`.


Version 1.1.0: July 8, 2021
-------------------------------------------------------------------------------

General:

+ Updated dependencies to alakazam >= 1.1.0 and ggplot2 >= 3.3.4.

Selection Analysis:

+ `observedMutations`, `expectedMutations`, and `calcBaseline` can analyze 
  mutations in all regions (CDR1, CDR2, CDR3, FWR1, FWR2, FWR3 and FWR4) by 
  specifying `regionDefinition=IMGT_VDJ` or 
  `regionDefinition=IMGT_VDJ_BY_REGIONS`.
+ Added the function `setRegionBoundaries` to build sequence-specific 
  `RegionDefinition` objects extending to CDR3 and FWR4.
+ Added the function `makeGraphDf` to facilitate mutational analysis on 
  lineage trees.

Distance Profiling:

+ Fixed a bug in `distToNearest` where TRB and TRD sequences where ignored in 
  distance calculation.
+ Fixed a bug in `distToNearest` causing a fatal error when `cross` was set.
+ Fixed a bug in `nearestDist` causing a fatal error when using `model="aa"` 
  and `crossGroups`.

Targeting Models:

+ Fixed an incompatibility with newer versions of ggplot2 in `plotMutability`.


Version 1.0.2: August 10, 2020
-------------------------------------------------------------------------------

Mutation Profiling:

+ Fixed a bug in `observedMutations` and `calcObservedMutations` causing 
  mutation counting to fail when there are gap (`-`) characters in the 
  germline sequence.

Targeting Models:

+ Fixed a bug in `createTargetingModel` causing empty counts in the 
  `numMutS` and `numMutR` slots.


Version 1.0.1: July 18, 2020
-------------------------------------------------------------------------------

Distance Profiling:

+ Added support for TCR genes to `distToNearest`. 
+ Renamed the `groupUsingOnlyIGH` argument of `distToNearest` to `onlyHeavy`.


Version 1.0.0 May 9, 2020
-------------------------------------------------------------------------------

Backwards Incompatible Changes:

+ Changed default expected data format from the Change-O data format to the
  AIRR Rearrangement standard. For example: where functions used the column 
  name `V_CALL` (Change-O) as the default to identify the field that stored 
  the V gene calls, they now use `v_call` (AIRR). That means, scripts that 
  relied on default values (previously, `v_call="V_CALL"`), will now fail if 
  calls to the functions are not updated to reflect the correct value for the 
  data. If data are in the Change-O format, the current default value 
  `v_call="v_call"` will fail to identify the column with the V gene calls
  as the column `v_call` doesn't exist. In this case, `v_call="V_CALL"` needs 
  to be specified in the function call.
+ `ExampleDb` converted to the AIRR Rearrangement standard and examples 
  updated accordingly.
+ For consistency with the style of the new data format default, other field 
  names have been updated to use the same capitalization. This change affects:
    - Region definitions. For example, the `labels` slot of `IMGT_V` has 
      changed from `CDR_R`, `CDR_S`, `FWR_R` and `FWR_S` to `cdr_r`, `cdr_s`, 
      `fwr_r` and `fwr_s`, respectively.
    - Mutations in `CODON_TABLE` and the different `MUTATION_SCHEMES` change
      from `R`, `S` and `Stop` to `r`, `s` and `stop`, respectively. 
    - Mutation profiling function output columns. For example, from 
      `MU_COUNT_SEQ` to `mu_count_seq`.
    - `calcBaseline` and related function output columns and S4 object slots. 
      For example, from `PVALUE`, `REGION` and `BASELINE_CI_PVALUE` to 
      `pvalue`, `region` and `baseline_ci_pvalue`, respectively.
 + Model names used by `createSubstitutionMatrix`, `createMutabilityMatrix` and 
   `createTargetingModel`, changed from `model=c("S","RS")` to 
   `model=c("s","rs")`.

General:

+ License changed to AGPL-3.

Targeting Models:

+ `createMutabilityMatrix`, `extendMutabilityMatrix`, `createTargetingMatrix`,
  and `createTargetingModel` now also returns the numbers of silent and
  replacement mutations used for estimating the 5-mer mutabilities. These 
  numbers are recorded in the `numMutS` and `numMutR` slots in the newly
  defined `MutabilityModel`, `MutabilityModelWithSource`, and `TargetingMatrix` 
  classes.
  
Mutation Profiling:

+ `shmulateSeq` now also supports specifying the frequency of mutations to be 
   introduced. (Previously, only the number of mutations was supported.)


Version 0.2.3 February 5, 2020
-------------------------------------------------------------------------------

General:

+ Removed SDMTools dependency.


Version 0.2.2 December 15, 2019
-------------------------------------------------------------------------------

General:

+ Fixed an incompatibility with R 4.0 matrix changes.

Version 0.2.1 July 19, 2019
-------------------------------------------------------------------------------

Distance Calculation:

+ Fixed a bug in `distToNearest` that could potentially cause sequences from
  different partitions to be used for distance calculation. 


Version 0.2.0 July 18, 2019
-------------------------------------------------------------------------------

General:

+ Upgraded to alakazam >= 0.3.0 and dplyr >= 0.8.1.

Distance Calculation:

+ Fixed a bug in `plotDensityThreshold` for negative densities.
+ Fixed a bug in `distToNearest` for performing subsampling while calculating
  cross-group nearest neighbor distances.
+ For partitioning sequences, `distToNearest` now supports, via a new argument
  `VJthenLen`, either a 2-stage partitioning (first by V gene and J gene, then 
  by junction length), or a 1-stage partitioning (simultaneously by V gene, J 
  gene, and junction length). For 1-stage partitioning, `distToNearest` supports
  export of the partitioning information as a new column via `keepVJLgroup`.
+ `distToNearest` now supports single-cell input data with the addition of new
  arguments `cellIdColumn`, `locusColumn`, and `groupUsingOnlyIGH`.

Mutation Profiling:

+ `shmulateTree` has new arguments, `start` and `end`, to specify the region
  in the sequence where mutations can be introduced. 

Selection Analysis:

+ Added the function `consensusSequence` which can be used to build a 
  consensus sequence using a variety of methods.


Version 0.1.11: January 27, 2019
-------------------------------------------------------------------------------

General:

+ Fixed a bug in the prototype declarations for the `TargetingModel` and 
  `RegionDefinition` S4 classes.


Version 0.1.10: September 19, 2018
-------------------------------------------------------------------------------

General:

+ Added `subsample` argument to `distToNearest` function.
+ Removed some internal utility functions in favor of importing them from 
  `alakazam`. Specifically, `progressBar`, `getBaseTheme` and `checkColumns`.
+ Removed `clearConsole`, `getnproc`, and `getPlatform` functions.

Distance Calculation:

+ Changed default `findThreshold` method to `density`.
+ Significantly reduced run time of the `density` method by retuning the 
  bandwidth detection process. The `density` method should now also yield more 
  consistent thresholds, on average.
+ The `subsample` argument to `findThreshold` now applies to both the 
  `density` and `gmm` methods. Subsampling of distance is not performed by 
  default.
+ Fixed a bug in `plotDensityThreshold` and `plotGmmThreshold` wherein the
  `breaks` argument was ignored when specifying `xmax` and/or `xmin`.

Selection Analysis:

+ Fixed a plotting bug in `plotBaselineDensity` arising when the `groupColumn` 
  and `idColumn` arguments were set to the same column.
+ Added the `sizeElement` argument to `plotBaselineDensity` to control 
  line size
+ Renamed the `field_name` argument to `field` in `editBaseline`.


Version 0.1.9: March 30, 2018
-------------------------------------------------------------------------------

Selection Analysis:

+ Fixed a bug in `plotBaselineDensity` which caused an empty plot to be
  generated if there was only a single value in the `idColumn`.
+ Fixed a bug in `calcBaseline` which caused a crash in `summarizeBaseline` 
  and `groupBaseline` when input `baseline` is based on only 1 sequence 
  (i.e. when `nrow(baseline@db)` is 1).
+ Set default `plot` call on a `Baseline` object to `plotBaselineDensity`.
+ Removed `getBaselineStats` function.
+ Added a `summary` method for `Baseline` objects that calls 
  `summarizeBaseline` and returns a data.frame.

Mutation Profiling:

+ Fixed a bug in `shmulateSeq` which caused a crash when the input 
  sequence contains gaps (`.`).
+ Renamed the argument `mutations` in `shmulateSeq` to `numMutations`.
+ Improved help documentation for `shmulateSeq` and `shmulateTree`.
+ Added vignette for simulating mutated sequences.
+ `calcExpectedMutations` will now treat non-ACTG characters as Ns rather
  than produce an error.
+ Added two new `RegionDefinition` objects for the full V segment as
  single region (`IMGT_V_BY_SEGMENTS`) and the V segment with each
  codon as a separate region (`IMGT_V_BY_CODONS`).

Targeting Models:

+ Added the `calculateMutability` function which computes the aggregate 
  mutability for sequences.
+ Fixed a bug that caused `createSubstitutionMatrix` to fail for data 
  containing only a single V family.
+ Changed the default model to silent mutations only (`model="S"`) in 
  `createSubstitutionMatrix`, `createSubstitutionMatrix` and 
  `createTargetingModel`
+ Set default `plot` call on a `TargetingModel` object to `plotMutability`.


Version 0.1.8: June 30, 2017
-------------------------------------------------------------------------------

General:

+ Corrected several functions so that they accept both tibbles and data.frames.

Distance Calculation:

+ Adding new fitting procedures to the `"gmm"` method of `findThreshold()` 
  that allows users to choose a mixture of two univariate density distribution 
  functions among four available combinations: `"norm-norm"`, `"norm-gamma"`,  
  `"gamma-norm"`, or `"gamma-gamma"`.
+ Added the ability to choose the threshold selection criteria in the `"gmm"`
  method of `findThreshold()` from the best average sensitivity and specificity, 
  the curve intersection or user defined sensitivity or specificity.
+ Renamed the `cutEdge` argument of `findThreshold()` to `edge`.

Mutation Profiling:

+ Redesigned `collapseClones()`, adding various deterministic and stochastic
  methods to obtain effective clonal sequences, support for including ambiguous 
  IUPAC characters in output, as well as extensive documentation. Removed 
  `calcClonalConsensus()` from exported functions.
+ Added support for including ambiguous IUPAC characters in input for 
  `observedMutations()` and `calcObservedMutations()`.
+ Fixed a minor bug in calculating the denominator for mutation frequency in 
  `calcObservedMutations()` for sequences with non-triplet overhang at the tail.
+ Renamed column names of observed mutations (previously `OBSERVED`) and 
  expected mutations (previously `EXPECTED`) returned by `observedMutations()`
  and `expectedMutations()` to `MU_COUNT` and `MU_EXPECTED` respectively.

Selection Analysis:

+ `calcBaseline()` no longer calls `collapseClones()` automatically if a `CLONE`
  column is present. As indicated by the documentation for `calcBaseline()` 
  users are advised to obtain effective clonal sequences (for example, calling
  `collapseClones()`) before running `calcBaseline()`.
+ Updated vignette to reflect changes in `calcBaseline()`.


Version 0.1.7: May 14, 2017
-------------------------------------------------------------------------------

Mutation Profiling:

+ Fixed a bug in `collapseClones()` that prevented it from running when `nproc`
  is greater than 1.


Version 0.1.6: May 12, 2017
-------------------------------------------------------------------------------

General:

+ Internal changes for compatibility with dplyr v0.6.0.
+ Removed data.table dependency.

Mutation Profiling:

+ Fixed a bug in `collapseClones()` that resulted in erroneous `CLONAL_SEQUENCE`
  and `CLONAL_GERMLINE` being returned.
+ Added a vignette describing basic mutational analysis.
+ Remove console notification that `observedMutations` was running.


Version 0.1.5: March 23, 2017
-------------------------------------------------------------------------------

General:

+ License changed to Creative Commons Attribution-ShareAlike 4.0 International
  (CC BY-SA 4.0).

Selection Analysis:

+ Fixed a bug in p-value calculation in `summarizeBaseline()`. The returned 
  p-value can now be either positive or negative. Its magnitude (without the 
  sign) should be interpreted as per normal. Its sign indicates the direction 
  of the seLicense chalection detected. A positive p-value indicates positive selection, 
  whereas a negative p-value indicates negative selection.
+ Added `editBaseline()` to exported functions, and a corresponding section 
  in the vignette. 
+ Fixed a bug in counting the total number of observed mutations when performing
  a local test for codon-by-codon selection analysis in `calcBaseline()`.

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

+ Fixed a bug in `distToNearest` wherein `normalize="length"` for 5-mer models
  was resulting in distances normalized by junction length squared instead of
  raw junction length.
+ Fixed a bug in `distToNearest` wherein `symmetry="min"` was calculating the 
  minimum of the total distance between two sequences instead of the minimum
  distance at each mutated position.
+ Added `findThreshold` function to infer clonal distance threshold from 
  nearest neighbor distances returned by `distToNearest`.
+ Renamed the `length` option for the `normalize` argument of `distToNearest`
  to `len` so it matches Change-O.
+ Deprecated the `HS1FDistance` and `M1NDistance` distance models, which have 
  been renamed to `hs1f_compat` and `m1n_compat` in the `model` argument of
  `distToNearest`. These deprecated models should be used for compatibility 
  with DefineClones in Change-O v0.3.3. These models have been replaced by 
  replaced by `hh_s1f` and `mk_rs1nf`, which are supported by Change-O v0.3.4. 
+ Renamed the `hs5f` model in `distToNearest` to `hh_s5f`.
+ Added support for `MK_RS5NF` models to `distToNearest`.
+ Updated `calcTargetingDistance()` to enable calculation of a symmetric distance
  matrix given a 1-mer substitution matrix normalized by row, such as `HH_S1F`.
+ Added a Gaussian mixture model (GMM) approach for threshold determination to 
  `findThreshold`. The previous smoothed density method is available via the 
  `method="density"` argument and the new GMM method is available via
  `method="gmm"`.
+ Added the functions `plotGmmThreshold` and `plotDensityThreshold` to plot 
  the threshold detection results from `findThreshold` for the `"gmm"` and
  `"density"` methods, respectively.

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
