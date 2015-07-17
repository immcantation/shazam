Version 0.1.1.beta:  July 17, 2015
-------------------------------------------------------------------------------

Distance Profiling:

+ Changed underlying method of calcTargetingDistance to be negative log10 of
  the probability that is then centered at one by dividing by the mean 
  distance.

Mutation Profiling:

+ Changed denominator in calcObservedMutations to be based on informative 
  (unambiguous) positions only.
+ Added nonTerminalOnly parameter to calcDBClonalConsensus indicating whether
  to consider mutations at leaves or not (defaults to false).


Version 0.1.0:  June 18, 2015
-------------------------------------------------------------------------------

Initial public release.

General:

+ Restructured the S4 class documentation.
+ Fixed bug wherein example `Influenza.tab` file did not load on Mac OS X.
+ Added citations for `citation("shm")` command.
+ Added dependency on data.table >= 1.9.4 to fix bug that occured with 
  earlier versions of data.table.

Distance Profiling:

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