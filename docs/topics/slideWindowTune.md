**slideWindowTune** - *Parameter tuning for sliding window approach*

Description
--------------------

Apply [slideWindowDb](slideWindowDb.md) over a search grid made of combinations of `mutThresh` and 
`windowSize` to help with picking a pair of values for these parameters. Parameter 
tuning can be performed by choosing a combination that gives a reasonable number of 
filtered/remaining sequences.


Usage
--------------------
```
slideWindowTune(
db,
sequenceColumn = "sequence_alignment",
germlineColumn = "germline_alignment_d_mask",
dbMutList = NULL,
mutThreshRange,
windowSizeRange,
verbose = TRUE,
nproc = 1
)
```

Arguments
-------------------

db
:   `data.frame` containing sequence data.

sequenceColumn
:   name of the column containing IMGT-gapped sample sequences.

germlineColumn
:   name of the column containing IMGT-gapped germline sequences.

dbMutList
:   if supplied, this should be a list consisting of `data.frame`s 
returned as `$pos` in the nested list produced by 
[calcObservedMutations](calcObservedMutations.md) with `returnRaw=TRUE`; otherwise, 
[calcObservedMutations](calcObservedMutations.md) is called on columns `sequenceColumn`
and `germlineColumn` of `db`. Default is `NULL`.

mutThreshRange
:   range of threshold on the number of mutations in `windowSize` 
consecutive nucleotides to try. Must be between 1 and 
maximum `windowSizeRange` inclusive.

windowSizeRange
:   range of length of consecutive nucleotides to try. The lower end
must be at least 2.

verbose
:   whether to print out messages indicating current progress. Default
is `TRUE`.

nproc
:   Number of cores to distribute the operation over. If the 
`cluster` has already been set earlier, then pass the 
`cluster`. This will ensure that it is not reset.




Value
-------------------

a list of logical matrices. Each matrix corresponds to a `windowSize` in 
`windowSizeRange`. Each column in a matrix corresponds to a `mutThresh` in
`mutThreshRange`. Each row corresponds to a sequence. `TRUE` values
mean the sequences has at least the number of mutations specified in the column name,
for that `windowSize`.


Details
-------------------

If, in a given combination of `mutThresh` and `windowSize`, `mutThresh` 
is greater than `windowSize`, `NA`s will be returned for that particular
combination. A message indicating that the combination has been "skipped" will be 
printed if `verbose=TRUE`.

If [calcObservedMutations](calcObservedMutations.md) was previously run on `db` and saved, supplying
`$pos` from the saved result as `dbMutList` could save time by skipping a
second call of [calcObservedMutations](calcObservedMutations.md). This could be helpful especially when 
`db` is large.



Examples
-------------------

```R
# Load and subset example data
data(ExampleDb, package="alakazam")
db <- ExampleDb[1:5, ]

# Try out thresholds of 2-4 mutations in window sizes of 7-9 nucleotides. 
# In this case, all combinations are legal.
slideWindowTune(db, mutThreshRange=2:4, windowSizeRange=7:9)

```


```
Identifying mutated positions
  |                                                                                                                                               |                                                                                                                                       |   0%  |                                                                                                                                               |===========================                                                                                                            |  20%  |                                                                                                                                               |======================================================                                                                                 |  40%  |                                                                                                                                               |=================================================================================                                                      |  60%  |                                                                                                                                               |============================================================================================================                           |  80%  |                                                                                                                                               |=======================================================================================================================================| 100%  |                                                                                                                                               |                                                                                                                                       |   0%
Analyzing combinations of windowSizeRange and mutThreshRange
  |                                                                                                                                               |===============                                                                                                                        |  11%  |                                                                                                                                               |==============================                                                                                                         |  22%  |                                                                                                                                               |=============================================                                                                                          |  33%  |                                                                                                                                               |============================================================                                                                           |  44%  |                                                                                                                                               |===========================================================================                                                            |  56%  |                                                                                                                                               |==========================================================================================                                             |  67%  |                                                                                                                                               |=========================================================================================================                              |  78%  |                                                                                                                                               |========================================================================================================================               |  89%  |                                                                                                                                               |=======================================================================================================================================| 100%
```


```
$`7`
         2     3     4
[1,]  TRUE FALSE FALSE
[2,] FALSE FALSE FALSE
[3,]  TRUE FALSE FALSE
[4,]  TRUE FALSE FALSE
[5,] FALSE FALSE FALSE

$`8`
         2     3     4
[1,]  TRUE FALSE FALSE
[2,] FALSE FALSE FALSE
[3,]  TRUE FALSE FALSE
[4,]  TRUE FALSE FALSE
[5,] FALSE FALSE FALSE

$`9`
         2     3     4
[1,]  TRUE FALSE FALSE
[2,] FALSE FALSE FALSE
[3,]  TRUE FALSE FALSE
[4,]  TRUE FALSE FALSE
[5,] FALSE FALSE FALSE


```


```R

# Illegal combinations are skipped, returning NAs.
slideWindowTune(db, mutThreshRange=2:4, windowSizeRange=2:4, 
verbose=FALSE)

```


```
  |                                                                                                                                               |                                                                                                                                       |   0%  |                                                                                                                                               |===========================                                                                                                            |  20%  |                                                                                                                                               |======================================================                                                                                 |  40%  |                                                                                                                                               |=================================================================================                                                      |  60%  |                                                                                                                                               |============================================================================================================                           |  80%  |                                                                                                                                               |=======================================================================================================================================| 100%  |                                                                                                                                               |                                                                                                                                       |   0%  |                                                                                                                                               |===============                                                                                                                        |  11%  |                                                                                                                                               |==============================                                                                                                         |  22%  |                                                                                                                                               |=============================================                                                                                          |  33%  |                                                                                                                                               |============================================================                                                                           |  44%  |                                                                                                                                               |===========================================================================                                                            |  56%  |                                                                                                                                               |==========================================================================================                                             |  67%  |                                                                                                                                               |=========================================================================================================                              |  78%  |                                                                                                                                               |========================================================================================================================               |  89%  |                                                                                                                                               |=======================================================================================================================================| 100%
```


```
$`2`
         2  3  4
[1,] FALSE NA NA
[2,] FALSE NA NA
[3,] FALSE NA NA
[4,] FALSE NA NA
[5,] FALSE NA NA

$`3`
         2     3  4
[1,] FALSE FALSE NA
[2,] FALSE FALSE NA
[3,]  TRUE FALSE NA
[4,]  TRUE FALSE NA
[5,] FALSE FALSE NA

$`4`
         2     3     4
[1,]  TRUE FALSE FALSE
[2,] FALSE FALSE FALSE
[3,]  TRUE FALSE FALSE
[4,]  TRUE FALSE FALSE
[5,] FALSE FALSE FALSE


```


```R

# Run calcObservedMutations separately
exDbMutList <- sapply(1:5, function(i) {
calcObservedMutations(inputSeq=db[["sequence_alignment"]][i],
germlineSeq=db[["germline_alignment_d_mask"]][i],
returnRaw=TRUE)$pos })
slideWindowTune(db, dbMutList=exDbMutList, 
mutThreshRange=2:4, windowSizeRange=2:4)
```


```
dbMutList supplied; skipped calling calcObservedMutations()
  |                                                                                                                                               |                                                                                                                                       |   0%
Analyzing combinations of windowSizeRange and mutThreshRange
  |                                                                                                                                               |===============                                                                                                                        |  11%  |                                                                                                                                               |==============================                                                                                                         |  22%  |                                                                                                                                               |=============================================                                                                                          |  33%  |                                                                                                                                               |============================================================                                                                           |  44%>>> mutThresh = 3 > windowSize = 2 (skipped)
  |                                                                                                                                               |===========================================================================                                                            |  56%  |                                                                                                                                               |==========================================================================================                                             |  67%  |                                                                                                                                               |=========================================================================================================                              |  78%>>> mutThresh = 4 > windowSize = 2 (skipped)
  |                                                                                                                                               |========================================================================================================================               |  89%>>> mutThresh = 4 > windowSize = 3 (skipped)
  |                                                                                                                                               |=======================================================================================================================================| 100%
```


```
$`2`
         2  3  4
[1,] FALSE NA NA
[2,] FALSE NA NA
[3,] FALSE NA NA
[4,] FALSE NA NA
[5,] FALSE NA NA

$`3`
         2     3  4
[1,] FALSE FALSE NA
[2,] FALSE FALSE NA
[3,]  TRUE FALSE NA
[4,]  TRUE FALSE NA
[5,] FALSE FALSE NA

$`4`
         2     3     4
[1,]  TRUE FALSE FALSE
[2,] FALSE FALSE FALSE
[3,]  TRUE FALSE FALSE
[4,]  TRUE FALSE FALSE
[5,] FALSE FALSE FALSE


```



See also
-------------------

[slideWindowDb](slideWindowDb.md) is called on `db` for tuning. See [slideWindowTunePlot](slideWindowTunePlot.md) 
for visualization. See [calcObservedMutations](calcObservedMutations.md) for generating `dbMutList`.






