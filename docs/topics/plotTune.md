**plotTune** - *Visualize parameter tuning for minNumMutations and minNumSeqMutations*

Description
--------------------

Visualize results from [minNumMutationsTune](minNumMutationsTune.md) and [minNumSeqMutationsTune](minNumSeqMutationsTune.md)


Usage
--------------------
```
plotTune(
tuneMtx,
thresh,
criterion = c("5mer", "3mer", "1mer", "3mer+1mer", "measured", "inferred"),
pchs = 1,
ltys = 2,
cols = 1,
plotLegend = TRUE,
legendPos = "topright",
legendHoriz = FALSE,
legendCex = 1
)
```

Arguments
-------------------

tuneMtx
:   a `matrix` or a `list` of matrices produced by either 
[minNumMutationsTune](minNumMutationsTune.md) or [minNumSeqMutationsTune](minNumSeqMutationsTune.md).
In the case of a list, it is assumed that each matrix corresponds
to a sample and that all matrices in the list were produced using
the same set of trial values of `minNumMutations` or 
`minNumSeqMutations`.

thresh
:   a number or a vector of indicating the value or the range of values
of `minNumMutations` or `minNumSeqMutations` to plot. 
Should correspond to the columns of `tuneMtx`.

criterion
:   one of `"5mer"`, `"3mer"`, `"1mer"`, or `"3mer+1mer"` 
(for `tuneMtx` produced by [minNumMutationsTune](minNumMutationsTune.md)), or either 
`"measured"` or `"inferred"` (for `tuneMtx` produced by 
[minNumSeqMutationsTune](minNumSeqMutationsTune.md)).

pchs
:   point types to pass on to [plot](http://www.rdocumentation.org/packages/graphics/topics/plot.default).

ltys
:   line types to pass on to [plot](http://www.rdocumentation.org/packages/graphics/topics/plot.default).

cols
:   colors to pass on to [plot](http://www.rdocumentation.org/packages/graphics/topics/plot.default).

plotLegend
:   whether to plot legend. Default is `TRUE`. Only applicable 
if `tuneMtx` is a named list with names of the matrices 
corresponding to the names of the samples.

legendPos
:   position of legend to pass on to [legend](http://www.rdocumentation.org/packages/graphics/topics/legend). Can be either a
numeric vector specifying x-y coordinates, or one of 
`"topright"`, `"center"`, etc. Default is `"topright"`.

legendHoriz
:   whether to make legend horizontal. Default is `FALSE`.

legendCex
:   numeric values by which legend should be magnified relative to 1.




Details
-------------------

For `tuneMtx` produced by [minNumMutationsTune](minNumMutationsTune.md), for each sample, depending on
`criterion`, the numbers of 5-mers for which substitution rates are directly computed
(`"5mer"`), inferred based on inner 3-mers (`"3mer"`), inferred based on 
central 1-mers (`"1mer"`), or inferred based on inner 3-mers and central 1-mers
(`"3mer+1mer"`) are plotted on the y-axis against values of `minNumMutations` 
on the x-axis.

For `tuneMtx` produced by [minNumSeqMutationsTune](minNumSeqMutationsTune.md), for each sample, depending on
`criterion`, the numbers of 5-mers for which mutability rates are directly measured
(`"measured"`) or inferred (`"inferred"`) are plotted on the y-axis against values
of `minNumSeqMutations` on the x-axis.

Note that legends will be plotted only if `tuneMtx` is a supplied as a named `list`
of matrices, ideally with names of each `matrix` corresponding to those of the samples 
based on which the matrices were produced, even if `plotLegend=TRUE`.



Examples
-------------------

```R
# Subset example data to one isotype and 200 sequences
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, c_call == "IGHA")
set.seed(112)
db <- dplyr::slice_sample(db, n=50)

tuneMtx = list()
for (i in 1:length(unique(db$sample_id))) {
# Get data corresponding to current sample
curDb = db[db[["sample_id"]] == unique(db[["sample_id"]])[i], ]

# Count the number of mutations per 5-mer
subCount = createSubstitutionMatrix(db=curDb, model="s", 
sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask",
vCallColumn="v_call",
multipleMutation="independent",
returnModel="5mer", numMutationsOnly=TRUE)

# Tune over minNumMutations = 5..50
subTune = minNumMutationsTune(subCount, seq(from=5, to=50, by=5))

tuneMtx = c(tuneMtx, list(subTune))
}

# Name tuneMtx after sample names 
names(tuneMtx) = unique(db[["sample_id"]])

# plot with legend for both samples for a subset of minNumMutations values
plotTune(tuneMtx, thresh=c(5, 15, 25, 40), criterion="3mer",
pchs=16:17, ltys=1:2, cols=2:3, 
plotLegend=TRUE, legendPos=c(25, 30))

```

![2](plotTune-2.png)

```R

# plot for only 1 sample for all the minNumMutations values (no legend)
plotTune(tuneMtx[[1]], thresh=seq(from=5, to=50, by=5), criterion="3mer")

```

![4](plotTune-4.png)


See also
-------------------

See [minNumMutationsTune](minNumMutationsTune.md) and [minNumSeqMutationsTune](minNumSeqMutationsTune.md) for generating 
`tuneMtx`.






