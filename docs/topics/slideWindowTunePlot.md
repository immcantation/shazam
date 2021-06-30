**slideWindowTunePlot** - *Visualize parameter tuning for sliding window approach*

Description
--------------------

Visualize results from [slideWindowTune](slideWindowTune.md)


Usage
--------------------
```
slideWindowTunePlot(
tuneList,
plotFiltered = TRUE,
percentage = FALSE,
jitter.x = FALSE,
jitter.x.amt = 0.1,
jitter.y = FALSE,
jitter.y.amt = 0.1,
pchs = 1,
ltys = 2,
cols = 1,
plotLegend = TRUE,
legendPos = "topright",
legendHoriz = FALSE,
legendCex = 1,
title = NULL
)
```

Arguments
-------------------

tuneList
:   a list of logical matrices returned by [slideWindowTune](slideWindowTune.md).

plotFiltered
:   whether to plot the number of filtered sequences (as opposed to
the number of remaining sequences). Default is `TRUE`.

percentage
:   whether to plot on the y-axis the percentage of filtered sequences
(as opposed to the absolute number). Default is `FALSE`.

jitter.x
:   whether to jitter x-axis values. Default is `FALSE`.

jitter.x.amt
:   amount of jittering to be applied on x-axis values if 
`jitter.x=TRUE`. Default is 0.1.

jitter.y
:   whether to jitter y-axis values. Default is `FALSE`.

jitter.y.amt
:   amount of jittering to be applied on y-axis values if 
`jitter.y=TRUE`. Default is 0.1.

pchs
:   point types to pass on to [plot](http://www.rdocumentation.org/packages/base/topics/plot).

ltys
:   line types to pass on to [plot](http://www.rdocumentation.org/packages/base/topics/plot).

cols
:   colors to pass on to [plot](http://www.rdocumentation.org/packages/base/topics/plot).

plotLegend
:   whether to plot legend. Default is `TRUE`.

legendPos
:   position of legend to pass on to [legend](http://www.rdocumentation.org/packages/graphics/topics/legend). Can be either a
numeric vector specifying x-y coordinates, or one of 
`"topright"`, `"center"`, etc. Default is `"topright"`.

legendHoriz
:   whether to make legend horizontal. Default is `FALSE`.

legendCex
:   numeric values by which legend should be magnified relative to 1.

title
:   plot main title. Default is NULL (no title)




Details
-------------------

For each `windowSize`, the numbers of sequences filtered or remaining after applying
the sliding window approach are plotted on the y-axis against thresholds on the number of
mutations in a window on the x-axis.

When plotting, a user-defined `amount` of jittering can be applied on values plotted
on either axis or both axes via adjusting `jitter.x`, `jitter.y`, 
`jitter.x.amt` and `jitter.y.amt`. This may be help with visually distinguishing
lines for different window sizes in case they are very close or identical to each other. 
If plotting percentages (`percentage=TRUE`) and using jittering on the y-axis values 
(`jitter.y=TRUE`), it is strongly recommended that `jitter.y.amt` be set very
small (e.g. 0.01). 

`NA` for a combination of `mutThresh` and `windowSize` where 
`mutThresh` is greater than `windowSize` will not be plotted.



Examples
-------------------

```R
# Use an entry in the example data for input and germline sequence
data(ExampleDb, package="alakazam")

# Try out thresholds of 2-4 mutations in window sizes of 3-5 nucleotides 
# on a subset of ExampleDb
tuneList <- slideWindowTune(db = ExampleDb[1:10, ], 
mutThreshRange = 2:4, windowSizeRange = 3:5,
verbose = FALSE)

# Visualize
# Plot numbers of sequences filtered without jittering y-axis values
slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3, 
plotFiltered=TRUE, jitter.y=FALSE)

```

![2](slideWindowTunePlot-2.png)

```R

# Notice that some of the lines overlap
# Jittering could help
slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3,
plotFiltered=TRUE, jitter.y=TRUE)

```

![4](slideWindowTunePlot-4.png)

```R

# Plot numbers of sequences remaining instead of filtered
slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3, 
plotFiltered=FALSE, jitter.y=TRUE, 
legendPos="bottomright")

```

![6](slideWindowTunePlot-6.png)

```R

# Plot percentages of sequences filtered with a tiny amount of jittering
slideWindowTunePlot(tuneList, pchs=1:3, ltys=1:3, cols=1:3,
plotFiltered=TRUE, percentage=TRUE, 
jitter.y=TRUE, jitter.y.amt=0.01)
```

![8](slideWindowTunePlot-8.png)


See also
-------------------

See [slideWindowTune](slideWindowTune.md) for how to get `tuneList`. See [jitter](http://www.rdocumentation.org/packages/base/topics/jitter) for 
use of `amount` of jittering.






