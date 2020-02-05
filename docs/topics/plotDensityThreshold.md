**plotDensityThreshold** - *Plot findThreshold results for the density method*

Description
--------------------

`plotDensityThreshold` plots the results from `"density"` method of 
[findThreshold](findThreshold.md), including the smoothed density estimate, input nearest neighbor 
distance histogram, and threshold selected.


Usage
--------------------
```
plotDensityThreshold(
data,
cross = NULL,
xmin = NULL,
xmax = NULL,
breaks = NULL,
binwidth = NULL,
title = NULL,
size = 1,
silent = FALSE,
...
)
```

Arguments
-------------------

data
:   [DensityThreshold](DensityThreshold-class.md) object output by the `"density"` method 
of [findThreshold](findThreshold.md).

cross
:   numeric vector of distances from [distToNearest](distToNearest.md) to draw as a
histogram below the `data` histogram for comparison purposes.

xmin
:   minimum limit for plotting the x-axis. If `NULL` the limit will 
be set automatically.

xmax
:   maximum limit for plotting the x-axis. If `NULL` the limit will 
be set automatically.

breaks
:   number of breaks to show on the x-axis. If `NULL` the breaks will 
be set automatically.

binwidth
:   binwidth for the histogram. If `NULL` the binwidth 
will be set automatically to the bandwidth parameter determined by
[findThreshold](findThreshold.md).

title
:   string defining the plot title.

size
:   numeric value for the plot line sizes.

silent
:   if `TRUE` do not draw the plot and just return the ggplot2 
object; if `FALSE` draw the plot.

...
:   additional arguments to pass to ggplot2::theme.




Value
-------------------

A ggplot object defining the plot.



Examples
-------------------

```R
# Subset example data to one sample as a demo
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, SAMPLE == "-1h")

# Use nucleotide Hamming distance and normalize by junction length
db <- distToNearest(db, model="ham", normalize="len", nproc=1)

# To find the threshold cut, call findThreshold function for "gmm" method.
output <- findThreshold(db$DIST_NEAREST, method="density")
print(output)

```


```
[1] NA

```


```R

# Plot
plotDensityThreshold(output)
```

*Warning*:Removed 1 rows containing missing values (geom_vline).![5](plotDensityThreshold-5.png)


See also
-------------------

See [DensityThreshold](DensityThreshold-class.md) for the the input object definition and 
[findThreshold](findThreshold.md) for generating the input object. See 
[distToNearest](distToNearest.md) calculating nearest neighbor distances.






