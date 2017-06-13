





**plotBaselineSummary** - *Plots BASELINe summary statistics*

Description
--------------------

`plotBaselineSummary` plots a summary of the results of selection analysis 
using the BASELINe method.


Usage
--------------------
```
plotBaselineSummary(baseline, idColumn, groupColumn = NULL,
groupColors = NULL, subsetRegions = NULL, facetBy = c("region",
"group"), title = NULL, style = c("mean"), size = 1, silent = FALSE,
...)
```

Arguments
-------------------

baseline
:   either a data.frame returned from [summarizeBaseline](summarizeBaseline.md)
or a `Baseline` object returned from [groupBaseline](groupBaseline.md)
containing selection probability density functions and summary 
statistics.

idColumn
:   name of the column in `baseline` containing primary identifiers. 
If the input is a `Baseline` object, then this will be a column
in the `stats` slot of `baseline`.

groupColumn
:   name of the column in `baseline` containing secondary grouping 
identifiers. If the input is a `Baseline` object, then this will 
be a column in the `stats` slot of `baseline`.

groupColors
:   named vector of colors for entries in `groupColumn`, with 
names defining unique values in the `groupColumn` and values
being colors. Also controls the order in which groups appear on the
plot. If `NULL` alphabetical ordering and a default color palette 
will be used. Has no effect if `facetBy="group"`.

subsetRegions
:   character vector defining a subset of regions to plot, correspoding 
to the regions for which the `baseline` data was calculated. If
`NULL` all regions in `baseline` are plotted.

facetBy
:   one of c("group", "region") specifying which category to facet the
plot by, either values in `groupColumn` ("group") or regions
defined in `baseline` ("region"). The data that is not used
for faceting will be color coded.

title
:   string defining the plot title.

style
:   type of plot to draw. One of:

+  `"mean"`:     plots the mean and confidence interval for
the selection scores of each value in 
`idColumn`. Faceting and coloring
are determine by values in `groupColumn`
and regions defined in `baseline`, 
depending upon the `facetBy` argument.


size
:   numeric scaling factor for lines, points and text in the plot.

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
# Subset example data
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")

# Calculate BASELINe
baseline <- calcBaseline(db, 
sequenceColumn="SEQUENCE_IMGT",
germlineColumn="GERMLINE_IMGT_D_MASK", 
testStatistic="focused",
regionDefinition=IMGT_V,
targetingModel=HH_S5F,
nproc=1)

```


```
Calculating the expected frequencies of mutations...
Calculating BASELINe probability density functions...

```


```R
 
# Grouping the PDFs by sample and isotype annotations
grouped <- groupBaseline(baseline, groupBy=c("SAMPLE", "ISOTYPE"))

```


```
Grouping BASELINe probability density functions...
Calculating BASELINe statistics...

```


```R

# Plot mean and confidence interval by region with custom group colors
isotype_colors <- c("IgM"="darkorchid", "IgD"="firebrick", 
"IgG"="seagreen", "IgA"="steelblue")
plotBaselineSummary(grouped, "SAMPLE", "ISOTYPE", 
groupColors=isotype_colors)

```

![6](plotBaselineSummary-6.png)

```R

# Facet by group instead of region
plotBaselineSummary(grouped, "SAMPLE", "ISOTYPE", facetBy="group")
```

![8](plotBaselineSummary-8.png)


See also
-------------------

Takes as input either a [Baseline](Baseline-class.md) object returned by [groupBaseline](groupBaseline.md) 
or a data.frame returned from [summarizeBaseline](summarizeBaseline.md).



