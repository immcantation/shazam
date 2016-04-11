





**plotBaselineDensity** - *Plots BASELINe probability density functions*

Description
--------------------

`plotBaselineDensity` plots the probability density functions resulting from selection 
analysis using the BASELINe method.

Usage
--------------------

```
plotBaselineDensity(baseline, idColumn, groupColumn = NULL,
groupColors = NULL, subsetRegions = NULL, sigmaLimits = c(-5, 5),
title = NULL, facetBy = c("region", "group"), style = c("density"),
size = 1, silent = FALSE, ...)
```

Arguments
-------------------

baseline
:   `Baseline` object containing selection probability 
density functions.

idColumn
:   name of the column in the `db` slot of `baseline` 
containing primary identifiers.

groupColumn
:   name of the column in the `db` slot of `baseline` 
containing secondary grouping identifiers. If `NULL`, 
organize the plot only on values in `idColumn`.

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

sigmaLimits
:   numeric vector containing two values defining the `c(lower, upper)`
bounds of the selection scores to plot.

title
:   string defining the plot title.

facetBy
:   one of c("group", "region") specifying which category to facet the
plot by, either values in `groupColumn` ("group") or regions
defined in `baseline` ("region"). The data that is not used
for faceting will be color coded.

style
:   type of plot to draw. One of:

+  `"density"`:  plots a set of curves for each probability 
density function in `baseline`, 
with the line type determined by `idColumn`.
Faceting and coloring are determine by values in 
`groupColumn` and regions defined in 
`baseline`, depending upon the 
`facetBy` argument.


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
db <- subset(InfluenzaDb, CPRIMER %in% c("IGHA","IGHM") & 
BARCODE %in% c("RL016","RL018","RL019","RL021"))

# Calculate BASELINe
# By default, calcBaseline collapses the sequences in the db by the column "CLONE",
# calculates the numbers of observed mutations and expected frequencies of mutations,
# as defined in the IMGT_V_NO_CDR3 and using the HS5FModel targeting model.
# Then, it calculates  the BASELINe posterior probability density functions (PDFs) for
# sequences in the updated db files; using the focused test statistic
db_baseline <- calcBaseline(db, 
sequenceColumn="SEQUENCE_IMGT",
germlineColumn="GERMLINE_IMGT_D_MASK", 
testStatistic="focused",
regionDefinition=IMGT_V_NO_CDR3,
targetingModel = HS5FModel,
nproc = 1)

```


```
Collapsing clonal sequences...
Calculating observed number of mutations...
Calculating the expected frequencies of mutations...
Calculating BASELINe probability density functions...

```


```R
 
# Grouping the PDFs by the BARCODE and CPRIMER columns in the db, corresponding 
# respectively to sample barcodes and the constant region isotype primers.
baseline <- groupBaseline(db_baseline, groupBy=c("BARCODE", "CPRIMER"))

```


```
Grouping BASELINe probability density functions...
Calculating BASELINe statistics...

```


```R

# Plot mean and confidence interval
plotBaselineDensity(baseline, "BARCODE", "CPRIMER", style="density")

```

![6](plotBaselineDensity-6.png)

```R
plotBaselineDensity(baseline, "BARCODE", "CPRIMER", subsetRegions="CDR", style="density")

```

![8](plotBaselineDensity-8.png)

```R
plotBaselineDensity(baseline, "BARCODE", "CPRIMER", facetBy="group", style="density")

```

![10](plotBaselineDensity-10.png)

```R

# Reorder and recolor groups
group_colors <- c("IGHM"="darkorchid", "IGHD"="firebrick", "IGHG"="seagreen", "IGHA"="steelblue")
plotBaselineDensity(baseline, "BARCODE", "CPRIMER", groupColors=group_colors, style="density")
```

![12](plotBaselineDensity-12.png)


See also
-------------------

Takes as input a [Baseline](Baseline-class.md) object returned from [groupBaseline](groupBaseline.md).



