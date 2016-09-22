





**groupBaseline** - *Group BASELINe PDFs*

Description
--------------------

`groupBaseline` convolves groups of BASELINe posterior probability density 
functions (PDFs) to get combined PDFs for each group.


Usage
--------------------
```
groupBaseline(baseline, groupBy, nproc = 1)
```

Arguments
-------------------

baseline
:   `Baseline` object containing the `db` and the 
BASELINe posterior probability density functions 
(PDF) for each of the sequences, as returned by
[calcBaseline](calcBaseline.md).

groupBy
:   The columns in the `db` slot of the `Baseline`
object by which to group the sequence PDFs.

nproc
:   number of cores to distribute the operation over. If 
`nproc` = 0 then the `cluster` has already been
set and will not be reset.



Value
-------------------

A `Baseline` object, containing the modified `db` and the BASELINe 
posterior probability density functions (PDF) for each of the groups.

Details
-------------------

While the selection strengths predicted by BASELINe perform well on average, 
the estimates for individual sequences can be highly variable, especially when the 
number of mutations is small. 

To overcome this, PDFs from sequences grouped by biological or experimental relevance,
are convolved to from a single PDF for the selection strength. For example, sequences
from each sample may be combined together, allowing you to compare selection  across 
samples. This is accomplished through a fast numerical convolution technique.

References
-------------------


1. Yaari G, et al. Quantifying selection in high-throughput immunoglobulin 
sequencing data sets. 
Nucleic Acids Res. 2012 40(17):e134.
 



Examples
-------------------

```R
# Subset example data from alakazam
library(alakazam)
db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG"))

# Calculate BASELINe
baseline <- calcBaseline(db, 
sequenceColumn="SEQUENCE_IMGT",
germlineColumn="GERMLINE_IMGT_D_MASK", 
testStatistic="focused",
regionDefinition=IMGT_V,
targetingModel=HS5FModel,
nproc=1)

```


```
Collapsing clonal sequences...
Calculating observed number of mutations...
Calculating the expected frequencies of mutations...
Calculating BASELINe probability density functions...

```


```R

# Group PDFs by sample
grouped1 <- groupBaseline(baseline, groupBy="SAMPLE")

```


```
Grouping BASELINe probability density functions...
Calculating BASELINe statistics...

```


```R
plotBaselineDensity(grouped1, idColumn="SAMPLE", colorElement="group", 
sigmaLimits=c(-1, 1))

```

![6](groupBaseline-6.png)

```R
 
# Group PDFs by both sample (between variable) and isotype (within variable)
grouped2 <- groupBaseline(baseline, groupBy=c("SAMPLE", "ISOTYPE"))

```


```
Grouping BASELINe probability density functions...
Calculating BASELINe statistics...

```


```R
plotBaselineDensity(grouped2, idColumn="SAMPLE", groupColumn="ISOTYPE",
colorElement="group", colorValues=IG_COLORS,
sigmaLimits=c(-1, 1))

```

![10](groupBaseline-10.png)

```R

# Collapse previous isotype (within variable) grouped PDFs into sample PDFs
grouped3 <- groupBaseline(grouped2, groupBy="SAMPLE")

```


```
Grouping BASELINe probability density functions...
Calculating BASELINe statistics...

```


```R
plotBaselineDensity(grouped3, idColumn="SAMPLE", colorElement="group",
sigmaLimits=c(-1, 1))
```

![14](groupBaseline-14.png)


See also
-------------------

To generate the baseline object see [calcBaseline](calcBaseline.md).
To calculate BASELINe statistics, such as the mean selection strength
and the 95% confidence interval, see [summarizeBaseline](summarizeBaseline.md).



