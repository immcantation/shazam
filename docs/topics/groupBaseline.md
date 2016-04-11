





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
# Subset example data
db <- subset(InfluenzaDb, CPRIMER %in% c("IGHA","IGHM") & 
BARCODE %in% c("RL016","RL018","RL019","RL021"))

# Calculate BASELINe
baseline <- calcBaseline(db, 
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

# Grouping the PDFs by the sample barcode column
baseline_grp1 <- groupBaseline(baseline, groupBy="BARCODE")

```


```
Grouping BASELINe probability density functions...
Calculating BASELINe statistics...

```


```R
 
# Grouping the PDFs by the sample barcode and C-region primer columns
baseline_grp2 <- groupBaseline(baseline, groupBy=c("BARCODE", "CPRIMER"))

```


```
Grouping BASELINe probability density functions...
Calculating BASELINe statistics...

```


```R

# Re-grouping the PDFs by the sample barcode from sample barcode and C-region groups
baseline_grp3 <- groupBaseline(baseline_grp2, groupBy=c("BARCODE"))
```


```
Grouping BASELINe probability density functions...
Calculating BASELINe statistics...

```



See also
-------------------

To generate the baseline object see [calcBaseline](calcBaseline.md).
To calculate BASELINe statistics, such as the mean selection strength
and the 95% confidence interval, see [summarizeBaseline](summarizeBaseline.md).



