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

A [Baseline](Baseline-class.md) object, containing the modified `db` and the BASELINe 
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
(Corrections at http://selection.med.yale.edu/baseline/correction/)
 



Examples
-------------------

```R
### Not run:
# Subset example data from alakazam as a demo
# data(ExampleDb, package="alakazam")
# db <- subset(ExampleDb, c_call %in% c("IGHM", "IGHG"))
# set.seed(112)
# db <- dplyr::slice_sample(db, n=200)
# 
# # Collapse clones
# db <- collapseClones(db, cloneColumn="clone_id",
# sequenceColumn="sequence_alignment",
# germlineColumn="germline_alignment_d_mask",
# method="thresholdedFreq", minimumFrequency=0.6,
# includeAmbiguous=FALSE, breakTiesStochastic=FALSE)
# 
# # Calculate BASELINe
# baseline <- calcBaseline(db, 
# sequenceColumn="clonal_sequence",
# germlineColumn="clonal_germline", 
# testStatistic="focused",
# regionDefinition=IMGT_V,
# targetingModel=HH_S5F,
# nproc=1)
# 
# # Group PDFs by sample
# grouped1 <- groupBaseline(baseline, groupBy="sample_id")
# sample_colors <- c("-1h"="steelblue", "+7d"="firebrick")
# plotBaselineDensity(grouped1, idColumn="sample_id", colorValues=sample_colors, 
# sigmaLimits=c(-1, 1))
#  
# # Group PDFs by both sample (between variable) and isotype (within variable)
# grouped2 <- groupBaseline(baseline, groupBy=c("sample_id", "c_call"))
# isotype_colors <- c("IGHM"="darkorchid", "IGHD"="firebrick", 
# "IGHG"="seagreen", "IGHA"="steelblue")
# plotBaselineDensity(grouped2, idColumn="sample_id", groupColumn="c_call",
# colorElement="group", colorValues=isotype_colors,
# sigmaLimits=c(-1, 1))
# # Collapse previous isotype (within variable) grouped PDFs into sample PDFs
# grouped3 <- groupBaseline(grouped2, groupBy="sample_id")
# sample_colors <- c("-1h"="steelblue", "+7d"="firebrick")
# plotBaselineDensity(grouped3, idColumn="sample_id", colorValues=sample_colors,
# sigmaLimits=c(-1, 1))

```



See also
-------------------

To generate the [Baseline](Baseline-class.md) object see [calcBaseline](calcBaseline.md).
To calculate BASELINe statistics, such as the mean selection strength
and the 95% confidence interval, see [summarizeBaseline](summarizeBaseline.md).






