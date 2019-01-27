**createBaseline** - *Creates a Baseline object*

Description
--------------------

`createBaseline` creates and initialize a `Baseline` object.


Usage
--------------------
```
createBaseline(description = "", db = data.frame(),
regionDefinition = createRegionDefinition(), testStatistic = "",
regions = NULL, numbOfSeqs = matrix(), binomK = matrix(),
binomN = matrix(), binomP = matrix(), pdfs = list(),
stats = data.frame())
```

Arguments
-------------------

description
:   `character` providing general information regarding the 
sequences, selection analysis and/or object.

db
:   `data.frame` containing annotation information about 
the sequences and selection results.

regionDefinition
:   [RegionDefinition](RegionDefinition-class.md) object defining the regions
and boundaries of the Ig sequences.

testStatistic
:   `character` indicating the statistical framework 
used to test for selection. For example, `"local"` or 
`"focused"` or `"imbalanced"`.

regions
:   `character` vector defining the regions the BASELINe 
analysis was carried out on. For `"CDR"` and `"FWR"` 
or `"CDR1"`, `"CDR2"`, `"CDR3"`, etc. If `NULL`
then regions will be determined automatically from `regionDefinition`.

numbOfSeqs
:   `matrix` of dimensions `r x c` containing the number of 
sequences or PDFs in each region, where:
`r` = number of rows = number of groups or sequences.
`c` = number of columns = number of regions.

binomK
:   `matrix` of dimensions `r x c` containing the number of 
successes in the binomial trials in each region, where:
`r` = number of rows = number of groups or sequences.
`c` = number of columns = number of regions.

binomN
:   `matrix` of dimensions `r x c` containing the total 
number of trials in the binomial in each region, where:
`r` = number of rows = number of groups or sequences.
`c` = number of columns = number of regions.

binomP
:   `matrix` of dimensions `r x c` containing the probability 
of success in one binomial trial in each region, where:
`r` = number of rows = number of groups or sequences.
`c` = number of columns = number of regions.

pdfs
:   `list` of matrices containing PDFs with one item for each 
defined region (e.g. "CDR" and "FWR"). Matrices have dimensions
`r x c` dementions, where:
`r` = number of rows = number of sequences or groups. 
`c` = number of columns = length of the PDF (default 4001).

stats
:   `data.frame` of BASELINe statistics, 
including: mean selection strength (mean Sigma), 95% confidence 
intervals, and p-values with positive signs for the presence of 
positive selection and/or p-values with negative signs for the
presence of negative selection.




Value
-------------------

A `Baseline` object.


Details
-------------------

Create and initialize a `Baseline` object. 

The `testStatistic` indicates the statistical framework used to test for selection. 
For example,

+ `local` = CDR_R / (CDR_R + CDR_S).
+ `focused` = CDR_R / (CDR_R + CDR_S + FWR_S).
+ `immbalance` = CDR_R + CDR_s / (CDR_R + CDR_S + FWR_S + FWR_R)

For `focused` the `regionDefinition` must only contain two regions. If more 
than two regions are defined, then the `local` test statistic will be used.
For further information on the frame of these tests see Uduman et al. (2011).


References
-------------------


1. Hershberg U, et al. Improved methods for detecting selection by mutation 
analysis of Ig V region sequences. 
Int Immunol. 2008 20(5):683-94.
1. Uduman M, et al. Detecting selection in immunoglobulin sequences. 
Nucleic Acids Res. 2011 39(Web Server issue):W499-504.
1. Yaari G, et al. Models of somatic hypermutation targeting and substitution based
on synonymous mutations from high-throughput immunoglobulin sequencing data.
Front Immunol. 2013 4(November):358.
 



Examples
-------------------

```R
# Creates an empty Baseline object
createBaseline()
```


```
An object of class "Baseline"
Slot "description":
[1] ""

Slot "db":
data frame with 0 columns and 0 rows

Slot "regionDefinition":
An object of class "RegionDefinition"
Slot "name":
[1] ""

Slot "description":
[1] ""

Slot "boundaries":
factor(0)
Levels: 

Slot "seqLength":
[1] 0

Slot "regions":
character(0)

Slot "labels":
character(0)

Slot "citation":
[1] ""


Slot "testStatistic":
[1] ""

Slot "regions":
character(0)

Slot "numbOfSeqs":
     [,1]
[1,]   NA

Slot "binomK":
     [,1]
[1,]   NA

Slot "binomN":
     [,1]
[1,]   NA

Slot "binomP":
     [,1]
[1,]   NA

Slot "pdfs":
list()

Slot "stats":
[1] GROUP              REGION             BASELINE_SIGMA     BASELINE_CI_LOWER 
[5] BASELINE_CI_UPPER  BASELINE_CI_PVALUE
<0 rows> (or 0-length row.names)


```



See also
-------------------

See [Baseline](Baseline-class.md) for the return object.



