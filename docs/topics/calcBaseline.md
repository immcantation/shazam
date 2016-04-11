





**calcBaseline** - *Calculate the BASELINe PDFs*

Description
--------------------

`calcBaseline` calculates the BASELINe posterior probability density 
functions (PDFs) for sequences in the given Change-O `data.frame`.

Usage
--------------------

```
calcBaseline(db, sequenceColumn = "SEQUENCE_IMGT",
germlineColumn = "GERMLINE_IMGT_D_MASK", testStatistic = c("local",
"focused", "imbalance"), regionDefinition = NULL,
targetingModel = HS5FModel, calcStats = FALSE, nproc = 1)
```

Arguments
-------------------

db
:   `data.frame` containing sequence data and annotations.

sequenceColumn
:   `character` name of the column in `db` 
containing input sequences.

germlineColumn
:   `character` name of the column in `db` 
containing germline sequences.

testStatistic
:   `character` indicating the statistical framework 
used to test for selection. One of `c("local", "focused", "imbalance")`.

regionDefinition
:   [RegionDefinition](RegionDefinition-class.md) object defining the regions
and boundaries of the Ig sequences.

targetingModel
:   [TargetingModel](TargetingModel-class.md) object. Default is  [HS5FModel](HS5FModel.md).

calcStats
:   `logical` indicating whether or not to calculate the 
summary statistics `data.frame` stored in the 
`stats` slot of a [Baseline](Baseline-class.md) object.

nproc
:   number of cores to distribute the operation over. If 
`nproc=0` then the `cluster` has already been
set and will not be reset.



Value
-------------------

A `Baseline` object containing the modified `db` and BASELINe 
posterior probability density functions (PDF) for each of the sequences.

Details
-------------------

Calculates the BASELINe posterior probability density function (PDF) for 
sequences in the provided `db`. 

If the `db` does not contain the 
required columns to calculate the PDFs (namely OBSERVED & EXPECTED mutations)
then the function will:

1. Collapse the sequences by the CLONE column (if present).
1. Calculate the numbers of observed mutations.
1. Calculate the expected frequencies of mutations and modify the provided 
`db`. The modified `db` will be included as part of the 
returned `Baseline` object.


The `testStatistic` indicates the statistical framework used to test for selection. 
E.g.

+ `local` = CDR_R / (CDR_R + CDR_S).
+ `focused` = CDR_R / (CDR_R + CDR_S + FWR_S).
+ `imbalance` = CDR_R + CDR_S / (CDR_R + CDR_S + FWR_S + FRW_R).

For `focused` the `regionDefinition` must only contain two regions. If more 
than two regions are defined the `local` test statistic will be used.
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
nproc=1)
```


```
Collapsing clonal sequences...
Calculating observed number of mutations...
Calculating the expected frequencies of mutations...
Calculating BASELINe probability density functions...

```




