





**testBaseline** - *Two-sided test of BASELINe PDFs*

Description
--------------------

`testBaseline` performs a two-sample signifance test of BASELINe 
posterior probability density functions (PDFs).


Usage
--------------------
```
testBaseline(baseline, groupBy)
```

Arguments
-------------------

baseline
:   `Baseline` object containing the `db` and grouped 
BASELINe PDFs returned by [groupBaseline](groupBaseline.md).

groupBy
:   string defining the column in the `db` slot of the 
`Baseline` containing sequence or group identifiers.



Value
-------------------

A data.frame with test results containing the following columns:

+ `REGION`:  sequence region, such as "CDR" and "FWR".
+ `TEST`:    string defining the two group values compared.
+ `PVALUE`:  two-sided p-value for the comparison.
+ `FDR`:     FDR corrected `PVALUE`.


References
-------------------


1. Yaari G, et al. Quantifying selection in high-throughput immunoglobulin 
sequencing data sets. 
Nucleic Acids Res. 2012 40(17):e134.
 



Examples
-------------------

```R
# Subset example data
db <- subset(InfluenzaDb, CPRIMER %in% c("IGHA") & 
BARCODE %in% c("RL016","RL019","RL021"))

# Calculate BASELINe
baseline <- calcBaseline(db, 
sequenceColumn="SEQUENCE_IMGT",
germlineColumn="GERMLINE_IMGT_D_MASK", 
testStatistic="focused",
regionDefinition=IMGT_V_NO_CDR3,
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

# Group PDFs by the sample barcode column
grouped <- groupBaseline(baseline, groupBy="BARCODE")

```


```
Grouping BASELINe probability density functions...
Calculating BASELINe statistics...

```


```R

# Perform test on barcode PDFs
testBaseline(grouped, groupBy="BARCODE")
```


```
  REGION           TEST     PVALUE       FDR
1    CDR RL019 != RL021 0.48366214 0.4836621
2    CDR RL019 != RL016 0.22876315 0.3853040
3    CDR RL021 != RL016 0.25686930 0.3853040
4    FWR RL019 != RL021 0.06239216 0.1871765
5    FWR RL019 != RL016 0.03626564 0.1871765
6    FWR RL021 != RL016 0.37966816 0.4556018

```



See also
-------------------

To generate the [Baseline](Baseline-class.md) input object see [groupBaseline](groupBaseline.md).



