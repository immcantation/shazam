





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
(Corretions at http://selection.med.yale.edu/baseline/correction/)
 



Examples
-------------------

```R
# Subset example data
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, ISOTYPE == "IgG")

# Collapse clones
db <- collapseClones(db, sequenceColumn="SEQUENCE_IMGT",
germlineColumn="GERMLINE_IMGT_D_MASK",
method="thresholdedFreq", minimumFrequency=0.6,
includeAmbiguous=FALSE, breakTiesStochastic=FALSE)

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

# Group PDFs by the sample identifier
grouped <- groupBaseline(baseline, groupBy="SAMPLE")

```


```
Grouping BASELINe probability density functions...
Calculating BASELINe statistics...

```


```R

# Perform test on sample PDFs
testBaseline(grouped, groupBy="SAMPLE")
```


```
  REGION       TEST    PVALUE       FDR
1    CDR -1h != +7d 0.1786408 0.1786408
2    FWR -1h != +7d 0.1155223 0.1786408

```



See also
-------------------

To generate the [Baseline](Baseline-class.md) input object see [groupBaseline](groupBaseline.md).



