





**getBaselineStats** - *Gets the summary statistics of a Baseline object*

Description
--------------------

`getBaselineStats` is an accessor method that returns the 
summary statistics `data.frame` stored in the `stats` slot of a 
[Baseline](Baseline-class.md) object - provided [groupBaseline](groupBaseline.md) has already been run.


Usage
--------------------
```
getBaselineStats(baseline)
```

Arguments
-------------------

baseline
:   `Baseline` object that has been run through
either [groupBaseline](groupBaseline.md) or [summarizeBaseline](summarizeBaseline.md).




Value
-------------------

A `data.frame` with the BASELINe selection strength scores (Sigma),
95% confidence intervals and P-values.



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
targetingModel = HH_S5F,
nproc = 1)

```


```
Collapsing clonal sequences...
Calculating observed number of mutations...
Calculating the expected frequencies of mutations...
Calculating BASELINe probability density functions...

```


```R

# Grouping the PDFs by the isotype and sample annotations.
grouped <- groupBaseline(baseline, groupBy=c("SAMPLE", "ISOTYPE"))

```


```
Grouping BASELINe probability density functions...
Calculating BASELINe statistics...

```


```R

# Get a data.frame of the summary statistics
getBaselineStats(grouped)
```


```
  SAMPLE ISOTYPE REGION BASELINE_SIGMA BASELINE_CI_LOWER BASELINE_CI_UPPER BASELINE_CI_PVALUE
1    +7d     IgA    CDR     -0.3066403        -0.4981192        -0.1460119      -5.499097e-02
2    +7d     IgA    FWR     -0.6990904        -0.8265410        -0.5797970      -2.842171e-14
3    +7d     IgG    CDR     -0.2353832        -0.3654914        -0.1244890      -2.420408e-02
4    +7d     IgG    FWR     -0.7065625        -0.8043858        -0.6189482      -5.595524e-14

```



See also
-------------------

For calculating the BASELINe summary statistics see [summarizeBaseline](summarizeBaseline.md).



