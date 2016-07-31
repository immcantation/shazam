





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
1    +7d     IgA    CDR     -0.3217321        -0.4973101        -0.1794578      -3.315013e-03
2    +7d     IgA    FWR     -0.7226396        -0.8341240        -0.6228244      -3.663736e-14
3    +7d     IgG    CDR     -0.2979155        -0.4029397        -0.2053538      -1.993574e-06
4    +7d     IgG    FWR     -0.7437273        -0.8309610        -0.6671739      -5.351275e-14

```



See also
-------------------

For calculating the BASELINe summary statistics see [summarizeBaseline](summarizeBaseline.md).



