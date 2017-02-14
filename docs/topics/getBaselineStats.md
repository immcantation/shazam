





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

A `data.frame` with the mean selection strength (mean Sigma), 95% 
confidence intervals, and p-values with positive signs for the presence 
of positive selection and/or p-values with negative signs for the presence 
of negative selection.



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
targetingModel=HH_S5F,
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
1    +7d     IgA    CDR     -0.3068555        -0.4983114        -0.1461950      -1.153623e-04
2    +7d     IgA    FWR     -0.6975822        -0.8249806        -0.5785985      -1.232348e-14
3    +7d     IgG    CDR     -0.2430318        -0.3730363        -0.1329690      -1.956581e-05
4    +7d     IgG    FWR     -0.7168326        -0.8145677        -0.6290923      -7.105427e-15

```



See also
-------------------

For calculating the BASELINe summary statistics see [summarizeBaseline](summarizeBaseline.md).



