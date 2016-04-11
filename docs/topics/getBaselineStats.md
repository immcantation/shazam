





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
nproc = 1)

```


```
Collapsing clonal sequences...
Calculating observed number of mutations...
Calculating the expected frequencies of mutations...
Calculating BASELINe probability density functions...

```


```R

# Grouping the PDFs by the BARCODE and CPRIMER columns in the db, corresponding 
# respectively to sample barcodes and the constant region isotype primers.
baseline_group <- groupBaseline(db_baseline, groupBy=c("BARCODE", "CPRIMER"))

```


```
Grouping BASELINe probability density functions...
Calculating BASELINe statistics...

```


```R

# Get a data.frame of the summary statistics
getBaselineStats(baseline_group)
```


```
   BARCODE CPRIMER REGION BASELINE_SIGMA BASELINE_CI_LOWER BASELINE_CI_UPPER BASELINE_CI_PVALUE
1    RL018    IGHA    CDR    -0.09930603        -0.5791206        0.38318372       1.499988e+01
2    RL018    IGHA    FWR    -0.31409259        -0.7076840        0.13094124       5.588523e+00
3    RL019    IGHM    CDR     0.29386108        -0.6382929        1.22729363       6.282453e+00
4    RL019    IGHM    FWR     0.25646516        -1.4405725        2.92533433       3.765002e+00
5    RL018    IGHM    CDR    -1.11620070        -3.5739877        0.67146613       2.450266e+00
6    RL018    IGHM    FWR    -0.70658387        -2.4171269        0.84514743       3.259269e+00
7    RL016    IGHM    CDR     2.58001794        -0.5851458        7.26794234      -5.694666e-03
8    RL016    IGHM    FWR     1.53241378        -2.6650480        8.53741764       1.122445e+00
9    RL019    IGHA    CDR    -0.30802196        -0.6562439        0.01435708       4.470269e+00
10   RL019    IGHA    FWR    -0.96736833        -1.2627322       -0.68476224       2.835981e-08
11   RL021    IGHM    CDR    -0.59639368        -2.6998350        1.29307900       3.448511e+00
12   RL021    IGHM    FWR     0.11198597        -1.2685964        1.64717439       5.237258e+00
13   RL021    IGHA    CDR    -0.34235651        -0.8791342        0.12815185       6.335781e+00
14   RL021    IGHA    FWR    -0.53679555        -1.0692065        0.16551298       2.340685e+00
15   RL016    IGHA    CDR    -0.11631575        -0.5934801        0.35003487       1.530050e+01
16   RL016    IGHA    FWR    -0.64435013        -0.9273246       -0.34483438       2.600630e-02

```



See also
-------------------

For calculating the BASELINe summary statistics see [summarizeBaseline](summarizeBaseline.md).



