





**summarizeBaseline** - *Calculate BASELINe summary statistics*

Description
--------------------

`summarizeBaseline` calculates BASELINe statistics such as the selection strength
(Sigma), the 95% confidence intervals and P-values.

Usage
--------------------

```
summarizeBaseline(baseline, returnType = c("baseline", "df"), nproc = 1)
```

Arguments
-------------------

baseline
:   `Baseline` object returned by [calcBaseline](calcBaseline.md) containing 
annotations and BASELINe posterior probability density functions 
(PDFs) for each sequence.

returnType
:   One of `c("baseline", "df")` defining whether
to return a `Baseline` object ("baseline") with an updated
`stats` slot or a data.frame ("df") of summary statistics.

nproc
:   number of cores to distribute the operation over. If 
`nproc` = 0 then the `cluster` has already been
set and will not be reset.



Value
-------------------

Either a modified `Baseline` object or data.frame containing the 
BASELINe selection strength, 95% confidence intervals and P-value.



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

# Grouping the PDFs by the sample barcode and C-region primer columns
baseline_grp <- groupBaseline(baseline, groupBy=c("BARCODE", "CPRIMER"))

```


```
Grouping BASELINe probability density functions...
Calculating BASELINe statistics...

```


```R

# Get a data.frame of the summary statistics
baseline_stats <- summarizeBaseline(baseline_grp, returnType="df")
```


```
Calculating BASELINe statistics...

```



See also
-------------------

See [calcBaseline](calcBaseline.md) for generating `Baseline` objects and
[groupBaseline](groupBaseline.md) for convolving groups of BASELINe PDFs.



