**summarizeBaseline** - *Calculate BASELINe summary statistics*

Description
--------------------

`summarizeBaseline` calculates BASELINe statistics such as the mean selection 
strength (mean Sigma), the 95% confidence intervals and p-values for the presence of
selection.


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
mean BASELINe selection strength, its 95% confidence intervals, and 
a p-value for the presence of selection.


Details
-------------------

The returned p-value can be either positive or negative. Its magnitude 
(without the sign) should be interpreted as per normal. Its sign indicates 
the direction of the selection detected. A positive p-value indicates positive
selection, whereas a negative p-value indicates negative selection.


References
-------------------


1. Uduman M, et al. Detecting selection in immunoglobulin sequences. 
Nucleic Acids Res. 2011 39(Web Server issue):W499-504.




Examples
-------------------

```R
# Subset example data
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, c_call == "IGHG")
set.seed(112)
db <- dplyr::slice_sample(db, n=100)

# Collapse clones
db <- collapseClones(db, cloneColumn="clone_id",
sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask",
method="thresholdedFreq", minimumFrequency=0.6,
includeAmbiguous=FALSE, breakTiesStochastic=FALSE)

# Calculate BASELINe
baseline <- calcBaseline(db, 
sequenceColumn="clonal_sequence",
germlineColumn="clonal_germline", 
testStatistic="focused",
regionDefinition=IMGT_V,
targetingModel=HH_S5F,
nproc = 1)

```

*calcBaseline will calculate observed and expected mutations for clonal_sequence using clonal_germline as a reference.*
```
Calculating BASELINe probability density functions...

```


```R

# Grouping the PDFs by the sample annotation
grouped <- groupBaseline(baseline, groupBy="sample_id")

```


```
Grouping BASELINe probability density functions...
Calculating BASELINe statistics...

```


```R

# Get a data.frame of the summary statistics
stats <- summarizeBaseline(grouped, returnType="df")

```


```
Calculating BASELINe statistics...

```



See also
-------------------

See [calcBaseline](calcBaseline.md) for generating `Baseline` objects and
[groupBaseline](groupBaseline.md) for convolving groups of BASELINe PDFs.






