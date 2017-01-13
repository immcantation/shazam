





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

# Grouping the PDFs by the sample and isotype annotations
grouped <- groupBaseline(baseline, groupBy=c("SAMPLE", "ISOTYPE"))

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



