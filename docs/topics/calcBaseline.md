**calcBaseline** - *Calculate the BASELINe PDFs (including for regions that include CDR3 and FWR4)*

Description
--------------------

`calcBaseline` calculates the BASELINe posterior probability density 
functions (PDFs) for sequences in the given Change-O `data.frame`.


Usage
--------------------
```
calcBaseline(
db,
sequenceColumn = "clonal_sequence",
germlineColumn = "clonal_germline",
testStatistic = c("local", "focused", "imbalanced"),
regionDefinition = NULL,
targetingModel = HH_S5F,
mutationDefinition = NULL,
calcStats = FALSE,
nproc = 1,
cloneColumn = NULL,
juncLengthColumn = NULL
)
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
used to test for selection. One of 
`c("local", "focused", "imbalanced")`.

regionDefinition
:   [RegionDefinition](RegionDefinition-class.md) object defining the regions
and boundaries of the Ig sequences.

targetingModel
:   [TargetingModel](TargetingModel-class.md) object. Default is  [HH_S5F](HH_S5F.md).

mutationDefinition
:   [MutationDefinition](MutationDefinition-class.md) object defining replacement
and silent mutation criteria. If `NULL` then 
replacement and silent are determined by exact 
amino acid identity. Note, if the input data.frame 
already contains observed and expected mutation frequency 
columns then mutations will not be recalculated and this
argument will be ignored.

calcStats
:   `logical` indicating whether or not to calculate the 
summary statistics `data.frame` stored in the 
`stats` slot of a [Baseline](Baseline-class.md) object.

nproc
:   number of cores to distribute the operation over. If 
`nproc=0` then the `cluster` has already been
set and will not be reset.

cloneColumn
:   `character` name of the column in `db` 
containing clonal identifiers. Relevant only for 
when regionDefinition includes CDR and FWR4 (else
this value can be `NULL`)

juncLengthColumn
:   `character` name of the column in `db` 
containing the junction length. Relevant only for 
when regionDefinition includes CDR and FWR4 (else
this value can be `NULL`)




Value
-------------------

A [Baseline](Baseline-class.md) object containing the modified `db` and BASELINe 
posterior probability density functions (PDF) for each of the sequences.


Details
-------------------

Calculates the BASELINe posterior probability density function (PDF) for 
sequences in the provided `db`. 

**Note**: Individual sequences within clonal groups are not, strictly speaking, 
independent events and it is generally appropriate to only analyze selection 
pressures on an effective sequence for each clonal group. For this reason,
it is strongly recommended that the input `db` contains one effective 
sequence per clone. Effective clonal sequences can be obtained by calling 
the [collapseClones](collapseClones.md) function.

If the `db` does not contain the 
required columns to calculate the PDFs (namely mu_count & mu_expected)
then the function will:

1. Calculate the numbers of observed mutations.
1. Calculate the expected frequencies of mutations and modify the provided 
`db`. The modified `db` will be included as part of the 
returned `Baseline` object.


The `testStatistic` indicates the statistical framework used to test for selection. 
E.g.

+ `local` = CDR_R / (CDR_R + CDR_S).
+ `focused` = CDR_R / (CDR_R + CDR_S + FWR_S).
+ `imbalanced` = CDR_R + CDR_S / (CDR_R + CDR_S + FWR_S + FRW_R).

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
# Load and subset example data
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, c_call == "IGHG" & sample_id == "+7d")

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
nproc=1)

```

*calcBaseline will calculate observed and expected mutations for clonal_sequence using clonal_germline as a reference.*
```
Calculating BASELINe probability density functions...

```



See also
-------------------

See [Baseline](Baseline-class.md) for the return object.
See [groupBaseline](groupBaseline.md) and [summarizeBaseline](summarizeBaseline.md) for further processing.
See [plotBaselineSummary](plotBaselineSummary.md) and [plotBaselineDensity](plotBaselineDensity.md) for plotting results.






