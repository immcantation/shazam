**createSubstitutionMatrix** - *Builds a substitution model*

Description
--------------------

`createSubstitutionMatrix` builds a 5-mer nucleotide substitution model by counting 
the number of substitution mutations occuring in the center position for all 5-mer 
motifs.


Usage
--------------------
```
createSubstitutionMatrix(
db,
model = c("S", "RS"),
sequenceColumn = "sequence_alignment",
germlineColumn = "germline_alignment_d_mask",
vCallColumn = "v_call",
multipleMutation = c("independent", "ignore"),
returnModel = c("5mer", "1mer", "1mer_raw"),
minNumMutations = 50,
numMutationsOnly = FALSE
)
```

Arguments
-------------------

db
:   data.frame containing sequence data.

model
:   type of model to create. The default model, "S", 
builds a model by counting only silent mutations. `model="S"`
should be used for data that includes functional sequences.
Setting `model="RS"` creates a model by counting both 
replacement and silent mutations and may be used on fully 
non-functional sequence data sets.

sequenceColumn
:   name of the column containing IMGT-gapped sample sequences.

germlineColumn
:   name of the column containing IMGT-gapped germline sequences.

vCallColumn
:   name of the column containing the V-segment allele call.

multipleMutation
:   string specifying how to handle multiple mutations occuring 
within the same 5-mer. If `"independent"` then multiple 
mutations within the same 5-mer are counted indepedently. 
If `"ignore"` then 5-mers with multiple mutations are 
excluded from the total mutation tally.

returnModel
:   string specifying what type of model to return; one of
`c("5mer", "1mer", "1mer_raw")`. If `"5mer"` 
(the default) then a 5-mer nucleotide context model is 
returned. If `"1mer"` or `"1mer_raw"` then a single 
nucleotide substitution matrix (no context) is returned;
where `"1mer_raw"` is the unnormalized version of the 
`"1mer"` model. Note, neither 1-mer model may be used
as input to [createMutabilityMatrix](createMutabilityMatrix.md).

minNumMutations
:   minimum number of mutations required to compute the 5-mer 
substitution rates. If the number of mutations for a 5-mer
is below this threshold, its substitution rates will be 
estimated from neighboring 5-mers. Default is 50. 
Not required if `numMutationsOnly=TRUE`.

numMutationsOnly
:   when `TRUE`, return counting information on the number
of mutations for each 5-mer, instead of building a substitution
matrix. This option can be used for parameter tuning for 
`minNumMutations` during preliminary analysis. 
Default is `FALSE`. Only applies when `returnModel` 
is set to `"5mer"`. The `data.frame` returned when
this argument is `TRUE` can serve as the input for
[minNumMutationsTune](minNumMutationsTune.md).




Value
-------------------

For `returnModel = "5mer"`: 

When `numMutationsOnly` is `FALSE`, a 4x1024 matrix of column 
normalized substitution rates for each 5-mer motif with row names defining 
the center nucleotide, one of `c("A", "C", "G", "T")`, and column names 
defining the 5-mer nucleotide sequence. 

When `numMutationsOnly` is 
`TRUE`, a 1024x4 data frame with each row providing information on 
counting the number of mutations for a 5-mer. Columns are named 
`fivemer.total`, `fivemer.every`, `inner3.total`, and
`inner3.every`, corresponding to, respectively,
the total number of mutations when counted as a 5-mer, 
whether there is mutation to every other base when counted as a 5-mer,
the total number of mutations when counted as an inner 3-mer, and
whether there is mutation to every other base when counted as an inner 3-mer.

For `returnModel = "1mer"` or `"1mer_raw"`:
a 4x4 normalized or un-normalized 1-mer substitution matrix respectively.


Details
-------------------

**Caution: The targeting model functions do NOT support ambiguous 
characters in their inputs. You MUST make sure that your input and germline
sequences do NOT contain ambiguous characters (especially if they are
clonal consensuses returned from `collapseClones`).**


References
-------------------


1. Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
on synonymous mutations from high-throughput immunoglobulin sequencing data. 
Front Immunol. 2013 4(November):358.
 



Examples
-------------------

```R
# Subset example data to one isotype and sample as a demo
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, c_call == "IGHA" & sample_id == "-1h")

# Count the number of mutations per 5-mer
subCount <- createSubstitutionMatrix(db, sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask",
vCallColumn="v_call",
model="S", multipleMutation="independent",
returnModel="5mer", numMutationsOnly=TRUE)

# Create model using only silent mutations
sub <- createSubstitutionMatrix(db, sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask",
vCallColumn="v_call",
model="S", multipleMutation="independent",
returnModel="5mer", numMutationsOnly=FALSE,
minNumMutations=20)
```



See also
-------------------

[extendSubstitutionMatrix](extendSubstitutionMatrix.md), [createMutabilityMatrix](createMutabilityMatrix.md), 
[createTargetingMatrix](createTargetingMatrix.md), [createTargetingModel](createTargetingModel.md),
[minNumMutationsTune](minNumMutationsTune.md).






