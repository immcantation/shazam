**createMutabilityMatrix** - *Builds a mutability model*

Description
--------------------

`createMutabilityMatrix` builds a 5-mer nucleotide mutability model by counting 
the number of mutations occurring in the center position for all 5-mer motifs.


Usage
--------------------
```
createMutabilityMatrix(
db,
substitutionModel,
model = c("s", "rs"),
sequenceColumn = "sequence_alignment",
germlineColumn = "germline_alignment_d_mask",
vCallColumn = "v_call",
multipleMutation = c("independent", "ignore"),
minNumSeqMutations = 500,
numSeqMutationsOnly = FALSE
)
```

Arguments
-------------------

db
:   data.frame containing sequence data.

substitutionModel
:   matrix of 5-mer substitution rates built by 
[createSubstitutionMatrix](createSubstitutionMatrix.md). Note, this model will
only impact mutability scores when `model="s"`
(using only silent mutations).

model
:   type of model to create. The default model, "s", 
builds a model by counting only silent mutations. `model="s"`
should be used for data that includes functional sequences.
Setting `model="rs"` creates a model by counting both 
replacement and silent mutations and may be used on fully 
non-functional sequence data sets.

sequenceColumn
:   name of the column containing IMGT-gapped sample sequences.

germlineColumn
:   name of the column containing IMGT-gapped germline sequences.

vCallColumn
:   name of the column containing the V-segment allele call.

multipleMutation
:   string specifying how to handle multiple mutations occurring 
within the same 5-mer. If `"independent"` then multiple 
mutations within the same 5-mer are counted independently. 
If `"ignore"` then 5-mers with multiple mutations are 
excluded from the total mutation tally.

minNumSeqMutations
:   minimum number of mutations in sequences containing each 5-mer
to compute the mutability rates. If the number is smaller 
than this threshold, the mutability for the 5-mer will be 
inferred. Default is 500. Not required if 
`numSeqMutationsOnly=TRUE`.

numSeqMutationsOnly
:   when `TRUE`, return only a vector counting the number of 
observed mutations in sequences containing each 5-mer. This 
option can be used for parameter tuning for `minNumSeqMutations` 
during preliminary analysis using [minNumSeqMutationsTune](minNumSeqMutationsTune.md). 
Default is `FALSE`.




Value
-------------------

When `numSeqMutationsOnly` is `FALSE`, a `MutabilityModel` containing a
named numeric vector of 1024 normalized mutability rates for each 5-mer motif with names 
defining the 5-mer nucleotide sequence.

When `numSeqMutationsOnly` is `TRUE`, a named numeric
vector of length 1024 counting the number of observed mutations in sequences containing 
each 5-mer.


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
# Subset example data to 50 sequences of one isotype and sample as a demo
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, c_call == "IGHA" & sample_id == "-1h")[1:50,]

# Create model using only silent mutations
sub_model <- createSubstitutionMatrix(db, sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask",
vCallColumn="v_call",model="s")
mut_model <- createMutabilityMatrix(db, sub_model, model="s", 
sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask",
vCallColumn="v_call",
minNumSeqMutations=200,
numSeqMutationsOnly=FALSE)

```

*Warning*:Insufficient number of mutations to infer some 5-mers. Filled with 0. 
```R

# View top 5 mutability estimates
head(sort(mut_model, decreasing=TRUE), 5)

```


```
    AAGTA     ACGTA     AGGTA     ATGTA     CAGTA 
0.0088995 0.0088995 0.0088995 0.0088995 0.0088995 

```


```R

# View the number of S mutations used for estimating mutabilities
mut_model@numMutS

```


```
[1] 410

```


```R

# Count the number of mutations in sequences containing each 5-mer
mut_count <- createMutabilityMatrix(db, sub_model, model="s", 
sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask",
vCallColumn="v_call",
numSeqMutationsOnly=TRUE)

```



See also
-------------------

[MutabilityModel](MutabilityModel-class.md), [extendMutabilityMatrix](extendMutabilityMatrix.md), [createSubstitutionMatrix](createSubstitutionMatrix.md), 
[createTargetingMatrix](createTargetingMatrix.md), [createTargetingModel](createTargetingModel.md),
[minNumSeqMutationsTune](minNumSeqMutationsTune.md)






