**createTargetingModel** - *Creates a TargetingModel*

Description
--------------------

`createTargetingModel` creates a 5-mer `TargetingModel`.


Usage
--------------------
```
createTargetingModel(
db,
model = c("S", "RS"),
sequenceColumn = "sequence_alignment",
germlineColumn = "germline_alignment_d_mask",
vCallColumn = "v_call",
multipleMutation = c("independent", "ignore"),
minNumMutations = 50,
minNumSeqMutations = 500,
modelName = "",
modelDescription = "",
modelSpecies = "",
modelCitation = "",
modelDate = NULL
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
:   name of the column containing the V-segment allele calls.

multipleMutation
:   string specifying how to handle multiple mutations occuring 
within the same 5-mer. If `"independent"` then multiple 
mutations within the same 5-mer are counted indepedently. 
If `"ignore"` then 5-mers with multiple mutations are 
excluded from the otal mutation tally.

minNumMutations
:   minimum number of mutations required to compute the 5-mer 
substitution rates. If the number of mutations for a 5-mer
is below this threshold, its substitution rates will be 
estimated from neighboring 5-mers. Default is 50.

minNumSeqMutations
:   minimum number of mutations in sequences containing each 5-mer
to compute the mutability rates. If the number is smaller 
than this threshold, the mutability for the 5-mer will be 
inferred. Default is 500.

modelName
:   name of the model.

modelDescription
:   description of the model and its source data.

modelSpecies
:   genus and species of the source sequencing data.

modelCitation
:   publication source.

modelDate
:   date the model was built. If `NULL` the current date
will be used.




Value
-------------------

A [TargetingModel](TargetingModel-class.md) object.


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

# Create model using only silent mutations and ignore multiple mutations
model <- createTargetingModel(db, model="S", sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask",
vCallColumn="v_call", multipleMutation="ignore")

```

*Warning*:Insufficient number of mutations to infer some 5-mers. Filled with 0. 
```R


# Access and view mutability estimates (not run)
# print(model@mutability)

# View the number of S mutations used for estimating mutabilities
model@mutability@numMutS
```


```
[1] 0

```



See also
-------------------

See [TargetingModel](TargetingModel-class.md) for the return object. 
See [plotMutability](plotMutability.md) plotting a mutability model.
See [createSubstitutionMatrix](createSubstitutionMatrix.md), [extendSubstitutionMatrix](extendSubstitutionMatrix.md), 
[createMutabilityMatrix](createMutabilityMatrix.md), [extendMutabilityMatrix](extendMutabilityMatrix.md) and 
[createTargetingMatrix](createTargetingMatrix.md) for component steps in building a model.






