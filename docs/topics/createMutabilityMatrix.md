





**createMutabilityMatrix** - *Builds a mutability model*

Description
--------------------

`createMutabilityMatrix` builds a 5-mer nucleotide mutability model by counting 
the number of mutations occuring in the center position for all 5-mer motifs.


Usage
--------------------
```
createMutabilityMatrix(db, substitutionModel, model = c("RS", "S"),
sequenceColumn = "SEQUENCE_IMGT", germlineColumn = "GERMLINE_IMGT_D_MASK",
vCallColumn = "V_CALL", multipleMutation = c("independent", "ignore"),
minNumSeqMutations = 500, returnSource = FALSE)
```

Arguments
-------------------

db
:   data.frame containing sequence data.

substitutionModel
:   matrix of 5-mer substitution rates built by 
[createSubstitutionMatrix](createSubstitutionMatrix.md).

model
:   type of model to create. The default model, "RS", creates 
a model by counting both replacement and silent mutations.
The "S" specification builds a model by counting only 
silent mutations.

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

minNumSeqMutations
:   minimum number of mutations in sequences containing each 5-mer
to compute the mutability rates. If the number is smaller 
than this threshold, the mutability for the 5-mer will be 
inferred. Default is 500.

returnSource
:   return the sources of 5-mer mutabilities (measured vs.
inferred). Default is false.



Value
-------------------

A named numeric vector of 1024 normalized mutability rates for each 5-mer 
motif with names defining the 5-mer nucleotide sequence.

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
db <- subset(ExampleDb, ISOTYPE == "IgG" & SAMPLE == "+7d")

# Create model using only silent mutations and ignore multiple mutations
sub_model <- createSubstitutionMatrix(db, model="S", multipleMutation="ignore")
mut_model <- createMutabilityMatrix(db, sub_model, model="S", 
multipleMutation="ignore", minNumSeqMutations=10)
```

*Warning*:Insufficient number of mutations to infer some 5-mers. Filled with 0. 

See also
-------------------

[extendMutabilityMatrix](extendMutabilityMatrix.md), [createSubstitutionMatrix](createSubstitutionMatrix.md), 
[createTargetingMatrix](createTargetingMatrix.md), [createTargetingModel](createTargetingModel.md)



