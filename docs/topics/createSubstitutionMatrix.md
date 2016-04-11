





**createSubstitutionMatrix** - *Builds a substitution model*

Description
--------------------

`createSubstitutionMatrix` builds a 5-mer nucleotide substitution model by counting 
the number of substitution mutations occuring in the center position for all 5-mer 
motifs.

Usage
--------------------

```
createSubstitutionMatrix(db, model = c("RS", "S"),
sequenceColumn = "SEQUENCE_IMGT", germlineColumn = "GERMLINE_IMGT_D_MASK",
vCallColumn = "V_CALL", multipleMutation = c("independent", "ignore"),
returnModel = c("5mer", "1mer", "1mer_raw"), minNumMutations = 50)
```

Arguments
-------------------

db
:   data.frame containing sequence data.

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



Value
-------------------

A 4x1024 matrix of column normalized substitution rates for each 5-mer motif with 
row names defining the center nucleotide, one of `c("A", "C", "G", "T")`, 
and column names defining the 5-mer nucleotide sequence.

References
-------------------


1. Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
on synonymous mutations from high-throughput immunoglobulin sequencing data. 
Front Immunol. 2013 4(November):358.
 



Examples
-------------------

```R
# Subset example data to one isotype and sample as a demo
db <- subset(InfluenzaDb, CPRIMER == "IGHA" & BARCODE == "RL014")

# Create model using only silent mutations and ignore multiple mutations
sub <- createSubstitutionMatrix(db, model="S", multipleMutation="ignore")
```



See also
-------------------

[extendSubstitutionMatrix](extendSubstitutionMatrix.md), [createMutabilityMatrix](createMutabilityMatrix.md), 
[createTargetingMatrix](createTargetingMatrix.md), [createTargetingModel](createTargetingModel.md)



