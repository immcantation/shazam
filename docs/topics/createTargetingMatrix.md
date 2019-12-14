**createTargetingMatrix** - *Calculates a targeting rate matrix*

Description
--------------------

`createTargetingMatrix` calculates the targeting model matrix as the
combined probability of mutability and substitution.


Usage
--------------------
```
createTargetingMatrix(substitutionModel, mutabilityModel)
```

Arguments
-------------------

substitutionModel
:   matrix of 5-mers substitution rates built by 
[createSubstitutionMatrix](createSubstitutionMatrix.md) or 
[extendSubstitutionMatrix](extendSubstitutionMatrix.md).

mutabilityModel
:   vector of 5-mers mutability rates built by 
[createMutabilityMatrix](createMutabilityMatrix.md) or 
[extendMutabilityMatrix](extendMutabilityMatrix.md).




Value
-------------------

A matrix with the same dimensions as the input `substitutionModel` 
containing normalized targeting probabilities for each 5-mer motif with 
row names defining the center nucleotide and column names defining the 
5-mer nucleotide sequence.


Details
-------------------

Targeting rates are calculated by multiplying the normalized mutability rate by the 
normalized substitution rates for each individual 5-mer.


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
db <- subset(ExampleDb, ISOTYPE == "IgA" & SAMPLE == "-1h")

# Create 4x1024 models using only silent mutations
sub_model <- createSubstitutionMatrix(db, model="S")
mut_model <- createMutabilityMatrix(db, sub_model, model="S")

```

*Warning*:Insufficient number of mutations to infer some 5-mers. Filled with 0. 
```R

# Extend substitution and mutability to including Ns (5x3125 model)
sub_model <- extendSubstitutionMatrix(sub_model)
mut_model <- extendMutabilityMatrix(mut_model)

# Create targeting model from substitution and mutability
tar_model <- createTargetingMatrix(sub_model, mut_model)
```



See also
-------------------

[createSubstitutionMatrix](createSubstitutionMatrix.md), [extendSubstitutionMatrix](extendSubstitutionMatrix.md), 
[createMutabilityMatrix](createMutabilityMatrix.md), [extendMutabilityMatrix](extendMutabilityMatrix.md), 
[createTargetingModel](createTargetingModel.md)






