





**extendMutabilityMatrix** - *Extends a mutability model to include Ns.*

Description
--------------------

`extendMutabilityMatrix` extends a 5-mer nucleotide mutability model 
with 5-mers that include Ns by averaging over all corresponding 5-mers without Ns.

Usage
--------------------

```
extendMutabilityMatrix(mutabilityModel)
```

Arguments
-------------------

mutabilityModel
:   vector of 5-mer mutability rates built by 
[createMutabilityMatrix](createMutabilityMatrix.md).



Value
-------------------

A 3125 vector of normalized mutability rates for each 5-mer motif with 
names defining the 5-mer nucleotide sequence.



Examples
-------------------

```R
# Subset example data to one isotype and sample as a demo
db <- subset(InfluenzaDb, CPRIMER == "IGHA" & BARCODE == "RL014")

# Create model using only silent mutations and ignore multiple mutations
sub_model <- createSubstitutionMatrix(db, model="S", multipleMutation="ignore")
mut_model <- createMutabilityMatrix(db, sub_model, model="S", multipleMutation="ignore",
minNumSeqMutations=10)

```

*Warning*:Insufficient number of mutations to infer some 5-mers. Filled with 0. 
```R
ext_model <- extendMutabilityMatrix(mut_model)
```



See also
-------------------

[createMutabilityMatrix](createMutabilityMatrix.md), [extendSubstitutionMatrix](extendSubstitutionMatrix.md)



