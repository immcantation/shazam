





**extendSubstitutionMatrix** - *Extends a substitution model to include Ns.*

Description
--------------------

`extendSubstitutionMatrix` extends a 5-mer nucleotide substitution model 
with 5-mers that include Ns by averaging over all corresponding 5-mers without Ns.


Usage
--------------------
```
extendSubstitutionMatrix(substitutionModel)
```

Arguments
-------------------

substitutionModel
:   matrix of 5-mers substitution counts built by 
[createSubstitutionMatrix](createSubstitutionMatrix.md).



Value
-------------------

A 5x3125 matrix of normalized substitution rate for each 5-mer motif with 
rows names defining the center nucleotide, one of `c("A", "C", "G", "T", "N")`, 
and column names defining the 5-mer nucleotide sequence.



Examples
-------------------

```R
# Subset example data to one isotype and sample as a demo
db <- subset(InfluenzaDb, CPRIMER == "IGHA" & BARCODE == "RL014")

# Create model using only silent mutations and ignore multiple mutations
sub_model <- createSubstitutionMatrix(db, model="S", multipleMutation="ignore")
ext_model <- extendSubstitutionMatrix(sub_model)
```



See also
-------------------

[createSubstitutionMatrix](createSubstitutionMatrix.md), [extendMutabilityMatrix](extendMutabilityMatrix.md)



