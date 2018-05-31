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
names defining the 5-mer nucleotide sequence. Note that "normalized" means
that the mutability rates for the 1024 5-mers that contain no "N" at any
position sums up to 1 (as opposed to the entire vector summing up to 1).



Examples
-------------------

```R
# Subset example data to one isotype and sample as a demo
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, ISOTYPE == "IgA" & SAMPLE == "-1h")

# Create model using only silent mutations and ignore multiple mutations
sub_model <- createSubstitutionMatrix(db, model="S")
mut_model <- createMutabilityMatrix(db, sub_model, model="S")

```

*Warning*:Insufficient number of mutations to infer some 5-mers. Filled with 0. 
```R
ext_model <- extendMutabilityMatrix(mut_model)
```



See also
-------------------

[createMutabilityMatrix](createMutabilityMatrix.md), [extendSubstitutionMatrix](extendSubstitutionMatrix.md)



