





**minNumSeqMutationsTune** - *Parameter tuning for minNumSeqMutations*

Description
--------------------

`minNumSeqMutationsTune` helps with picking a threshold value for `minNumSeqMutations`
in [createMutabilityMatrix](createMutabilityMatrix.md) by tabulating the number of 5-mers for which 
mutability would be computed directly or inferred at various threshold values.


Usage
--------------------
```
minNumSeqMutationsTune(mutCount, minNumSeqMutationsRange)
```

Arguments
-------------------

mutCount
:   a `vector` of length 1024 returned by 
[createMutabilityMatrix](createMutabilityMatrix.md) with `numSeqMutationsOnly=TRUE`.

minNumSeqMutationsRange
:   a number or a vector indicating the value or the range of values 
of `minNumSeqMutations` to try.



Value
-------------------

A 2xn `matrix`, where n is the number of trial values of `minNumSeqMutations`
supplied in `minNumSeqMutationsRange`. Each column corresponds to a value
in `minNumSeqMutationsRange`. The rows correspond to the number of 5-mers
for which mutability would be computed directly (`"measured"`) and inferred
(`"inferred"`), respectively.

Details
-------------------

At a given threshold value of `minNumSeqMutations`, for a given 5-mer,
if the total number of mutations is greater than the threshold, mutability 
is computed directly. Otherwise, mutability is inferred.

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

# Create model using only silent mutations
sub <- createSubstitutionMatrix(db, model="S", multipleMutation="independent",
returnModel="5mer", numMutationsOnly=FALSE,
minNumMutations=20)

# Count the number of mutations in sequences containing each 5-mer
mutCount <- createMutabilityMatrix(db, substitutionModel = sub,
model="S", multipleMutation="independent",
numSeqMutationsOnly=TRUE)

# Tune minNumSeqMutations
minNumSeqMutationsTune(mutCount, seq(from=100, to=300, by=50))
```


```
         100 150 200 250 300
measured 237 200 140 107  97
inferred 787 824 884 917 927

```



See also
-------------------

See argument `numSeqMutationsOnly` in [createMutabilityMatrix](createMutabilityMatrix.md) 
for generating the required input `vector` `mutCount`. 
See argument `minNumSeqMutations` in [createMutabilityMatrix](createMutabilityMatrix.md)
for what it does.



