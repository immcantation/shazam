





**calcTargetingDistance** - *Calculates a 5-mer distance matrix from a TargetingModel object*

Description
--------------------

`calcTargetingDistance` converts either the targeting rates in a `TargetingModel`
 model to a matrix of 5-mer to single-nucleotide mutation distances, or the substitution 
 rates in a 1-mer substitution model to a symmetric distance matrix.


Usage
--------------------
```
calcTargetingDistance(model, places = 2)
```

Arguments
-------------------

model
:   [TargetingModel](TargetingModel-class.md) object with mutation likelihood information, or
a 4x4 1-mer substitution matrix normalized by row with rownames and 
colnames consisting of "A", "T", "G", and "C".

places
:   decimal places to round distances to.




Value
-------------------

For input of [TargetingModel](TargetingModel-class.md), a matrix of distances for each 5-mer motif with 
rows names defining the center nucleotide and column names defining the 5-mer 
nucleotide sequence. For input of 1-mer substitution matrix, a 4x4 symmetric distance
matrix.


Details
-------------------

The targeting model is transformed into a distance matrix by:

1. Converting the likelihood of being mutated <code class = 'eq'>p=mutability*substitution</code> to 
distance <code class = 'eq'>d=-log10(p)</code>.
1. Dividing this distance by the mean of the distances.
1. Converting all infinite, no change (e.g., A->A), and NA distances to 
zero.


The 1-mer substitution matrix is transformed into a distance matrix by:

1. Symmetrize the 1-mer substitution matrix.
1. Converting the rates to distance <code class = 'eq'>d=-log10(p)</code>.
1. Dividing this distance by the mean of the distances.
1. Converting all infinite, no change (e.g., A -> A), and NA distances to 
zero.




Examples
-------------------

```R
# Calculate targeting distance of HH_S5F
dist <- calcTargetingDistance(HH_S5F)

# Calculate targeting distance of HH_S1F
dist <- calcTargetingDistance(HH_S1F)
```



See also
-------------------

See [TargetingModel](TargetingModel-class.md) for this class of objects and
[createTargetingModel](createTargetingModel.md) for building one.



