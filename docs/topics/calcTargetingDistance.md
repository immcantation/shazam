





**calcTargetingDistance** - *Calculates a 5-mer distance matrix from a TargetingModel object*

Description
--------------------

`calcTargetingDistance` converts the targeting rates in a TargetingModel model 
to a matrix of 5-mer to single-nucleotide mutation distances.

Usage
--------------------

```
calcTargetingDistance(model)
```

Arguments
-------------------

model
:   [TargetingModel](TargetingModel-class.md) object with mutation likelihood information.



Value
-------------------

A matrix of distances for each 5-mer motif with rows names defining 
the center nucleotide and column names defining the 5-mer nucleotide 
sequence.

Details
-------------------

The targeting model is transformed into a distance matrix by:

1. Converting the likelihood of being mutated <code class = 'eq'>p=mutability*substitution</code> to 
distance <code class = 'eq'>d=-log10(p)</code>.
1. Dividing this distance by the mean of the distances
1. Converting all infinite, no change (e.g., A->A), and NA distances to 
zero.




Examples
-------------------

```R
# Calculate targeting distance of HS5FModel
dist <- calcTargetingDistance(HS5FModel)
```



See also
-------------------

Takes as input a [TargetingModel](TargetingModel-class.md) object. See [createTargetingModel](createTargetingModel.md)
for building a model.



