**writeTargetingDistance** - *Write targeting model distances to a file*

Description
--------------------

`writeTargetingDistance` writes a 5-mer targeting distance matrix 
to a tab-delimited file.


Usage
--------------------
```
writeTargetingDistance(model, file)
```

Arguments
-------------------

model
:   [TargetingModel](TargetingModel-class.md) object with 
mutation likelihood information.

file
:   name of file to write.




Details
-------------------

The targeting distance write as a tab-delimited 5x3125 matrix. Rows define the mutated 
nucleotide at the center of each 5-mer, one of `c("A", "C", "G", "T", "N")`, 
and columns define the complete 5-mer of the unmutated nucleotide sequence. 
`NA` values in the distance matrix are replaced with distance 0.



Examples
-------------------

```R
### Not run:
# Write HS5F targeting model to working directory as hs5f.tab
# writeTargetingModel(HH_S5F, "hh_s5f.tsv")
```



See also
-------------------

Takes as input a [TargetingModel](TargetingModel-class.md) object and calculates  
distances using [calcTargetingDistance](calcTargetingDistance.md).



