**minNumMutationsTune** - *Parameter tuning for minNumMutations*

Description
--------------------

`minNumMutationsTune` helps with picking a threshold value for `minNumMutations`
in [createSubstitutionMatrix](createSubstitutionMatrix.md) by tabulating the number of 5-mers for which 
substitution rates would be computed directly or inferred at various threshold values.


Usage
--------------------
```
minNumMutationsTune(subCount, minNumMutationsRange)
```

Arguments
-------------------

subCount
:   `data.frame` returned by [createSubstitutionMatrix](createSubstitutionMatrix.md)
with `numMutationsOnly=TRUE`.

minNumMutationsRange
:   a number or a vector indicating the value or range of values
of `minNumMutations` to try.




Value
-------------------

A 3xn `matrix`, where n is the number of trial values of `minNumMutations`
supplied in `minNumMutationsRange`. Each column corresponds to a value
in `minNumMutationsRange`. The rows correspond to the number of 5-mers
for which substitution rates would be computed directly using the 5-mer itself 
(`"5mer"`), using its inner 3-mer (`"3mer"`), and using the central 
1-mer (`"1mer"`), respectively.


Details
-------------------

At a given threshold value of `minNumMutations`, for a given 5-mer,
if the total number of mutations is greater than the threshold and there
are mutations to every other base, substitution rates are computed directly
for the 5-mer using its mutations. Otherwise, mutations from 5-mers with 
the same inner 3-mer as the 5-mer of interest are aggregated. If the number 
of such mutations is greater than the threshold and there are mutations to 
every other base, these mutations are used for inferring the substitution 
rates for the 5-mer of interest; if not, mutations from all 5-mers with the 
same center nucleotide are aggregated and used for inferring the substitution
rates for the 5-mer of interest (i.e. the 1-mer model).


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
db <- subset(ExampleDb, c_call == "IGHA" & sample_id == "-1h")

# Count the number of mutations per 5-mer
subCount <- createSubstitutionMatrix(db, sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask",
vCallColumn="v_call",
model="s", multipleMutation="independent",
returnModel="5mer", numMutationsOnly=TRUE)

# Tune minNumMutations
minNumMutationsTune(subCount, seq(from=10, to=80, by=10))

```


```
      10  20  30   40   50   60   70   80
5mer   5   2   1    1    0    0    0    0
3mer 143 111  53   22   17    2    1    1
1mer 876 911 970 1001 1007 1022 1023 1023

```



See also
-------------------

See argument `numMutationsOnly` in [createSubstitutionMatrix](createSubstitutionMatrix.md) 
for generating the required input `data.frame` `subCount`. 
See argument `minNumMutations` in [createSubstitutionMatrix](createSubstitutionMatrix.md)
for what it does.






