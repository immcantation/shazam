





**shmulateSeq** - *Simulate mutations in a single sequence*

Description
--------------------

Generates random mutations in a sequence iteratively using a targeting model.
Targeting probabilities at each position are updated after each iteration.


Usage
--------------------
```
shmulateSeq(sequence, mutations, targetingModel = HH_S5F)
```

Arguments
-------------------

sequence
:   sequence string in which mutations are to be introduced.

mutations
:   number of mutations to be introduced into `sequence`.

targetingModel
:   5-mer [TargetingModel](TargetingModel-class.md) object to be used for computing 
probabilities of mutations at each position. Defaults to
[HH_S5F](HH_S5F.md).




Value
-------------------

A string defining the mutated sequence.



Examples
-------------------

```R
# Define example input sequence
sequence <- "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA"

# Simulate using the default human 5-mer targeting model
shmulateSeq(sequence, mutations=6)
```


```
[1] "NCATATGACGACACGGCCGTCTATTTTTGTGCGAGAGATAGTTCA"

```



See also
-------------------

See [shmulateTree](shmulateTree.md) for imposing mutations on a lineage tree. 
See [HH_S5F](HH_S5F.md) and [MK_RS5NF](MK_RS5NF.md) for predefined 
[TargetingModel](TargetingModel-class.md) objects.



