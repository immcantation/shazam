**shmulateSeq** - *Simulate mutations in a single sequence*

Description
--------------------

Generates random mutations in a sequence iteratively using a targeting model.
Targeting probabilities at each position are updated after each iteration.


Usage
--------------------
```
shmulateSeq(sequence, numMutations, targetingModel = HH_S5F)
```

Arguments
-------------------

sequence
:   sequence string in which mutations are to be introduced.
Accepted alphabet: `{A, T, G, C, N, .}`. Note
that `-` is not accepted.

numMutations
:   number of mutations to be introduced into `sequence`.

targetingModel
:   5-mer [TargetingModel](TargetingModel-class.md) object to be used for computing 
probabilities of mutations at each position. Defaults to
[HH_S5F](HH_S5F.md).




Value
-------------------

A string defining the mutated sequence.


Details
-------------------

If the input `sequence` has a non-triplet overhang at the end, it will be trimmed
to the last codon. For example, `ATGCATGC` will be trimmed to `ATGCAT`.

Mutations are not introduced to positions in the input `sequence` that contain 
`.` or `N`.



Examples
-------------------

```R
# Define example input sequence
sequence <- "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATA.TTTA"

# Simulate using the default human 5-mer targeting model
shmulateSeq(sequence, numMutations=6)
```


```
[1] "NGGTCTGACGACACGGCTGTCTATTATTGTGCGAGGGACA.TTTA"

```



See also
-------------------

See [shmulateTree](shmulateTree.md) for imposing mutations on a lineage tree. 
See [HH_S5F](HH_S5F.md) and [MK_RS5NF](MK_RS5NF.md) for predefined 
[TargetingModel](TargetingModel-class.md) objects.



