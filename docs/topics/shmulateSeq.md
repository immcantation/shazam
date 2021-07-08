**shmulateSeq** - *Simulate mutations in a single sequence*

Description
--------------------

Generates random mutations in a sequence iteratively using a targeting model.
Targeting probabilities at each position are updated after each iteration.


Usage
--------------------
```
shmulateSeq(
sequence,
numMutations,
targetingModel = HH_S5F,
start = 1,
end = nchar(sequence),
frequency = FALSE
)
```

Arguments
-------------------

sequence
:   sequence string in which mutations are to be introduced.
Accepted alphabet: `{A, T, G, C, N, .}`. Note
that `-` is not accepted.

numMutations
:   a whole number indicating the number of mutations to be 
introduced into `sequence`, if `frequency=FALSE`.
A fraction bewteen 0 and 1 indicating the mutation frequency
if `frequency=TRUE`.

targetingModel
:   5-mer [TargetingModel](TargetingModel-class.md) object to be used for computing 
probabilities of mutations at each position. Defaults to
[HH_S5F](HH_S5F.md).

start
:   Initial position in `sequence` where mutations can 
be introduced. Default: 1

end
:   Last position in `sequence` where mutations can 
be introduced. Default: last position (sequence length).

frequency
:   If `TRUE`, treat `numMutations` as a frequency.




Value
-------------------

A string defining the mutated sequence.


Details
-------------------

If the input `sequence` has a non-triplet overhang at the end, it will be trimmed
to the last codon. For example, `ATGCATGC` will be trimmed to `ATGCAT`.

Mutations are not introduced to positions in the input `sequence` that contain 
`.` or `N`.

With `frequency=TRUE`, the number of mutations introduced is the `floor` of 
the length of the sequence multiplied by the mutation frequency specified via
`numMutations`.



Examples
-------------------

```R
# Define example input sequence
sequence <- "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATA.TTTA"

# Simulate using the default human 5-mer targeting model
# Introduce 6 mutations
shmulateSeq(sequence, numMutations=6, frequency=FALSE)

```


```
[1] "NCATCTGACGTCACGGCCGAGTATTACTATGCGAGAGACA.TTCA"

```


```R

# Introduction 5% mutations
shmulateSeq(sequence, numMutations=0.05, frequency=TRUE)
```


```
[1] "NGATCTGACGGCACGGCCGTGTATTACTGTGCGAGTGATA.TTTA"

```



See also
-------------------

See [shmulateTree](shmulateTree.md) for imposing mutations on a lineage tree. 
See [HH_S5F](HH_S5F.md) and [MK_RS5NF](MK_RS5NF.md) for predefined 
[TargetingModel](TargetingModel-class.md) objects.






