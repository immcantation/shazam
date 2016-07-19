





**shmulateSeq** - *Simulate mutations in a single sequence*

Description
--------------------

Generates mutations in sequence one by one while updating targeting
probability of each position after each mutation.


Usage
--------------------
```
shmulateSeq(input_seq, num_muts, targeting_model = HS5FModel)
```

Arguments
-------------------

input_seq
:   sequence in which mutations are to be introduced

num_muts
:   number of mutations to be introduced into `input_seq`

targeting_model
:   targeting model of class `TargetingModel` to be used for 
computing probabilities of mutations at each position. Default is
`HS5FModel` from `SHazaM`.



Value
-------------------

A mutated sequence.



Examples
-------------------

```R
# Example input sequence
input_seq <- "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA"

# Simulate using the default human S5F targeting model
shmulateSeq(input_seq, num_muts = 6)
```


```
[1] "NGGTCTGACTCCACGGCCGTGTATCATTGTGCGAGAGATAGTATA"

```



See also
-------------------

[shmulateTree](shmulateTree.md), [HS5FModel](HS5FModel.md), [TargetingModel](TargetingModel-class.md)



