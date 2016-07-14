





**shmulateSeq** - *Simulate mutations in a single sequence*

Description
--------------------

Simulate mutations in a single sequence


Usage
--------------------
```
shmulateSeq(input_seq, num_muts)
```

Arguments
-------------------

input_seq
:   sequence in which mutations are to be introduced

num_muts
:   number of mutations to be introduced into `input_seq`



Value
-------------------

A mutated sequence.

Details
-------------------

Generates mutations in sequence one by one while updating targeting
probability of each position after each mutation.





