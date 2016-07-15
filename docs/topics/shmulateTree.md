





**shmulateTree** - *Simulate sequences to populate a tree*

Description
--------------------

`shmulateTree` returns a set of simulated sequences generated from an input sequence and an
igraph object. The input sequence is used to replace the founder node of the `igraph` lineage
tree and sequences are simulated with mutations corresponding to edge weights in the tree.
Sequences will not be generated for groups of nodes that are specified to be excluded.


Usage
--------------------
```
shmulateTree(input_seq, graph, targeting_model = HS5FModel, field = NULL,
exclude = NULL, jun_frac = NULL)
```

Arguments
-------------------

input_seq
:   sequence in which mutations are to be introduced.

graph
:   `igraph` object with vertex annotations whose edges are to be recreated.

targeting_model
:   targeting model of class `TargetingModel` to be used for 
computing probabilities of mutations at each position. Default is
`HS5FModel` from `SHazaM`.

field
:   annotation field to use for both unweighted path length exclusion and
consideration as a founder node. If `NULL` do not exclude any nodes.

exclude
:   vector of annotation values in the given field to exclude from potential
founder set. If `NULL` do not exclude any nodes. Has no effect if `field=NULL`.

jun_frac
:   fraction of characters in the junction region to add proportional number
of trunk mutations to the sequence.



Value
-------------------

A `data.frame` of simulated sequences.



Examples
-------------------

```R
# Example input
input_seq <- "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA"

# Load example graph
library(alakazam)
graph <- ExampleTrees[[17]]

# Simulate using the mouse RS5NF targeting model
shmulateTree(input_seq, graph, targeting_model = MRS5NFModel)
```


```
            name                                      sequence distance
1      Inferred1 NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA        0
2 GN5SHBT07JDYW5 NGATCTGACGACACGGCTGTGTATTACTGTGCGAAAGATAGTTTA        2
3 GN5SHBT01AKANC NGTTCTGACGACACAGTCGTGTATTACTGTGCGAGAGATAGTTTG        4
4 GN5SHBT03EP4KC NGATCTGACGACACGGCTGTGTGTTACTGGGCGAAATATAGTTTA        3
5 GN5SHBT01A3SFZ NGATCTGACGACACGACTATATATTACTTTACGAGAGATAGTTTA        6
6 GN5SHBT04CEA6I NGATCTGACGACACGGCTGTATATTACTGTGCGAAAGATAGTTTA        1
7 GN5SHBT06IXJIH NGATCTGACGACACGGCTGTGTATTACTGTGCGAAAAATAGTTTA        1
8 GN5SHBT08HUU7M NGATCTGACGATACGATTGTGTGTTACTGGGCGAAATATAGTTTA        3

```



See also
-------------------

[shmulateSeq](shmulateSeq.md), [HS5FModel](HS5FModel.md), [TargetingModel](TargetingModel-class.md), [ExampleTrees](http://www.inside-r.org/packages/cran/alakazam/docs/ExampleTrees)



