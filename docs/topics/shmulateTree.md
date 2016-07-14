





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
shmulateTree(input_seq, graph, field = NULL, exclude = NULL,
jun_frac = NULL)
```

Arguments
-------------------

input_seq
:   sequence in which mutations are to be introduced.

graph
:   `igraph` object with vertex annotations whose edges are to be recreated.

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





