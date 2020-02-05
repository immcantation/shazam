**shmulateTree** - *Simulate mutations in a lineage tree*

Description
--------------------

`shmulateTree` returns a set of simulated sequences generated from an input 
sequence and a lineage tree. The input sequence is used to replace the most recent 
common ancestor (MRCA) node of the `igraph` object defining the lineage tree. 
Sequences are then simulated with mutations corresponding to edge weights in the tree. 
Sequences will not be generated for groups of nodes that are specified to be excluded.


Usage
--------------------
```
shmulateTree(
sequence,
graph,
targetingModel = HH_S5F,
field = NULL,
exclude = NULL,
junctionWeight = NULL,
start = 1,
end = nchar(sequence)
)
```

Arguments
-------------------

sequence
:   string defining the MRCA sequence to seed mutations from.

graph
:   `igraph` object defining the seed lineage tree, with 
vertex annotations, whose edges are to be recreated.

targetingModel
:   5-mer [TargetingModel](TargetingModel-class.md) object to be used for computing 
probabilities of mutations at each position. Defaults to
[HH_S5F](HH_S5F.md).

field
:   annotation to use for both unweighted path length exclusion 
and consideration as the MRCA node. If `NULL` do not 
exclude any nodes.

exclude
:   vector of annotation values in `field` to exclude from 
potential MRCA set. If `NULL` do not exclude any nodes.
Has no effect if `field=NULL`.

junctionWeight
:   fraction of the nucleotide sequence that is within the 
junction region. When specified this adds a proportional 
number of mutations to the immediate offspring nodes of the 
MRCA. Requires a value between 0 and 1. If `NULL` then 
edge weights are unmodified from the input `graph`.

start
:   Initial position in `sequence` where mutations can 
be introduced. Default: 1

end
:   Last position in `sequence` where mutations can 
be introduced. Default: last position (sequence length).




Value
-------------------

A `data.frame` of simulated sequences with columns:

+  `NAME`:      name of the corresponding node in the input 
`graph`.  
+  `SEQUENCE`:  mutated sequence.
+  `DISTANCE`:  Hamming distance of the mutated sequence from 
the seed `sequence`.




Examples
-------------------

```R
# Load example lineage and define example MRCA
data(ExampleTrees, package="alakazam")
graph <- ExampleTrees[[17]]
sequence <- "NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA"

# Simulate using the default human 5-mer targeting model
shmulateTree(sequence, graph)

```


```
            NAME                                      SEQUENCE DISTANCE
1      Inferred1 NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA        0
2 GN5SHBT07JDYW5 NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGAAAGTGTA        2
3 GN5SHBT03EP4KC NGATCTGACGACACGGCCGTATATTACTGTGCGAGAGAAAGTACA        3
4 GN5SHBT01AKANC NGATCCGACGACACGGCCGTCTATTATTGTGTGAGAGATAGTTTA        4
5 GN5SHBT01A3SFZ NGATCTGACGACACGGCCGTGTCCTACTGTGTCAGAGAAAGTTCA        6
6 GN5SHBT08HUU7M NGATCTGACGACACGGCCCTATATTACTTTGCGAGAGAAAGTTCA        3
7 GN5SHBT04CEA6I NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGAACGTGTA        1
8 GN5SHBT06IXJIH NGATCTGACGACACGGCCGTGAATTACTGTGCGAGAGAAAGTGTA        1

```


```R

# Simulate using the mouse 5-mer targeting model
# Exclude nodes without a sample identifier
# Add 20% mutation rate to the immediate offsprings of the MRCA
shmulateTree(sequence, graph, targetingModel=MK_RS5NF,
field="SAMPLE", exclude=NA, junctionWeight=0.2)
```


```
            NAME                                      SEQUENCE DISTANCE
1 GN5SHBT07JDYW5 NGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA        0
2 GN5SHBT03EP4KC NGATCTGTCGATACGGCCGTATGTTACTGTGCGAGAGATAGTTTA        4
3 GN5SHBT01A3SFZ NGGTCTGGCGACACGGTCGTATATGACTGTGTGAGATATAGTTTA        7
4 GN5SHBT08HUU7M NGATCTGTCGATACGGCCGTGTATTACTGTGCGAGAGATAATTTA        3
5 GN5SHBT04CEA6I NGATCTGACGACACGGCCGTGTATTCCTGTGCGAGAGATAGTTTA        1
6 GN5SHBT06IXJIH NGATCTGATGACACGGCCGTGTATTACTGTGCGAGAGATAGTTTA        1

```



See also
-------------------

See [shmulateSeq](shmulateSeq.md) for imposing mutations on a single sequence. 
See [HH_S5F](HH_S5F.md) and [MK_RS5NF](MK_RS5NF.md) for predefined 
[TargetingModel](TargetingModel-class.md) objects.






