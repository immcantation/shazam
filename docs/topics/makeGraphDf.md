**makeGraphDf** - *Generating a data frame from graph object, and mering it with clone data frame*

Description
--------------------

Generating a data frame from graph object, and mering it with clone data frame


Usage
--------------------
```
makeGraphDf(
curCloneGraph,
curCloneObj,
objSeqId = "sequence_id",
objSeq = "sequence"
)
```

Arguments
-------------------

curCloneGraph
:   an [igraph](http://www.rdocumentation.org/packages/igraph/topics/aaa-igraph-package) object of the specific clone.

curCloneObj
:   `ChangeoClone` object of the specific clone.

objSeqId
:   sequence id field name of `curCloneObj`

objSeq
:   sequence field name of `curCloneObj`




Value
-------------------

A `ChangeoClone` object with additional columns (parent_sequence and parent)
and additional rows (for germline and inferred sequences)


Details
-------------------

`makeGraphDf` adds columns and rows to the clones database: 

Additional **columns** are added for parent_sequence and parent 
(which is the parent sequence id).

Additional **rows** are added for inferred sequences and the germline of the clone graph.

[makeGraphDf](makeGraphDf.md) also renames sequence_id content according to the following 
(assume clone number is 34):  

34_Germline, 34_Inferred1, 34_1, 34_2, 34_3, 34_Inferred2, 34_4, etc.

The original sequence id is kept under a new column named `orig_sequence_id`, 
and the original parent sequence id is kept under a new column named `orig_parent`.

Note: sequence field name of `curCloneGraph` argument must be "sequence".



Examples
-------------------

```R
### Not run:
library("igraph")

```

*
Attaching package: ‘igraph’
**The following objects are masked from ‘package:stats’:

    decompose, spectrum
**The following object is masked from ‘package:base’:

    union
*
```R
# library("dplyr")
# # Load and subset example data:
# data(ExampleDb, package="alakazam")
# clone_3170_db <- subset(ExampleDb, clone_id == 3170)
# clone_3170_obj <- makeChangeoClone(clone_3170_db, 
# seq="sequence_alignment",
# germ="germline_alignment")
# dnapars_exec <- "~/apps/phylip-3.69/dnapars"
# clone_3170_graph <- buildPhylipLineage(clone_3170_obj, 
# dnapars_exec, rm_temp = TRUE)  
# clone_3170_GraphDf <- makeGraphDf(clone_3170_graph, clone_3170_obj)
```








