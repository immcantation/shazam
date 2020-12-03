**plotCloneLineageTree** - *Plotting a tree-plot of a specific clone:*

Description
--------------------

Plotting a tree-plot of a specific clone:


Usage
--------------------
```
plotCloneLineageTree(
db,
curClone = NULL,
id = "sequence_id",
seq = "sequence_alignment",
germ = "germline_alignment",
v_call = "v_call",
j_call = "j_call",
junc_len = "junction_length",
clone = "clone_id",
mask_char = "N",
max_mask = 0,
pad_end = FALSE,
dnapars_exec
)
```

Arguments
-------------------

db
:   a `ChangeoClone` database

curClone
:   clone number for which to plot the tree plot.

id
:   name of the column containing sequence identifiers.

seq
:   name of the column containing observed DNA sequences. 
All sequences in this column must be multiple aligned.

germ
:   name of the column containing germline DNA sequences. 
All entries in this column should be identical for any given clone, 
and they must be multiple aligned with the data in the seq column.

v_call
:   name of the column containing V-segment allele assignments. 
All entries in this column should be identical to the gene level.

j_call
:   name of the column containing J-segment allele assignments. 
All entries in this column should be identical to the gene level.

junc_len
:   name of the column containing the length of the junction as a numeric value. 
All entries in this column should be identical for any given clone.

clone
:   name of the column containing the identifier for the clone. 
All entries in this column should be identical.

mask_char
:   character to use for masking and padding.

max_mask
:   maximum number of characters to mask at the leading and trailing sequence ends.
If NULL then the upper masking bound will be automatically determined from the
maximum number of observed leading or trailing Ns amongst all sequences. 
If set to 0 (default) then masking will not be performed

pad_end
:   if TRUE pad the end of each sequence with mask_char to make every sequence the same length.

dnapars_exec
:   absolute path to the PHYLIP dnapars executable




Value
-------------------

a plot of the specific clone linegae tree.


Details
-------------------

This function will take as input any `ChangeoClone` db, and a specific
clone number. It will build a `ChangeoClone` object and an [igraph](http://www.rdocumentation.org/packages/igraph/topics/aaa-igraph-package)
object for the specific clone, and will plot a tree plot for it.

Notes: 

1. This function will give an error in case that all the sequences in the 
db of the specific clone are the same (as no lineage tree can be built on that).

2. In case one of the sequences is equal to the germline sequence - then the 
graph may show that sequence ID as the graph root (and not the "Germline" sequence ID).



Examples
-------------------

```R
library ("igraph")
library("alakazam")
data(ExampleDb)
dnapars_exec <- "~/apps/phylip-3.69/dnapars"
plotCloneLineageTree(db=ExampleDb, curClone=3139, seq="sequence_alignment",
id="sequence_id", germ="germline_alignment", dnapars_exec = dnapars_exec)
```

**Error in buildPhylipLineage(curCloneObj, dnapars_exec, rm_temp = TRUE)**: The file ~/apps/phylip-3.69/dnapars cannot be executed.






