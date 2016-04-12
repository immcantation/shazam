





**collapseByClone** - *Identifies clonal consensus sequences*

Description
--------------------

Identifies effective/consensus sequences collapsed by clone


Usage
--------------------
```
collapseByClone(db, cloneColumn = "CLONE", sequenceColumn = "SEQUENCE_IMGT",
germlineColumn = "GERMLINE_IMGT_D_MASK", expandedDb = FALSE,
regionDefinition = NULL, nonTerminalOnly = FALSE, nproc = 1)
```

Arguments
-------------------

db
:   `data.frame` containing sequence data.

cloneColumn
:   `character` name of the column containing clonal 
identifiers.

sequenceColumn
:   `character` name of the column containing input 
sequences.

germlineColumn
:   `character` name of the column containing germline 
sequences.

expandedDb
:   `logical` indicating whether or not to return the 
expanded `db`, containing all the sequences (as opposed
to returning just one sequence per clone collapsed by )

regionDefinition
:   [RegionDefinition](RegionDefinition-class.md) object defining the regions
and boundaries of the Ig sequences.

nonTerminalOnly
:   `logical` indicating whether to include mutations
at the leaves.

nproc
:   Number of cores to distribute the operation over. If the 
`cluster` has already been set earlier, then pass the 
`cluster`. This will ensure that it is not reset.



Value
-------------------

A modified `db` data.frame with clonal consensus sequences in the
CLONAL_SEQUENCE column.

Details
-------------------

`collapseByClone` identifies the consensus sequence of each clonal 
group and appends a column to the input `data.frame` containing the clonal 
consensus for each sequence.


For sequences identified to be part of the same clone, this function defines an 
effective sequence that will be representative for all mutations in the clone. Each 
position in this consensus (or effective) sequence is created by a weighted sampling 
of each mutated base (and non "N", "." or "-" characters) from all the sequences in 
the clone. 

For example, in a clone with 5 sequences that have a C at position 1, and 5 sequences
with a T at this same position, the consensus sequence will have a C 50%  and T 50% 
of the time it is called.

The function returns an updated `db` that collpases all the sequences by clones 
defined in the `cloneColumn` column argument.

Non-terminal branch mutations are defined as the set of mutations that occur on 
branches of the lineage tree that are not connected to a leaf. For computational 
efficiency, the set of non-terminal branch mutations is approximated as those that are
shared between more than one sequence in a clone. In this case the terminal branch 
mutations are filtered out.

This function can be parallelized if `db` contains thousands of sequences. 
Specify the number of cores available using the `nproc` parameter.



Examples
-------------------

```R
# Subset example data
db <- subset(InfluenzaDb, CPRIMER %in% c("IGHA","IGHM") & 
BARCODE %in% c("RL016","RL018","RL019","RL021"))

# Run collapseByClone
db_new <- collapseByClone(db, cloneColumn="CLONE", 
sequenceColumn="SEQUENCE_IMGT",
germlineColumn="GERMLINE_IMGT_D_MASK",
expandedDb=FALSE,
regionDefinition=IMGT_V_NO_CDR3,
nproc=1)
```


```
Collapsing clonal sequences...

```



See also
-------------------

See [IMGT_SCHEMES](IMGT_SCHEMES.md) for a set of predefined [RegionDefinition](RegionDefinition-class.md) objects.



