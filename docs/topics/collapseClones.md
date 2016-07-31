





**collapseClones** - *Constructs clonal consensus sequences*

Description
--------------------

`collapseClones` identifies the consensus sequence of each clonal 
group and appends columns to the input `data.frame` containing the clonal 
consensus and germline for each sequence.


Usage
--------------------
```
collapseClones(db, cloneColumn = "CLONE", sequenceColumn = "SEQUENCE_IMGT",
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

A modified `db` with clonal consensus sequences added 
in the following columns:

+  `CLONAL_SEQUENCE`:  consensus input sequence for the clone.
+  `CLONAL_GERMLINE`:  consensus germline sequence for the clone.
Generally, this will be unchanged from
the data in `germlineColumn`, but
may be truncated when the input sequence
is truncaated due to inconsistencies 
in the lengths of the input sequences or
`regionDefinition` limits.


Details
-------------------

For sequences identified to be part of the same clone, this function defines an 
effective sequence that will be representative for all mutations in the clone. Each 
position in this consensus (or effective) sequence is created by a weighted sampling 
of each mutated base (and non "N", "." or "-" characters) from all the sequences in 
the clone. For example, in a clone with 5 sequences that have a C at position 1, and 
5 sequences with a T at this same position, the consensus sequence will have a C 50%  
and T 50% of the time it is called.

Non-terminal branch mutations are defined as the set of mutations that occur on 
branches of the lineage tree that are not connected to a leaf. For computational 
efficiency, the set of non-terminal branch mutations is approximated as those that are
shared between more than one sequence in a clone. In this case the terminal branch 
mutations are filtered out.



Examples
-------------------

```R
# Subset example data
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")

# Build clonal consensus for the full sequence
clones <- collapseClones(db, nproc=1)

```


```
Collapsing clonal sequences...

```


```R

# Build clonal consensus for V-region only 
# Return the same number of rows as the input
clones <- collapseClones(db, regionDefinition=IMGT_V_NO_CDR3, 
expandedDb=TRUE, nproc=1)
```


```
Collapsing clonal sequences...

```



See also
-------------------

See [IMGT_SCHEMES](IMGT_SCHEMES.md) for a set of predefined [RegionDefinition](RegionDefinition-class.md) objects.



