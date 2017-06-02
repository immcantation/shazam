





**collapseClones** - *Constructs effective clonal sequences for all clones*

Description
--------------------

`collapseClones` creates effective input and germline sequences for each clonal 
group and appends columns containing the consensus sequences to the input 
`data.frame`.


Usage
--------------------
```
collapseClones(db, cloneColumn = "CLONE", sequenceColumn = "SEQUENCE_IMGT",
germlineColumn = "GERMLINE_IMGT_D_MASK", expandedDb = FALSE,
regionDefinition = NULL, method = c("thresholdedFreq", "mostCommon",
"catchAll"), minimumFrequency = 0.6, nproc = 1)
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
to returning just one sequence per clone).

regionDefinition
:   [RegionDefinition](RegionDefinition-class.md) object defining the regions
and boundaries of the Ig sequences.

method
:   method for calculating input consensus sequence. One of 
`"thresholdedFreq"`, `"mostCommon"`, or 
`"catchAll"`. See [calcClonalConsensusHelper](calcClonalConsensusHelper.md) 
for details.

minimumFrequency
:   frequency threshold for calculating input consensus sequence.
Required for the `"thresholdedFreq"` method. Default is
0.6. See [calcClonalConsensusHelper](calcClonalConsensusHelper.md) for details.

nproc
:   Number of cores to distribute the operation over. If the 
`cluster` has already been set earlier, then pass the 
`cluster`. This will ensure that it is not reset.




Value
-------------------

A modified `db` with clonal consensus sequences added 
in the following columns:

+  `CLONAL_SEQUENCE`:  effective sequence for the clone.
+  `CLONAL_GERMLINE`:  germline sequence for the clone.



Details
-------------------

See `Details` for [calcClonalConsensus](calcClonalConsensus.md) and [calcClonalConsensusHelper](calcClonalConsensusHelper.md)
for explanation on how the different methods work.



Examples
-------------------

```R
# Subset example data
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")

# Build clonal consensus for the full sequence
clones <- collapseClones(db, method="thresholdedFreq", 
minimumFrequency=0.65, nproc=1)

```


```
Collapsing clonal sequences...

```


```R

# Build clonal consensus for V-region only 
# Return the same number of rows as the input
clones <- collapseClones(db, method="mostCommon", 
regionDefinition=IMGT_V, 
expandedDb=TRUE, nproc=1)
```


```
Collapsing clonal sequences...

```



See also
-------------------

See [IMGT_SCHEMES](IMGT_SCHEMES.md) for a set of predefined [RegionDefinition](RegionDefinition-class.md) objects.
See [calcClonalConsensus](calcClonalConsensus.md) and [calcClonalConsensusHelper](calcClonalConsensusHelper.md) for helper functions.



