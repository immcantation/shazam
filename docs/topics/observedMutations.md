





**observedMutations** - *Calculate observed numbers of mutations*

Description
--------------------

`observedMutations` calculates the observed number of mutations for each 
sequence in the input `data.frame`.


Usage
--------------------
```
observedMutations(db, sequenceColumn = "SEQUENCE_IMGT",
germlineColumn = "GERMLINE_IMGT_D_MASK", frequency = FALSE,
regionDefinition = NULL, mutationDefinition = NULL, nproc = 1)
```

Arguments
-------------------

db
:   `data.frame` containing sequence data.

sequenceColumn
:   `character` name of the column containing input 
sequences.

germlineColumn
:   `character` name of the column containing 
the germline or reference sequence.

frequency
:   `logical` indicating whether or not to calculate
mutation frequencies. Default is `FALSE`.

regionDefinition
:   [RegionDefinition](RegionDefinition-class.md) object defining the regions
and boundaries of the Ig sequences. If NULL, mutations 
are counted for entire sequence.

mutationDefinition
:   [MutationDefinition](MutationDefinition-class.md) object defining replacement
and silent mutation criteria. If `NULL` then 
replacement and silent are determined by exact 
amino acid identity.

nproc
:   number of cores to distribute the operation over. If the 
cluster has already been set the call function with 
`nproc` = 0 to not reset or reinitialize. Default is 
`nproc` = 1.



Value
-------------------

A modified `db` `data.frame` with observed mutation counts for each 
sequence listed. The columns names are dynamically created based on the
regions in the `regionDefinition`. For example, when using the
[IMGT_V_NO_CDR3](IMGT_SCHEMES.md) definition, which defines positions for CDR and
FWR, the following columns are added:

+ `OBSERVED_CDR_R`:  number of replacement mutations in CDR1 and 
CDR2 of the V-segment.
+ `OBSERVED_CDR_S`:  number of silent mutations in CDR1 and CDR2 
of the V-segment.
+ `OBSERVED_FWR_R`:  number of replacement mutations in FWR1, 
FWR2 and FWR3 of the V-segment.
+ `OBSERVED_FWR_S`:  number of silent mutations in FWR1, FWR2 and
FWR3 of the V-segment.


Details
-------------------

Mutation count are determined by comparing the input sequences (in the column specified 
by `sequenceColumn`) to the germline sequence (in the column specified by 
`germlineColumn`). 

The mutations are binned as either replacement (R) or silent (S) across the different 
regions of the sequences as defined by `regionDefinition`. Typically, this would 
be the framework (FWR) and complementarity determining (CDR) regions of IMGT-gapped 
nucleotide sequences. Mutation counts are appended to the input `db` as 
additional columns.



Examples
-------------------

```R
# Subset example data
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")

# Calculate mutation frequency over the entire sequence
db_obs <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
germlineColumn="GERMLINE_IMGT_D_MASK",
frequency=TRUE,
nproc=1)

```


```
Calculating observed number of mutations...

```


```R

# Count of V-region mutations split by FWR and CDR
# With mutations only considered replacement if charge changes
db_obs <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
germlineColumn="GERMLINE_IMGT_D_MASK",
regionDefinition=IMGT_V_NO_CDR3,
mutationDefinition=CHARGE_MUTATIONS,
nproc=1)
```


```
Calculating observed number of mutations...

```



See also
-------------------

[calcObservedMutations](calcObservedMutations.md) is called by this function to get the list of mutations 
in each sequence grouped by the [RegionDefinition](RegionDefinition-class.md). 
See [IMGT_SCHEMES](IMGT_SCHEMES.md) for a set of predefined [RegionDefinition](RegionDefinition-class.md) objects.
See [expectedMutations](expectedMutations.md) for calculating expected mutation frequencies.



