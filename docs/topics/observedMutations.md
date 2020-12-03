**observedMutations** - *Calculate observed numbers of mutations*

Description
--------------------

`observedMutations` calculates the observed number of mutations for each 
sequence in the input `data.frame`.


Usage
--------------------
```
observedMutations(
db,
sequenceColumn = "sequence_alignment",
germlineColumn = "germline_alignment_d_mask",
regionDefinition = NULL,
mutationDefinition = NULL,
ambiguousMode = c("eitherOr", "and"),
frequency = FALSE,
combine = FALSE,
nproc = 1,
cloneColumn = "clone_id",
juncLengthColumn = "junction_length"
)
```

Arguments
-------------------

db
:   `data.frame` containing sequence data.

sequenceColumn
:   `character` name of the column containing input 
sequences. IUPAC ambiguous characters for DNA are 
supported.

germlineColumn
:   `character` name of the column containing 
the germline or reference sequence. IUPAC ambiguous 
characters for DNA are supported.

regionDefinition
:   [RegionDefinition](RegionDefinition-class.md) object defining the regions
and boundaries of the Ig sequences. If NULL, mutations 
are counted for entire sequence.

mutationDefinition
:   [MutationDefinition](MutationDefinition-class.md) object defining replacement
and silent mutation criteria. If `NULL` then 
replacement and silent are determined by exact 
amino acid identity.

ambiguousMode
:   whether to consider ambiguous characters as 
`"either or"` or `"and"` when determining and 
counting the type(s) of mutations. Applicable only if
`sequenceColumn` and/or `germlineColumn` 
contain(s) ambiguous characters. One of 
`c("eitherOr", "and")`. Default is `"eitherOr"`.

frequency
:   `logical` indicating whether or not to calculate
mutation frequencies. Default is `FALSE`.

combine
:   `logical` indicating whether for each sequence should
the mutation counts for the different regions (CDR, FWR) and 
mutation types be combined and return one value of 
count/frequency per sequence instead of 
multiple values. Default is `FALSE`.

nproc
:   number of cores to distribute the operation over. If the 
cluster has already been set the call function with 
`nproc` = 0 to not reset or reinitialize. Default is 
`nproc` = 1.

cloneColumn
:   clone id column name in `db`

juncLengthColumn
:   junction length column name in `db`




Value
-------------------

A modified `db` `data.frame` with observed mutation counts for each 
sequence listed. The columns names are dynamically created based on the
regions in the `regionDefinition`. For example, when using the
[IMGT_V](IMGT_SCHEMES.md) definition, which defines positions for CDR and
FWR, the following columns are added:

+ `mu_count_cdr_r`:  number of replacement mutations in CDR1 and 
CDR2 of the V-segment.
+ `mu_count_cdr_s`:  number of silent mutations in CDR1 and CDR2 
of the V-segment.
+ `mu_count_fwr_r`:  number of replacement mutations in FWR1, 
FWR2 and FWR3 of the V-segment.
+ `mu_count_fwr_s`:  number of silent mutations in FWR1, FWR2 and
FWR3 of the V-segment.

If `frequency=TRUE`, R and S mutation frequencies are
calculated over the number of non-N positions in the specified regions.

+ `mu_freq_cdr_r`:  frequency of replacement mutations in CDR1 and 
CDR2 of the V-segment.
+ `mu_freq_cdr_s`:  frequency of silent mutations in CDR1 and CDR2 
of the V-segment.
+ `mu_freq_fwr_r`:  frequency of replacement mutations in FWR1, 
FWR2 and FWR3 of the V-segment.
+ `mu_freq_fwr_s`:  frequency of silent mutations in FWR1, FWR2 and
FWR3 of the V-segment.
 
If `frequency=TRUE` and `combine=TRUE`, the mutations and non-N positions
are aggregated and a single `mu_freq` value is returned

+ `mu_freq`:  frequency of replacement and silent mutations in the 
specified region



Details
-------------------

Mutation counts are determined by comparing the input sequences (in the column specified 
by `sequenceColumn`) to a reference sequence. If `db` includes lineage information,
e.g. in the field `parent_sequence`, the reference sequence can be set to  
use that field as reference sequence. Ssee more details in [makeGraphDf](makeGraphDf.md)).
See [calcObservedMutations](calcObservedMutations.md) for more technical details, 
**including criteria for which sequence differences are included in the mutation 
counts and which are not**.

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
db <- subset(ExampleDb, c_call == "IGHG" & sample_id == "+7d")

# Calculate mutation frequency over the entire sequence
db_obs <- observedMutations(db, sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask",
frequency=TRUE,
nproc=1)

# Count of V-region mutations split by FWR and CDR
# With mutations only considered replacement if charge changes
db_obs <- observedMutations(db, sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask",
regionDefinition=IMGT_V,
mutationDefinition=CHARGE_MUTATIONS,
nproc=1)

# Count of VDJ-region mutations, split by FWR and CDR
db_obs <- observedMutations(db, sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask",
regionDefinition=IMGT_VDJ,
nproc=1)    
### Not run:
# Count of VDJ-region mutations, split by FWR and CDR
# # This doesn't work because 'parent_sequence' doesn't exist,
# # it should be calculated before
# Update example to include how to create that column.
# db_obs <- observedMutations(db, sequenceColumn="parent_sequence",
# germlineColumn="germline_alignment_d_mask",
# regionDefinition=IMGT_VDJ,
# nproc=1)
```



See also
-------------------

[calcObservedMutations](calcObservedMutations.md) is called by this function to get the number of mutations 
in each sequence grouped by the [RegionDefinition](RegionDefinition-class.md). 
See [IMGT_SCHEMES](IMGT_SCHEMES.md) for a set of predefined [RegionDefinition](RegionDefinition-class.md) objects.
See [expectedMutations](expectedMutations.md) for calculating expected mutation frequencies.






