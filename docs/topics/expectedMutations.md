**expectedMutations** - *Calculate expected mutation frequencies*

Description
--------------------

`expectedMutations` calculates the expected mutation frequencies for each 
sequence in the input `data.frame`.


Usage
--------------------
```
expectedMutations(
db,
sequenceColumn = "sequence_alignment",
germlineColumn = "germline_alignment",
targetingModel = HH_S5F,
regionDefinition = NULL,
mutationDefinition = NULL,
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
sequences.

germlineColumn
:   `character` name of the column containing 
the germline or reference sequence.

targetingModel
:   [TargetingModel](TargetingModel-class.md) object. Default is [HH_S5F](HH_S5F.md).

regionDefinition
:   [RegionDefinition](RegionDefinition-class.md) object defining the regions
and boundaries of the Ig sequences. To use regions definitions,
sequences in `sequenceColum` and `germlineColumn`
must be aligned, following the IMGT schema.

mutationDefinition
:   [MutationDefinition](MutationDefinition-class.md) object defining replacement
and silent mutation criteria. If `NULL` then 
replacement and silent are determined by exact 
amino acid identity.

nproc
:   `numeric` number of cores to distribute the operation
over. If the cluster has already been set the call function with 
`nproc` = 0 to not reset or reinitialize. Default is 
`nproc` = 1.

cloneColumn
:   clone id column name in `db`

juncLengthColumn
:   junction length column name in `db`




Value
-------------------

A modified `db` `data.frame` with expected mutation frequencies 
for each region defined in `regionDefinition`.

The columns names are dynamically created based on the regions in  
`regionDefinition`. For example, when using the [IMGT_V](IMGT_SCHEMES.md)
definition, which defines positions for CDR and FWR, the following columns are
added:  

+ `mu_expected_cdr_r`:  number of replacement mutations in CDR1 and 
CDR2 of the V-segment.
+ `mu_expected_cdr_s`:  number of silent mutations in CDR1 and CDR2 
of the V-segment.
+ `mu_expected_fwr_r`:  number of replacement mutations in FWR1, 
FWR2 and FWR3 of the V-segment.
+ `mu_expected_fwr_s`:  number of silent mutations in FWR1, FWR2 and
FWR3 of the V-segment.



Details
-------------------

Only the part of the sequences defined in `regionDefinition` are analyzed. 
For example, when using the [IMGT_V](IMGT_SCHEMES.md) definition, mutations in
positions beyond 312 will be ignored.



Examples
-------------------

```R
# Subset example data
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, c_call %in% c("IGHA", "IGHG") & sample_id == "+7d")

# Calculate expected mutations over V region
db_exp <- expectedMutations(db,
sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask",
regionDefinition=IMGT_V,
nproc=1)

# Calculate hydropathy expected mutations over V region
db_exp <- expectedMutations(db,
sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask",
regionDefinition=IMGT_V,
mutationDefinition=HYDROPATHY_MUTATIONS,
nproc=1)
```



See also
-------------------

[calcExpectedMutations](calcExpectedMutations.md) is called by this function to calculate the expected 
mutation frequencies. See [observedMutations](observedMutations.md) for getting observed 
mutation counts. See [IMGT_SCHEMES](IMGT_SCHEMES.md) for a set of predefined 
[RegionDefinition](RegionDefinition-class.md) objects.






