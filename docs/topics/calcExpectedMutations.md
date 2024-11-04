**calcExpectedMutations** - *Calculate expected mutation frequencies of a sequence*

Description
--------------------

`calcExpectedMutations` calculates the expected mutation
frequencies of a given sequence. This is primarily a helper function for
[expectedMutations](expectedMutations.md).


Usage
--------------------
```
calcExpectedMutations(
germlineSeq,
inputSeq = NULL,
targetingModel = HH_S5F,
regionDefinition = NULL,
mutationDefinition = NULL
)
```

Arguments
-------------------

germlineSeq
:   germline (reference) sequence.

inputSeq
:   input (observed) sequence. If this is not `NULL`, 
then `germlineSeq` will be processed to be the same
same length as `inputSeq` and positions in 
`germlineSeq` corresponding to positions with Ns in 
`inputSeq` will also be assigned an N.

targetingModel
:   [TargetingModel](TargetingModel-class.md) object. Default is [HH_S5F](HH_S5F.md).

regionDefinition
:   [RegionDefinition](RegionDefinition-class.md) object defining the regions
and boundaries of the Ig sequences.

mutationDefinition
:   [MutationDefinition](MutationDefinition-class.md) object defining replacement
and silent mutation criteria. If `NULL` then 
replacement and silent are determined by exact 
amino acid identity.




Value
-------------------

A `numeric` vector of the expected frequencies of mutations in the 
regions in the `regionDefinition`. For example, when using the default 
[IMGT_V](IMGT_SCHEMES.md) definition, which defines positions for CDR and 
FWR, the following columns are calculated:

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

`calcExpectedMutations` calculates the expected mutation frequencies of a 
given sequence and its germline. 

Note, only the part of the sequences defined in `regionDefinition` are analyzed. 
For example, when using the default [IMGT_V](IMGT_SCHEMES.md) definition, mutations in
positions beyond 312 will be ignored.



Examples
-------------------

```R
# Load example data
data(ExampleDb, package="alakazam")

# Use first entry in the exampled data for input and germline sequence
in_seq <- ExampleDb[["sequence_alignment"]][1]
germ_seq <-  ExampleDb[["germline_alignment_d_mask"]][1]

# Identify all mutations in the sequence
calcExpectedMutations(germ_seq,in_seq)

```


```
    seq_r     seq_s 
0.7636446 0.2363554 

```


```R

# Identify only mutations the V segment minus CDR3
calcExpectedMutations(germ_seq, in_seq, regionDefinition=IMGT_V)

```


```
     cdr_r      cdr_s      fwr_r      fwr_s 
0.20544721 0.04081758 0.56090228 0.19283293 

```


```R

# Define mutations based on hydropathy
calcExpectedMutations(germ_seq, in_seq, regionDefinition=IMGT_V,
mutationDefinition=HYDROPATHY_MUTATIONS)

```


```
    cdr_r     cdr_s     fwr_r     fwr_s 
0.1209459 0.1253189 0.3169116 0.4368236 

```



See also
-------------------

[expectedMutations](expectedMutations.md) calls this function.
To create a custom `targetingModel` see [createTargetingModel](createTargetingModel.md).
See [calcObservedMutations](calcObservedMutations.md) for getting observed mutation counts.






