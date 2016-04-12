





**calcExpectedMutations** - *Calculate expected mutation frequencies of a sequence*

Description
--------------------

`calcExpectedMutations` calculates the expected mutation
frequencies of a given sequence. This is primarily a helper function for
[calcDBExpectedMutations](calcDBExpectedMutations.md).


Usage
--------------------
```
calcExpectedMutations(germlineSeq, inputSeq = NULL,
targetingModel = HS5FModel, regionDefinition = NULL,
mutationDefinition = NULL)
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
:   [TargetingModel](TargetingModel-class.md) object. Default is [HS5FModel](HS5FModel.md).

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
[IMGT_V_NO_CDR3](IMGT_SCHEMES.md) definition, which defines positions for CDR and 
FWR, the following columns are calculated:

+ `EXPECTED_CDR_R`:  number of replacement mutations in CDR1 and 
CDR2 of the V-segment.
+ `EXPECTED_CDR_S`:  number of silent mutations in CDR1 and CDR2 
of the V-segment.
+ `EXPECTED_FWR_R`:  number of replacement mutations in FWR1, 
FWR2 and FWR3 of the V-segment.
+ `EXPECTED_FWR_S`:  number of silent mutations in FWR1, FWR2 and
FWR3 of the V-segment.


Details
-------------------

`calcExpectedMutations` calculates the expected mutation frequencies of a 
given sequence and its germline. 

Note, only the part of the sequences defined in `regionDefinition` are analyzed. 
For example, when using the default [IMGT_V_NO_CDR3](IMGT_SCHEMES.md) definition, mutations in
positions beyond 312 will be ignored.



Examples
-------------------

```R
# Extracting the first entry in the exampled data to use for input and germline sequences.
inputSeq <- InfluenzaDb[1, "SEQUENCE_IMGT"]
germlineSeq <-  InfluenzaDb[1, "GERMLINE_IMGT_D_MASK"]

# Identify all mutations in the sequence
calcExpectedMutations(inputSeq, germlineSeq)

```


```
    SEQ_R     SEQ_S 
0.7379522 0.2620478 

```


```R

# Identify only mutations the V segment minus CDR3
calcExpectedMutations(inputSeq, germlineSeq, regionDefinition=IMGT_V_NO_CDR3)

```


```
     CDR_R      CDR_S      FWR_R      FWR_S 
0.16184071 0.04872069 0.57766105 0.21177756 

```


```R

# Define mutations based on hydropathy
calcExpectedMutations(inputSeq, germlineSeq, regionDefinition=IMGT_V_NO_CDR3,
mutationDefinition=HYDROPATHY_MUTATIONS)
```


```
    CDR_R     CDR_S     FWR_R     FWR_S 
0.1043777 0.1061837 0.3317107 0.4577279 

```



See also
-------------------

[calcDBExpectedMutations](calcDBExpectedMutations.md) calls this function.
To create a custom `targetingModel` see [createTargetingModel](createTargetingModel.md).
See [calcObservedMutations](calcObservedMutations.md) for getting observed mutation counts.



