





**calcObservedMutations** - *Count the number of observed mutations in a sequence.*

Description
--------------------

`calcObservedMutations` determines all the mutations in a given input seqeunce compared
to its germline sequence.


Usage
--------------------
```
calcObservedMutations(inputSeq, germlineSeq, frequency = FALSE,
regionDefinition = NULL, mutationDefinition = NULL)
```

Arguments
-------------------

inputSeq
:   input sequence.

germlineSeq
:   germline sequence.

frequency
:   `logical` indicating whether or not to calculate
mutation frequencies. Default is `FALSE`.

regionDefinition
:   [RegionDefinition](RegionDefinition-class.md) object defining the regions
and boundaries of the Ig sequences. Note, only the part of
sequences defined in `regionDefinition` are analyzed.
If NULL, mutations are counted for entire sequence.

mutationDefinition
:   [MutationDefinition](MutationDefinition-class.md) object defining replacement
and silent mutation criteria. If `NULL` then 
replacement and silent are determined by exact 
amino acid identity.



Value
-------------------

An `array` of the mutations, replacement (R) or silent(S), with the 
names indicating the nucleotide postion of the mutations in the sequence.

Details
-------------------

Each mutation is considered independently in its codon context. Note, only the part of 
`inputSeq` defined in `regionDefinition` is analyzed. For example, when using 
the default [IMGT_V_NO_CDR3](IMGT_SCHEMES.md) definition, then mutations in positions beyond 
312 will be ignored.



Examples
-------------------

```R
# Use first entry in the exampled data for input and germline sequence
data(ExampleDb, package="alakazam")
in_seq <- ExampleDb[1, "SEQUENCE_IMGT"]
germ_seq <-  ExampleDb[1, "GERMLINE_IMGT_D_MASK"]

# Identify all mutations in the sequence
calcObservedMutations(in_seq, germ_seq)

```


```
SEQ_R SEQ_S 
    0     2 

```


```R

# Identify only mutations the V segment minus CDR3
calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V_NO_CDR3)

```


```
[1] NA

```


```R
 
# Identify mutations by change in hydropathy class
calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V_NO_CDR3,
mutationDefinition=HYDROPATHY_MUTATIONS, frequency=TRUE)
```


```
[1] NA

```



See also
-------------------

See [calcDBObservedMutations](calcDBObservedMutations.md) for counting the number of observed mutations.



