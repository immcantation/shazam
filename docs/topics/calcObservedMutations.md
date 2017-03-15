





**calcObservedMutations** - *Count the number of observed mutations in a sequence.*

Description
--------------------

`calcObservedMutations` determines all the mutations in a given input seqeunce compared
to its germline sequence.


Usage
--------------------
```
calcObservedMutations(inputSeq, germlineSeq, frequency = FALSE,
regionDefinition = NULL, mutationDefinition = NULL, returnRaw = FALSE)
```

Arguments
-------------------

inputSeq
:   input sequence.

germlineSeq
:   germline sequence.

frequency
:   `logical` indicating whether or not to calculate
mutation frequencies. The denominator used is the number of bases
that are non-N in both the input and the germline sequences.
Default is `FALSE`.

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

returnRaw
:   return the positions of point mutations and their corresponding
mutation types, as opposed to counts of mutations.
Also returns the number of non-N bases used as the denominator when
calculating frequency. Default is `FALSE`.




Value
-------------------

For `returnRaw=FALSE`, an `array` with the number of replacement (R) 
and silent(S) mutations. For `returnRaw=TRUE`, a list containing a data 
frame (`$pos`) whose columns (`position`, `type`, and `region`) 
indicate the position, mutation type (R or S), and region of each mutation; and a 
vector (`$nonN`) indicating the number of non-N bases in regions defined by
`regionDefinition`.


Details
-------------------

Each mutation is considered independently in the germline context. Note, only the part of 
`inputSeq` defined in `regionDefinition` is analyzed. For example, when using 
the default [IMGT_V](IMGT_SCHEMES.md) definition, then mutations in positions beyond 
312 will be ignored. 

Note that only replacement (R) and silent (S) mutations are included in the 
results. Stop mutations and mutations such as the case in which NNN in the germline
sequence is observed as NNC in the input sequence are excluded. In other words,
a result that is `NA` or zero indicates absence of R and S mutations, not 
necessarily all types of mutations, such as the excluded ones mentioned above.



Examples
-------------------

```R
# Use an entry in the example data for input and germline sequence
data(ExampleDb, package="alakazam")
in_seq <- ExampleDb[100, "SEQUENCE_IMGT"]
germ_seq <-  ExampleDb[100, "GERMLINE_IMGT_D_MASK"]

# Identify all mutations in the sequence
ex1_raw = calcObservedMutations(in_seq, germ_seq, returnRaw=TRUE)
# Count all mutations in the sequence
ex1_count = calcObservedMutations(in_seq, germ_seq, returnRaw=FALSE)
ex1_freq = calcObservedMutations(in_seq, germ_seq, returnRaw=FALSE, frequency=TRUE)
# Compare this with ex1_count
table(ex1_raw$pos$region, ex1_raw$pos$type)

```


```
     
       R  S
  SEQ 11  7

```


```R
# Compare this with ex1_freq
table(ex1_raw$pos$region, ex1_raw$pos$type) / ex1_raw$nonN

```


```
     
               R          S
  SEQ 0.03353659 0.02134146

```


```R

# Identify only mutations the V segment minus CDR3
ex2_raw = calcObservedMutations(in_seq, germ_seq, 
regionDefinition=IMGT_V, returnRaw=TRUE)
# Count only mutations the V segment minus CDR3
ex2_count = calcObservedMutations(in_seq, germ_seq, 
regionDefinition=IMGT_V, returnRaw=FALSE)
ex2_freq = calcObservedMutations(in_seq, germ_seq, 
regionDefinition=IMGT_V, returnRaw=FALSE,
frequency=TRUE)
# Compare this with ex2_count
table(ex2_raw$pos$region, ex2_raw$pos$type)                                 

```


```
     
      R S
  CDR 4 1
  FWR 7 4

```


```R
# Compare this with ex2_freq
table(ex2_raw$pos$region, ex2_raw$pos$type) / ex2_raw$nonN                                        

```


```
     
               R          S
  CDR 0.08333333 0.02083333
  FWR 0.02916667 0.01666667

```


```R

# Identify mutations by change in hydropathy class
ex3_raw = calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V,
mutationDefinition=HYDROPATHY_MUTATIONS, returnRaw=TRUE)
# Count mutations by change in hydropathy class
ex3_count = calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V,
mutationDefinition=HYDROPATHY_MUTATIONS, returnRaw=FALSE)
ex3_freq = calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V,
mutationDefinition=HYDROPATHY_MUTATIONS, returnRaw=FALSE, 
frequency=TRUE)
# Compre this with ex3_count
table(ex3_raw$pos$region, ex3_raw$pos$type)                                        

```


```
     
      R S
  CDR 3 2
  FWR 4 7

```


```R
# Compare this with ex3_freq
table(ex3_raw$pos$region, ex3_raw$pos$type) / ex3_raw$nonN
```


```
     
               R          S
  CDR 0.06250000 0.04166667
  FWR 0.01666667 0.02916667

```



See also
-------------------

See [observedMutations](observedMutations.md) for counting the number of observed mutations 
in a `data.frame`.



