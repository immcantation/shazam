





**calcObservedMutations** - *Count the number of observed mutations in a sequence.*

Description
--------------------

`calcObservedMutations` determines all the mutations in a given input seqeunce compared
to its germline sequence.


Usage
--------------------
```
calcObservedMutations(inputSeq, germlineSeq, regionDefinition = NULL,
mutationDefinition = NULL, returnRaw = FALSE, frequency = FALSE)
```

Arguments
-------------------

inputSeq
:   input sequence. IUPAC ambiguous characters for DNA are supported.

germlineSeq
:   germline sequence. Germline should **not** contain ambiguous
characters.

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
mutation types, as opposed to counts of mutations across positions.
Also returns the number of bases used as the denominator when 
calculating frequency. Default is `FALSE`.

frequency
:   `logical` indicating whether or not to calculate
mutation frequencies. The denominator used is the number of bases
that are not one of "N", "-", or "." in either the input or the 
germline sequences. If set, this overwrites `returnRaw`. 
Default is `FALSE`.




Value
-------------------

For `returnRaw=FALSE`, an `array` with the numbers of replacement (R) 
and silent (S) mutations. 

For `returnRaw=TRUE`, a list containing 

+  A data frame (`$pos`) whose columns (`position`, `R`, `S`, 
and `region`) indicate the nucleotide position, the number of R mutations, 
the number of S mutations, and the region in which the nucleotide position is 
in.
+  A vector (`$nonN`) indicating the number of bases in regions defined by 
`regionDefinition` (excluding non-triplet overhang, if any) that are not 
one of "N", "-", or "." in either the observed or the germline.


For `frequency=TRUE`, regardless of `returnRaw`, an `array` with the 
frequencies of replacement (R) and silent (S) mutations.


Details
-------------------

Each mutation is considered independently in the germline context. Note, only the part of 
`inputSeq` defined in `regionDefinition` is analyzed. For example, when using 
the default [IMGT_V](IMGT_SCHEMES.md) definition, then mutations in positions beyond 
312 will be ignored. Additionally, non-triplet overhang at the sequence end is ignored.

Only replacement (R) and silent (S) mutations are included in the results. Excluded are: 

+  Stop mutations
+  Mutations occurring in codons where one or both of the observed and the 
germline involve(s) one or more of "N", "-", or ".".

E.g.: the case in which NNN in the germline sequence is observed as NNC in 
the input sequence.

In other words, a result that is `NA` or zero indicates absence of R and S mutations, 
not necessarily all types of mutations, such as the excluded ones mentioned above.



Examples
-------------------

```R
# Use an entry in the example data for input and germline sequence
data(ExampleDb, package="alakazam")
in_seq <- ExampleDb[["SEQUENCE_IMGT"]][100]
germ_seq <-  ExampleDb[["GERMLINE_IMGT_D_MASK"]][100]

# Identify all mutations in the sequence
ex1_raw = calcObservedMutations(in_seq, germ_seq, returnRaw=TRUE)
# Count all mutations in the sequence
ex1_count = calcObservedMutations(in_seq, germ_seq, returnRaw=FALSE)
ex1_freq = calcObservedMutations(in_seq, germ_seq, returnRaw=FALSE, frequency=TRUE)
# Compare this with ex1_count
table(ex1_raw$pos$region, ex1_raw$pos$R)[, "1"]

```


```
[1] 11

```


```R
table(ex1_raw$pos$region, ex1_raw$pos$S)[, "1"]

```


```
[1] 7

```


```R
# Compare this with ex1_freq
table(ex1_raw$pos$region, ex1_raw$pos$R)[, "1"] / ex1_raw$nonN

```


```
       SEQ 
0.03363914 

```


```R
table(ex1_raw$pos$region, ex1_raw$pos$S)[, "1"] / ex1_raw$nonN

```


```
       SEQ 
0.02140673 

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
table(ex2_raw$pos$region, ex2_raw$pos$R)[, "1"]

```


```
CDR FWR 
  4   7 

```


```R
table(ex2_raw$pos$region, ex2_raw$pos$S)[, "1"]                              

```


```
CDR FWR 
  1   4 

```


```R
# Compare this with ex2_freq
table(ex2_raw$pos$region, ex2_raw$pos$R)[, "1"] / ex2_raw$nonN     

```


```
       CDR        FWR 
0.08333333 0.02916667 

```


```R
table(ex2_raw$pos$region, ex2_raw$pos$S)[, "1"] / ex2_raw$nonN                                       

```


```
       CDR        FWR 
0.02083333 0.01666667 

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
table(ex3_raw$pos$region, ex3_raw$pos$R)[, "1"]

```


```
CDR FWR 
  3   4 

```


```R
table(ex3_raw$pos$region, ex3_raw$pos$S)[, "1"]

```


```
CDR FWR 
  2   7 

```


```R
# Compare this with ex3_freq
table(ex3_raw$pos$region, ex3_raw$pos$R)[, "1"] / ex3_raw$nonN                                        

```


```
       CDR        FWR 
0.06250000 0.01666667 

```


```R
table(ex3_raw$pos$region, ex3_raw$pos$S)[, "1"] / ex3_raw$nonN
```


```
       CDR        FWR 
0.04166667 0.02916667 

```



See also
-------------------

See [observedMutations](observedMutations.md) for counting the number of observed mutations 
in a `data.frame`.



