**calcObservedMutations** - *Count the number of observed mutations in a sequence.*

Description
--------------------

`calcObservedMutations` determines all the mutations in a given input sequence 
compared to its germline sequence.


Usage
--------------------
```
calcObservedMutations(
inputSeq,
germlineSeq,
regionDefinition = NULL,
mutationDefinition = NULL,
ambiguousMode = c("eitherOr", "and"),
returnRaw = FALSE,
frequency = FALSE
)
```

Arguments
-------------------

inputSeq
:   input sequence. IUPAC ambiguous characters for DNA are 
supported.

germlineSeq
:   germline sequence. IUPAC ambiguous characters for DNA 
are supported.

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

ambiguousMode
:   whether to consider ambiguous characters as 
`"either or"` or `"and"` when determining and 
counting the type(s) of mutations. Applicable only if
`inputSeq` and/or `germlineSeq` 
contain(s) ambiguous characters. One of 
`c("eitherOr", "and")`. Default is `"eitherOr"`.

returnRaw
:   return the positions of point mutations and their 
corresponding mutation types, as opposed to counts of 
mutations across positions. Also returns the number of 
bases used as the denominator when calculating frequency. 
Default is `FALSE`.

frequency
:   `logical` indicating whether or not to calculate
mutation frequencies. The denominator used is the number 
of bases that are not one of "N", "-", or "." in either 
the input or the germline sequences. If set, this 
overwrites `returnRaw`. Default is `FALSE`.




Value
-------------------

For `returnRaw=FALSE`, an `array` with the numbers of replacement (R) 
and silent (S) mutations. 

For `returnRaw=TRUE`, a list containing 

+  `$pos`: A data frame whose columns (`position`, `R`, 
`S`, and `region`) indicate, respecitively, the nucleotide 
position, the number of R mutations at that position, the number of S 
mutations at that position, and the region in which that nucleotide
is in.
+  `$nonN`: A vector indicating the number of bases in regions 
defined by `regionDefinition` (excluding non-triplet overhang, 
if any) that are not one of "N", "-", or "." in either the 
`inputSeq` or `germlineSeq`.


For `frequency=TRUE`, regardless of `returnRaw`, an `array` 
with the frequencies of replacement (R) and silent (S) mutations.


Details
-------------------

**Each mutation is considered independently in the germline context**. For illustration,
consider the case where the germline is `TGG` and the observed is `TAC`.
When determining the mutation type at position 2, which sees a change from `G` to 
`A`, we compare the codon `TGG` (germline) to `TAG` (mutation at position
2 independent of other mutations in the germline context). Similarly, when determining 
the mutation type at position 3, which sees a change from `G` to `C`, we 
compare the codon `TGG` (germline) to `TGC` (mutation at position 3 independent 
of other mutations in the germline context).

If specified, only the part of `inputSeq` defined in `regionDefinition` is 
analyzed. For example, when using the default [IMGT_V](IMGT_SCHEMES.md) definition, then mutations 
in positions beyond 312 will be ignored. Additionally, non-triplet overhang at the 
sequence end is ignored.

Only replacement (R) and silent (S) mutations are included in the results. **Excluded**
are: 

+  Stop mutations

E.g.: the case where `TAGTGG` is observed for the germline `TGGTGG`.

+  Mutations occurring in codons where one or both of the observed and the 
germline involve(s) one or more of "N", "-", or ".".

E.g.: the case where `TTG` is observed for the germline being any one of 
`TNG`, `.TG`, or `-TG`. Similarly, the case where any one of 
`TTN`, `TT.`, or `TT-` is observed for the germline `TTG`.


In other words, a result that is `NA` or zero indicates absence of R and S mutations, 
not necessarily all types of mutations, such as the excluded ones mentioned above.

`NA` is also returned if `inputSeq` or `germlineSeq` is shorter than 3
nucleotides.


Ambiguous characters
-------------------


When there are ambiguous characters present, the user could choose how mutations involving
ambiguous characters are counted through `ambiguousMode`. The two available modes 
are `"eitherOr"` and `"and"`. 

+  With `"eitherOr"`, ambiguous characters are each expanded but only 
1 mutation is recorded. When determining the type of mutation, the 
priority for different types of mutations, in decreasing order, is as follows:
no mutation, replacement mutation, silent mutation, and stop mutation. 

When counting the number of non-N, non-dash, and non-dot positions, each
position is counted only once, regardless of the presence of ambiguous
characters.

As an example, consider the case where `germlineSeq` is `"TST"` and
`inputSeq` is `"THT"`. Expanding `"H"` at position 2 in 
`inputSeq` into `"A"`, `"C"`, and `"T"`, as well as 
expanding `"S"` at position 2 in `germlineSeq` into `"C"` and 
`"G"`, one gets:


+  `"TCT"` (germline) to `"TAT"` (observed): replacement
+  `"TCT"` (germline) to `"TCT"` (observed): no mutation
+  `"TCT"` (germline) to `"TTT"` (observed): replacement 
+  `"TGT"` (germline) to `"TAT"` (observed): replacement 
+  `"TGT"` (germline) to `"TCT"` (observed): replacement
+  `"TGT"` (germline) to `"TTT"` (observed): replacement


Because "no mutation" takes priority over replacement mutation, the final 
mutation count returned for this example is `NA` (recall that only R and 
S mutations are returned). The number of non-N, non-dash, and non-dot 
positions is 3.

+  With `"and"`, ambiguous characters are each expanded and mutation(s)
from all expansions are recorded.

When counting the number of non-N, non-dash, and non-dot positions, if a 
position contains ambiguous character(s) in `inputSeq` and/or 
`germlineSeq`, the count at that position is taken to be the total 
number of combinations of germline and observed codons after expansion.

Using the same example from above, the final result returned for this example
is that there are 5 R mutations at position 2. The number of non-N, non-dash,
and non-dot positions is 8, since there are 6 combinations stemming from 
position 2 after expanding the germline codon (`"TST"`) and the observed 
codon (`"THT"`).




Examples
-------------------

```R
# Use an entry in the example data for input and germline sequence
data(ExampleDb, package="alakazam")
in_seq <- ExampleDb[["sequence_alignment"]][100]
germ_seq <-  ExampleDb[["germline_alignment_d_mask"]][100]

# Identify all mutations in the sequence
ex1_raw <- calcObservedMutations(in_seq, germ_seq, returnRaw=TRUE)
# Count all mutations in the sequence
ex1_count <- calcObservedMutations(in_seq, germ_seq, returnRaw=FALSE)
ex1_freq <- calcObservedMutations(in_seq, germ_seq, returnRaw=FALSE, frequency=TRUE)
# Compare this with ex1_count
table(ex1_raw$pos$region, ex1_raw$pos$r)[, "1"]

```


```
[1] 11

```


```R
table(ex1_raw$pos$region, ex1_raw$pos$s)[, "1"]

```


```
[1] 7

```


```R
# Compare this with ex1_freq
table(ex1_raw$pos$region, ex1_raw$pos$r)[, "1"]/ex1_raw$nonN

```


```
       seq 
0.03363914 

```


```R
table(ex1_raw$pos$region, ex1_raw$pos$s)[, "1"]/ex1_raw$nonN

```


```
       seq 
0.02140673 

```


```R

# Identify only mutations the V segment minus CDR3
ex2_raw <- calcObservedMutations(in_seq, germ_seq, 
regionDefinition=IMGT_V, returnRaw=TRUE)
# Count only mutations the V segment minus CDR3
ex2_count <- calcObservedMutations(in_seq, germ_seq, 
regionDefinition=IMGT_V, returnRaw=FALSE)
ex2_freq <- calcObservedMutations(in_seq, germ_seq, 
regionDefinition=IMGT_V, returnRaw=FALSE,
frequency=TRUE)
# Compare this with ex2_count
table(ex2_raw$pos$region, ex2_raw$pos$r)[, "1"]

```


```
cdr fwr 
  4   7 

```


```R
table(ex2_raw$pos$region, ex2_raw$pos$s)[, "1"]                              

```


```
cdr fwr 
  1   4 

```


```R
# Compare this with ex2_freq
table(ex2_raw$pos$region, ex2_raw$pos$r)[, "1"]/ex2_raw$nonN     

```


```
       cdr        fwr 
0.08333333 0.02916667 

```


```R
table(ex2_raw$pos$region, ex2_raw$pos$s)[, "1"]/ex2_raw$nonN                                       

```


```
       cdr        fwr 
0.02083333 0.01666667 

```


```R

# Identify mutations by change in hydropathy class
ex3_raw <- calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V,
mutationDefinition=HYDROPATHY_MUTATIONS, 
returnRaw=TRUE)
# Count mutations by change in hydropathy class
ex3_count <- calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V,
mutationDefinition=HYDROPATHY_MUTATIONS, 
returnRaw=FALSE)
ex3_freq <- calcObservedMutations(in_seq, germ_seq, regionDefinition=IMGT_V,
mutationDefinition=HYDROPATHY_MUTATIONS, 
returnRaw=FALSE, frequency=TRUE)
# Compre this with ex3_count
table(ex3_raw$pos$region, ex3_raw$pos$r)[, "1"]

```


```
cdr fwr 
  3   4 

```


```R
table(ex3_raw$pos$region, ex3_raw$pos$s)[, "1"]

```


```
cdr fwr 
  2   7 

```


```R
# Compare this with ex3_freq
table(ex3_raw$pos$region, ex3_raw$pos$r)[, "1"]/ex3_raw$nonN                                        

```


```
       cdr        fwr 
0.06250000 0.01666667 

```


```R
table(ex3_raw$pos$region, ex3_raw$pos$s)[, "1"]/ex3_raw$nonN
```


```
       cdr        fwr 
0.04166667 0.02916667 

```



See also
-------------------

See [observedMutations](observedMutations.md) for counting the number of observed mutations 
in a `data.frame`.






