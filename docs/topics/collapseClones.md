





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
germlineColumn = "GERMLINE_IMGT_D_MASK", muFreqColumn = NULL,
regionDefinition = NULL, method = c("mostCommon", "thresholdedFreq",
"catchAll", "mostMutated", "leastMutated"), minimumFrequency = NULL,
includeAmbiguous = FALSE, breakTiesStochastic = FALSE,
breakTiesByColumns = NULL, expandedDb = FALSE, nproc = 1)
```

Arguments
-------------------

db
:   `data.frame` containing sequence data. Required.

cloneColumn
:   `character` name of the column containing clonal 
identifiers. Required.

sequenceColumn
:   `character` name of the column containing input 
sequences. Required. The length of each input sequence should 
match that of its corresponding germline sequence.

germlineColumn
:   `character` name of the column containing germline 
sequences. Required. The length of each germline sequence should 
match that of its corresponding input sequence.

muFreqColumn
:   `character` name of the column containing mutation
frequency. Optional. Applicable to the `"mostMutated"`
and `"leastMutated"` methods. If not supplied, mutation
frequency is computed by calling `observedMutations`.
Default is `NULL`. See Cautions for note on usage.

regionDefinition
:   [RegionDefinition](RegionDefinition-class.md) object defining the regions
and boundaries of the Ig sequences. Optional. Default is 
`NULL`.

method
:   method for calculating input consensus sequence. Required. One of 
`"thresholdedFreq"`, `"mostCommon"`, `"catchAll"`,
`"mostMutated"`, or `"leastMutated"`. See "Methods" for 
details.

minimumFrequency
:   frequency threshold for calculating input consensus sequence.
Applicable to and required for the `"thresholdedFreq"` method.
A canonical choice is 0.6. Default is `NULL`.

includeAmbiguous
:   whether to use ambiguous characters to represent positions at
which there are multiple characters with frequencies that are at least
`minimumFrequency` or that are maximal (i.e. ties). Applicable to 
and required for the `"thresholdedFreq"` and `"mostCommon"` 
methods. Default is `FALSE`. See "Choosing ambiguous characters" 
for rules on choosing ambiguous characters.

breakTiesStochastic
:   In case of ties, whether to randomly pick a sequence from sequences that
fulfill the criteria as consensus. Applicable to and required for all methods
except for `"catchAll"`. Default is `FALSE`. See "Methods" for 
details.

breakTiesByColumns
:   A list of the form `list(c(col_1, col_2, ...), c(fun_1, fun_2, ...))`, 
where `col_i` is a `character` name of a column in `db`,
and `fun_i` is a function to be applied on that column. Currently, 
only `max` and `min` are supported. Note that the two `c()`'s
in `list()` are essential (i.e. if there is only 1 column, the list 
should be of the form `list(c(col_1), c(func_1))`. Applicable to and 
optional for the `"mostMutated"` and `"leastMutated"` methods. 
If supplied, `fun_i`'s are applied on `col_i`'s to help break 
ties. Default is `NULL`. See "Methods" for details.

expandedDb
:   `logical` indicating whether or not to return the 
expanded `db`, containing all the sequences (as opposed
to returning just one sequence per clone).

nproc
:   Number of cores to distribute the operation over. If the 
`cluster` has already been set earlier, then pass the 
`cluster`. This will ensure that it is not reset.




Value
-------------------

A modified `db` with the following additional columns: 

+  `CLONAL_SEQUENCE`:  effective sequence for the clone.
+  `CLONAL_GERMLINE`:  germline sequence for the clone.
+  `CLONAL_SEQUENCE_MUFREQ`:  mutation frequency of 
`CLONAL_SEQUENCE`; only added for the `"mostMutated"`
and `"leastMutated"` methods.


`CLONAL_SEQUENCE` is generated with the method of choice indicated 
by `method`, and `CLONAL_GERMLINE` is generated with the 
`"mostCommon"` method, along with, where applicable, user-defined 
parameters such as `minimumFrequency`, `includeAmbiguous`, 
`breakTiesStochastic`, and `breakTiesByColumns`.


Consensus lengths
-------------------

 For each clone, `CLONAL_SEQUENCE` and 
`CLONAL_GERMLINE` have the same length. 


+  For the `"thresholdedFreq"`, `"mostCommon"`, and `"catchAll"`
methods:

The length of the consensus sequences is determined by the longest possible
consensus sequence (baesd on `inputSeq` and `germlineSeq`) and 
`regionDefinition@seqLength` (if supplied), whichever is shorter.

Given a set of sequences of potentially varying lengths, the longest possible 
length of their consensus sequence is taken to be the longest length along 
which there is information contained at every nucleotide position across 
majority of the sequences. Majority is defined to be greater than 
`floor(n/2)`, where `n` is the number of sequences. If the longest 
possible consensus length is 0, there will be a warning and an empty string 
(`""`) will be returned. 

If a length limit is defined by supplying a `regionDefinition` via 
`regionDefinition@seqLength`, the consensus length will be further 
restricted to the shorter of the longest possible length and 
`regionDefinition@seqLength`.

+  For the `"mostMutated"` and `"leastMutated"` methods:

The length of the consensus sequences depends on that of the most/least 
mutated input sequence, and, if supplied, the length limit defined by 
`regionDefinition@seqLength`, whichever is shorter. If the germline 
consensus computed using the `"mostCommon"` method is longer than 
the most/least mutated input sequence, the germline consensus is trimmed 
to be of the same length as the input consensus.




Methods
-------------------

 The descriptions below use "sequences" as a generalization of input sequences
and germline sequences. 


+  `method="thresholdedFreq"`

A threshold must be supplied to the argument `minimumFrequency`. At 
each position along the length of the consensus sequence, the frequency 
of each nucleotide/character across sequences is tabulated. The 
nucleotide/character whose frequency is at least (i.e. `>=`) 
`minimumFrequency` becomes the consensus; if there is none, the
consensus nucleotide will be `"N"`.

When there are ties (frequencies of multiple nucleotides/characters 
are at least `minimumFrequency`), this method can be deterministic 
or stochastic, depending on additional parameters.


+  With `includeAmbiguous=TRUE`, ties are resolved 
deterministically by representing ties using ambiguous characters. 
See "Choosing ambiguous characters" for how ambiguous characters 
are chosen.
+  With `breakTiesStochastic=TRUE`, ties are resolved 
stochastically by randomly picking a character amongst the ties.
+  When both `TRUE`, `includeAmbiguous` takes precedence 
over `breakTiesStochastic`.
+  When both `FALSE`, the first character from the ties is 
taken to be the consensus following the order of `"A"`, 
`"T"`, `"G"`, `"C"`, `"N"`, `"."`, and 
`"-"`.


Below are some examples looking at a single position based on 5 sequences 
with `minimumFrequency=0.6`, `includeAmbiguous=FALSE`, and 
`breakTiesStochastic=FALSE`:


+  If the sequences have `"A"`, `"A"`, `"A"`, 
`"T"`, `"C"`, the consensus will be `"A"`, because 
`"A"` has frequency 0.6, which is at least 
`minimumFrequency`.
+  If the sequences have `"A"`, `"A"`, `"T"`, 
`"T"`, `"C"`, the consensus will be `"N"`, because 
none of `"A"`, `"T"`, or `"C"` has frequency that 
is at least `minimumFrequency`.


+  `method="mostCommon"`

The most frequent nucleotide/character across sequences at each position 
along the length of the consensus sequence makes up the consensus.

When there are ties (multiple nucleotides/characters with equally maximal 
frequencies), this method can be deterministic or stochastic, depending on 
additional parameters. The same rules for breaking ties for 
`method="thresholdedFreq"` apply.

Below are some examples looking at a single position based on 5 sequences
with `includeAmbiguous=FALSE`, and `breakTiesStochastic=FALSE`:


+  If the sequences have `"A"`, `"A"`, `"T"`, 
`"A"`, `"C"`, the consensus will be `"A"`.
+  If the sequences have `"T"`, `"T"`, `"C"`, 
`"C"`, `"G"`, the consensus will be `"T"`, 
because `"T"` is before `"C"` in the order of 
`"A"`, `"T"`, `"G"`, `"C"`, `"N"`, 
`"."`, and `"-"`. 



+  `method="catchAll"`

This method returns a consensus sequence capturing most of the information 
contained in the sequences. Ambiguous characters are used where applicable.
See "Choosing ambiguous characters" for how ambiguous characters are chosen.
This method is deterministic and does not involve breaking ties.

Below are some examples for `method="catchAll"` looking at a single 
position based on 5 sequences:


+  If the sequences have `"N"`, `"N"`, `"N"`, 
`"N"`, `"N"`, the consensus will be `"N"`.
+  If the sequences have `"N"`, `"A"`, `"A"`, 
`"A"`, `"A"`, the consensus will be `"A"`.
+  If the sequences have `"N"`, `"A"`, `"G"`, 
`"A"`, `"A"`, the consensus will be `"R"`.
+  If the sequences have `"-"`, `"-"`, `"."`, 
`"."`, `"."`, the consensus will be `"-"`.
+  If the sequences have `"-"`, `"-"`, `"-"`, 
`"-"`, `"-"`, the consensus will be `"-"`.
+  If the sequences have `"."`, `"."`, `"."`, 
`"."`, `"."`, the consensus will be `"."`.


+  `method="mostMutated"` and `method="leastMutated"`

These methods return the most/least mutated sequence as the consensus 
sequence. 

When there are ties (multple sequences have the maximal/minimal mutation
frequency), this method can be deterministic or stochastic, depending on 
additional parameters.


+  With `breakTiesStochastic=TRUE`, ties are resolved 
stochastically by randomly picking a sequence out of sequences 
with the maximal/minimal mutation frequency.
+  When `breakTiesByColumns` is supplied, ties are resolved
deterministically. Column by column, a function is applied on 
the column and sequences with column value matching the 
functional value are retained, until ties are resolved or columns
run out. In the latter case, the first remaining sequence is 
taken as the consensus.
+  When `breakTiesStochastic=TRUE` and `breakTiesByColumns` 
is also supplied, `breakTiesStochastic` takes precedence 
over `breakTiesByColumns`.
+  When `breakTiesStochastic=FALSE` and `breakTiesByColumns` 
is not supplied (i.e. `NULL`), the sequence that appears first
amongst the ties is taken as the consensus.





Choosing ambiguous characters
-------------------

 

Ambiguous characters may be present in the returned consensuses when using the
`"catchAll"` method and when using the `"thresholdedFreq"` or 
`"mostCommon"` methods with `includeAmbiguous=TRUE`. 

The rules on choosing ambiguous characters are as follows:


+  If a position contains only `"N"` across sequences, the consensus 
at that position is `"N"`.
+  If a position contains one or more of `"A"`, `"T"`, `"G"`, 
or `"C"`, the consensus will be an IUPAC character representing all 
of the characters present, regardless of whether `"N"`, `"-"`, or
`"."` is present.
+  If a position contains only `"-"` and `"."` across sequences, 
the consensus at thatp osition is taken to be `"-"`. 
+  If a position contains only one of `"-"` or `"."` across 
sequences, the consensus at that position is taken to be the character 
present. 



Cautions
-------------------

 


+ Note that this function does not perform multiple sequence alignment. 
As a prerequisite, it is assumed that the sequences in 
`sequenceColumn` and `germlineColumn` have been aligned somehow. 
In the case of immunoglobulin repertoire analysis, this usually means that 
the sequences are IMGT-gapped.
+ When using the `"mostMutated"` and `"leastMutated"` methods, 
if you supply both `muFreqColumn` and `regionDefinition`,
it is your responsibility to ensure that the mutation frequency in
`muFreqColumn` was calculated with sequence lengths restricted 
to the **same** `regionDefinition` you are supplying. Otherwise, 
the "most/least mutated" sequence you obtain might not be the most/least 
mutated given the `regionDefinition` supplied, because your mutation
frequency was based on a `regionDefinition` different from the one 
supplied.
+ If you intend to run `collapseClones` before 
building a 5-mer targeting model, you **must** choose 
parameters such that your collapsed clonal consensuses do 
**not** include ambiguous characters. This is because the 
targeting model functions do NOT support ambiguous characters 
in their inputs.




Examples
-------------------

```R
# Subset example data
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d" &
CLONE %in% c("3100", "3141", "3184"))
# make a copy of db that has a mutation frequency column
db2 <- observedMutations(db, frequency=TRUE, combine=TRUE)

# Build clonal consensus for the full sequence

# thresholdedFreq method, resolving ties deterministically without using ambiguous characters
clones1 <- collapseClones(db, method="thresholdedFreq", minimumFrequency=0.6,
includeAmbiguous=FALSE, breakTiesStochastic=FALSE)

```

*When both includeAmbiguous and breakTiesStochastic are FALSE, ties are broken in the order of 'A', 'T', 'G', 'C', 'N', '.', and '-'.*
```
Collapsing clonal sequences...

```


```R
# thresholdedFreq method, resolving ties deterministically using ambiguous characters
clones2 <- collapseClones(db, method="thresholdedFreq", minimumFrequency=0.6,
includeAmbiguous=TRUE, breakTiesStochastic=FALSE)

```


```
Collapsing clonal sequences...

```


```R
# thresholdedFreq method, resolving ties stochastically
clones3 <- collapseClones(db, method="thresholdedFreq", minimumFrequency=0.6,
includeAmbiguous=FALSE, breakTiesStochastic=TRUE)

```


```
Collapsing clonal sequences...

```


```R

# mostCommon method, resolving ties deterministically without using ambiguous characters
clones4 <- collapseClones(db, method="mostCommon", 
includeAmbiguous=FALSE, breakTiesStochastic=FALSE)

```

*When both includeAmbiguous and breakTiesStochastic are FALSE, ties are broken in the order of 'A', 'T', 'G', 'C', 'N', '.', and '-'.*
```
Collapsing clonal sequences...

```


```R
# mostCommon method, resolving ties deterministically using ambiguous characters
clones5 <- collapseClones(db, method="mostCommon", 
includeAmbiguous=TRUE, breakTiesStochastic=FALSE)

```


```
Collapsing clonal sequences...

```


```R
# mostCommon method, resolving ties stochastically
clones6 <- collapseClones(db, method="mostCommon", 
includeAmbiguous=FALSE, breakTiesStochastic=TRUE)

```


```
Collapsing clonal sequences...

```


```R

# catchAll method
clones7 <- collapseClones(db, method="catchAll")

```


```
Collapsing clonal sequences...

```


```R

# mostMutated method, resolving ties stochastically
clones8 <- collapseClones(db, method="mostMutated", muFreqColumn=NULL, 
breakTiesStochastic=TRUE, breakTiesByColumns=NULL)

```

*Calculating observed mutation frequency...*
```
Collapsing clonal sequences...

```


```R
clones9 <- collapseClones(db2, method="mostMutated", muFreqColumn="MU_FREQ", 
breakTiesStochastic=TRUE, breakTiesByColumns=NULL)

```


```
Collapsing clonal sequences...

```


```R
# mostMutated method, resolving ties deterministically using additional columns
clones10 <- collapseClones(db, method="mostMutated", muFreqColumn=NULL, 
breakTiesStochastic=FALSE, breakTiesByColumns=NULL)

```

*When breakTiesStochastic is FALSE and breakTiesByColumns is NULL, ties are broken by taking the sequence that appears earlier in the data.frame.**Calculating observed mutation frequency...*
```
Collapsing clonal sequences...

```


```R
clones11 <- collapseClones(db2, method="mostMutated", muFreqColumn="MU_FREQ", 
breakTiesStochastic=FALSE, 
breakTiesByColumns=list(c("DUPCOUNT"), c(max)))

```


```
Collapsing clonal sequences...

```


```R
# mostMutated method, resolving ties deterministically without using additional columns
clones12 <- collapseClones(db, method="mostMutated", muFreqColumn=NULL, 
breakTiesStochastic=FALSE, breakTiesByColumns=NULL)

```

*When breakTiesStochastic is FALSE and breakTiesByColumns is NULL, ties are broken by taking the sequence that appears earlier in the data.frame.**Calculating observed mutation frequency...*
```
Collapsing clonal sequences...

```


```R
clones13 <- collapseClones(db2, method="mostMutated", muFreqColumn="MU_FREQ", 
breakTiesStochastic=FALSE, breakTiesByColumns=NULL)

```

*When breakTiesStochastic is FALSE and breakTiesByColumns is NULL, ties are broken by taking the sequence that appears earlier in the data.frame.*
```
Collapsing clonal sequences...

```


```R


# Build clonal consensus for V-region only
clones14 <- collapseClones(db, method="mostCommon", regionDefinition=IMGT_V)

```

*When both includeAmbiguous and breakTiesStochastic are FALSE, ties are broken in the order of 'A', 'T', 'G', 'C', 'N', '.', and '-'.*
```
Collapsing clonal sequences...

```


```R

# Return the same number of rows as the input
clones15 <- collapseClones(db, method="mostCommon", expandedDb=TRUE)
```

*When both includeAmbiguous and breakTiesStochastic are FALSE, ties are broken in the order of 'A', 'T', 'G', 'C', 'N', '.', and '-'.*
```
Collapsing clonal sequences...

```



See also
-------------------

See [IMGT_SCHEMES](IMGT_SCHEMES.md) for a set of predefined [RegionDefinition](RegionDefinition-class.md) objects.



