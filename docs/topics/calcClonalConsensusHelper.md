





**calcClonalConsensusHelper** - *Helper function for calcClonalConsensus*

Description
--------------------

Helper function for calcClonalConsensus


Usage
--------------------
```
calcClonalConsensusHelper(seqs, minFreq = 0.6, lenLimit = NULL,
mtd = c("thresholdedFreq", "mostCommon", "catchAll"))
```

Arguments
-------------------

seqs
:   a character vector of sequences.

minFreq
:   minimum frequency. Required if `mtd` is `thresholdedFreq`.
Default is 0.6.

lenLimit
:   limit on consensus length. Length of consensus returned will depend on
`lenLimit` and the longest possible length based on `seqs`, 
whichever that is shorter.

mtd
:   method to calculate consensus sequence. One of `thresholdedFreq`, 
`mostCommon`, or `catchAll`.




Value
-------------------

A character string that is the consensus sequence for `seqs`.


Details
-------------------

Note that this function does not perform multiple sequence alignment. As a 
prerequisite, it is assumed that sequences in `seqs` have been aligned
somehow. In the case of immunoglobulin repertoire analysis, this usually means
that the sequences are IMGT-gapped.

Given a set of sequences of potentially varying lengths, the longest possible 
length of their consensus sequence is taken to be the longest length along 
which there is information contained at every nucleotide position across 
majority of the sequences. Majority is defined to be greater than 
`floor(n/2)`, where `n` is the number of sequences. If the longest 
possible consensus length is 0, there will be a warning and an empty string 
(`""`) will be returned. 

If a length limit is defined by supplying a `lenLimit`, the consensus 
length will be further restricted to the shorter of the longest possible 
length and `lenLimit`.

For the method `"thresholdedFreq"`, a value must be supplied to the 
argument `minFreq`. At each position along the length of the consensus 
sequence, the frequency of each nucleotide/character across sequences is 
tabulated. The nucleotide/character whose frequency is at least (i.e. 
`>=`) `minFreq` becomes the consensus. If frequencies of multiple 
nucleotides/characters are at least `minFreq`, the first one is taken 
to be the consensus following the order of `"A"`, `"T"`, `"G"`, 
`"C"`, `"N"`, `"."`, and `"-"`. Here are some examples 
looking at a single position based on 5 sequences with `minFreq=0.6`:


+  If the sequences have `"A"`, `"A"`, `"A"`, `"T"`, 
`"C"`, the consensus will be `"A"`, because `"A"` has frequency 0.6, which is at least `minFreq`.
+  If the sequences have `"A"`, `"A"`, `"T"`, `"T"`, 
`"C"`, the consensus will be `"N"`, because none of 
`"A"`, `"T"`, or `"C"` has frequency that is at least 
`minFreq`.


For the method `"mostCommon"`, the most frequent nucleotide/character 
across sequences at each position along the length of the consensus sequence 
makes up the consensus. If there are multiple nucleotides/characters with 
equally maximal frequencies, the first one is taken to be the consensus 
following the order of `"A"`, `"T"`, `"G"`, `"C"`, 
`"N"`, `"."`, and `"-"`. Here are some examples looking at 
a single position based on 5 sequences:


+  If the sequences have `"A"`, `"A"`, `"T"`, `"A"`,
`"C"`, the consensus will be `"A"`.
+  If the sequences have `"T"`, `"T"`, `"C"`, `"C"`, 
`"G"`, the consensus will be `"T"`, because `"T"` is 
before `"C"` in the order of `"A"`, `"T"`, `"G"`, 
`"C"`, `"N"`, `"."`, and `"-"`. 


For the method `"catchAll"`, a consensus sequence capturing most of the 
information contained in the sequences is obtained. If a position contains only 
`"N"` across sequences, the consensus at that position is `"N"`. If a 
position contains one or more of `"A"`, `"T"`, `"G"`, or 
`"C"`, the consensus will be an IUPAC character representing all of the 
characters present, regardless of whether `"N"` is present. If a position 
contains only `"-"` and `"."` across sequences, the consensus at that 
position is taken to be `"-"`. If a position contains only one of 
`"-"` or `"."` across sequences, the consensus at that position is 
taken to be the character present. Here are some examples looking at a single 
position based on 5 sequences:


+  If the sequences have `"N"`, `"N"`, `"N"`, `"N"`, 
`"N"`, the consensus will be `"N"`.
+  If the sequences have `"N"`, `"A"`, `"A"`, `"A"`, 
`"A"`, the consensus will be `"A"`.
+  If the sequences have `"N"`, `"A"`, `"G"`, `"A"`, 
`"A"`, the consensus will be `"R"`.
+  If the sequences have `"-"`, `"-"`, `"."`, `"."`, 
`"."`, the consensus will be `"-"`.
+  If the sequences have `"-"`, `"-"`, `"-"`, `"-"`, 
`"-"`, the consensus will be `"-"`.
+  If the sequences have `"."`, `"."`, `"."`, `"."`, 
`"."`, the consensus will be `"."`.




Examples
-------------------

```R
# Subset example data
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")

# Data corresponding to a single clone
clone <- db[db$CLONE=="3192", ]
# Number of sequences in this clone
nrow(clone)

```


```
[1] 19

```


```R

# Get consensus input sequence
consInput <- calcClonalConsensusHelper(seqs=clone$SEQUENCE_IMGT, 
mtd="thresholdedFreq", minFreq=0.65)
```



See also
-------------------

[calcClonalConsensus](calcClonalConsensus.md) and [collapseClones](collapseClones.md).



