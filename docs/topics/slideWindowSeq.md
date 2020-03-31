**slideWindowSeq** - *Sliding window approach towards filtering a single sequence*

Description
--------------------

`slideWindowSeq` determines whether an input sequence contains equal to or more than 
a given number of mutations in a given length of consecutive nucleotides (a "window") 
when compared to a germline sequence.


Usage
--------------------
```
slideWindowSeq(inputSeq, germlineSeq, mutThresh, windowSize)
```

Arguments
-------------------

inputSeq
:   input sequence.

germlineSeq
:   germline sequence.

mutThresh
:   threshold on the number of mutations in `windowSize` 
consecutive nucleotides. Must be between 1 and `windowSize` 
inclusive.

windowSize
:   length of consecutive nucleotides. Must be at least 2.




Value
-------------------

`TRUE` if there are equal to or more than `mutThresh` number of mutations
in any window of `windowSize` consecutive nucleotides (i.e. the sequence should
be filtered); `FALSE` if otherwise.



Examples
-------------------

```R
# Use an entry in the example data for input and germline sequence
data(ExampleDb, package="alakazam")
in_seq <- ExampleDb[["sequence_alignment"]][100]
germ_seq <-  ExampleDb[["germline_alignment_d_mask"]][100]

# Determine if in_seq has 6 or more mutations in 10 consecutive nucleotides
slideWindowSeq(inputSeq=in_seq, germlineSeq=germ_seq, mutThresh=6, windowSize=10)
```


```
[1] FALSE

```



See also
-------------------

[calcObservedMutations](calcObservedMutations.md) is called by `slideWindowSeq` to identify observed 
mutations. See [slideWindowDb](slideWindowDb.md) for applying the sliding window approach on a 
`data.frame`. See [slideWindowTune](slideWindowTune.md) for parameter tuning for `mutThresh`
and `windowSize`.






