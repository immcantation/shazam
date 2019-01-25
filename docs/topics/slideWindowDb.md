**slideWindowDb** - *Sliding window approach towards filtering sequences in a `data.frame`*

Description
--------------------

`slideWindowDb` determines whether each input sequence in a `data.frame` 
contains equal to or more than a given number of mutations in a given length of 
consecutive nucleotides (a "window") when compared to their respective germline 
sequence.


Usage
--------------------
```
slideWindowDb(db, sequenceColumn = "SEQUENCE_IMGT",
germlineColumn = "GERMLINE_IMGT_D_MASK", mutThresh, windowSize)
```

Arguments
-------------------

db
:   `data.frame` containing sequence data.

sequenceColumn
:   name of the column containing IMGT-gapped sample sequences.

germlineColumn
:   name of the column containing IMGT-gapped germline sequences.

mutThresh
:   threshold on the number of mutations in `windowSize` 
consecutive nucleotides. Must be between 1 and `windowSize` 
inclusive.

windowSize
:   length of consecutive nucleotides. Must be at least 2.




Value
-------------------

a logical vector. The length of the vector matches the number of input sequences in 
`db`. Each entry in the vector indicates whether the corresponding input sequence
should be filtered based on the given parameters.



Examples
-------------------

```R
# Use an entry in the example data for input and germline sequence
data(ExampleDb, package="alakazam")

# Apply the sliding window approach on a subset of ExampleDb
slideWindowDb(db = ExampleDb[1:10, ], mutThresh=6, windowSize=10)
```


```
 [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE

```



See also
-------------------

See [slideWindowSeq](slideWindowSeq.md) for applying the sliding window approach on a single sequence. 
See [slideWindowTune](slideWindowTune.md) for parameter tuning for `mutThresh` and `windowSize`.



