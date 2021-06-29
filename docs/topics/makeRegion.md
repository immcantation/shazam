**makeRegion** - *Build a RegionDefinition object that includes CDR3 and FWR4.*

Description
--------------------

`makeRegion` takes as input a junction length and an IMGT-numbered sequence
and outputs a custom `RegionDefinition` object that includes the boundary definitions of 
CDR1-3 and FWR1-4 for that sequence. In contrast the universal `RegionDefinition` object 
that end with FWR3, the returned definition is per-sequence due to variable junction lengths.


Usage
--------------------
```
makeRegion(juncLength, sequenceImgt, regionDefinition = NULL)
```

Arguments
-------------------

juncLength
:   junction length of the sequence.

sequenceImgt
:   IMGT-numbered sequence.

regionDefinition
:   `RegionDefinition` type to calculate the region definition for. 
Can be one of `IMGT_VDJ_BY_REGIONS` or `IMGT_VDJ`,
which are template definitions that include CDR1-3 and FWR1-4. 
Only these two regions include all CDR1-3 and FWR1-4 regions.
If this argument is set to `NULL` then an empty 
`RegionDefinition` will be returned.




Value
-------------------

A `RegionDefinition` object that includes CDR1-3 and FWR1-4 for the  
`sequenceImgt`, `juncLength`, and `regionDefinition` specified.

For `regionDefinition=IMGT_VDJ_BY_REGIONS` the returned `RegionDefinition` 
includes:


+ `fwr1`:   Positions 1 to 78.
+ `cdr1`:   Positions 79 to 114.
+ `fwr2`:   Positions 115 to 165.
+ `cdr2`:   Positions 166 to 195.
+ `fwr3`:   Positions 196 to 312.
+ `cdr3`:   Positions 313 to (313 + juncLength - 6) - since junction 
sequence includes (on the left) the last codon from FWR3, and 
(on the right) the first codon from FWR4.  
+ `fwr4`:   Positions (313 + juncLength - 6 + 1) to the end of the sequence.


For `regionDefinition=IMGT_VDJ` the returned `RegionDefinition` includes:


+ `fwr`:   Positions belonging to a FWR.
+ `cdr`:   Positions belonging to a CDR.



Note
-------------------

In case the `regionDefinition` argument is not one of the extended
regions (`IMGT_VDJ_BY_REGIONS` or `IMGT_VDJ`) - then this
function will return the input `regionDefinition` as is.



Examples
-------------------

```R
# Load and subset example data
data(ExampleDb, package="alakazam")  
len <- ExampleDb$junction_length[1]
sequence <- ExampleDb$sequence_alignment[1]
region <- makeRegion(len, sequence, regionDefinition=IMGT_VDJ)
```



See also
-------------------

See [RegionDefinition](RegionDefinition-class.md) for the return object. 
See [IMGT_SCHEMES](IMGT_SCHEMES.md) for a set of predefined `RegionDefinition` objects.






