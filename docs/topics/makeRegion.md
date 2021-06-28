**makeRegion** - *Defining a Region that will include also CDR3 and FWR4 based on junction length and sequence*

Description
--------------------

Predefined `RegionDefinition` objects don't include the CDR3 and FWR4 regions because their
boundaries depend on the junction length, and need to be calculated for each sequence. 
`makeRegion` gets as input a junction length and an IMGT aligned sequence
and outputs a `RegionDefinition` object that includes the segments CDR1/2/3 
and FWR1/2/3/4.


Usage
--------------------
```
makeRegion(juncLength, sequenceImgt, regionDefinition = NULL)
```

Arguments
-------------------

juncLength
:   The junction length of the sequence

sequenceImgt
:   The imgt aligned sequence

regionDefinition
:   The `RegionDefinition` type to calculate
the regionDefinition for. Can be one of 2: 
`"IMGT_VDJ_BY_REGIONS"` or `"IMGT_VDJ"`. 
Only these 2 regions include all
CDR1/2/3 and FWR1/2/3/4 regions.




Value
-------------------

a `RegionDefinition` object that includes CDR1/2/3 and 
FWR1/2/3/4 for the specific `sequenceImgt`, 
`juncLength` and `regionDefinition`.


Details
-------------------

For `regionDefinition="IMGT_VDJ_BY_REGIONS"` the function returns a `RegionDefinition` 
object with regions:


+ `fwr1`:   Positions 1 to 78.
+ `cdr1`:   Positions 79 to 114.
+ `fwr2`:   Positions 115 to 165.
+ `cdr2`:   Positions 166 to 195.
+ `fwr3`:   Positions 196 to 312.
+ `cdr3`:   Positions 313 to (313 + juncLength - 6) - since junction 
sequence includes (on the left) the last codon from fwr3, and 
(on the right) the first codon from fwr4.  
+ `fwr4`:   Positions (313 + juncLength - 6 + 1) to the end of the sequence.


For `regionDefinition="IMGT_VDJ"` the function returns a `RegionDefinition` 
object with regions:


+ `fwr`:   Positions belonging to a framework region.
+ `cdr`:   Positions belonging to a cdr.



Note
-------------------

In case the `regionDefinition` argument is not one of the extended
regions (`IMGT_VDJ_BY_REGIONS` or `IMGT_VDJ`) - then this
function will return the `regionDefinition` as is.



Examples
-------------------

```R
# Load and subset example data
data(ExampleDb, package="alakazam")  
juncLength <-ExampleDb[['junction_length']][1]
sequenceImgt<-ExampleDb[['sequence_alignment']][1]
seq_1_reg_def<-makeRegion(juncLength = juncLength, 
sequenceImgt = sequenceImgt, 
regionDefinition = IMGT_VDJ)
```



See also
-------------------

See [RegionDefinition](RegionDefinition-class.md) for the return object. 
See [IMGT_SCHEMES](IMGT_SCHEMES.md) for a set of predefined `RegionDefinition` objects.






