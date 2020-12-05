**makeRegion** - *Defining a Region that will include also CDR3 and FWR4 based on junction length and sequence*

Description
--------------------

Defining a Region that will include also CDR3 and FWR4 based on junction length and sequence


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
:   The [RegionDefinition](RegionDefinition-class.md) type to calculate
the regionDefinition for. Can be one of 2: 
`"IMGT_VDJ_BY_REGIONS"` or `"IMGT_VDJ"`. 
Only these 2 regions include all
CDR1/2/3 and FWR1/2/3/4 regions.




Value
-------------------

a [RegionDefinition](RegionDefinition-class.md) object that includes CDR1/2/3 and 
FWR1/2/3/4 for the specific `sequenceImgt`, 
`juncLength` and `regionDefinition`.


Details
-------------------

This function gets as input a junction length and an imgt aligned sequence
and outputs a [RegionDefinition](RegionDefinition-class.md) object that includes following regions:   

**For `regionDefinition="IMGT_VDJ_BY_REGIONS"`:**

- **fwr1**: Bases 1 to 78 (based on [IMGT_V_BY_REGIONS](IMGT_SCHEMES.md) definitions)  

- **cdr1**: Bases 79 to 114 (based on [IMGT_V_BY_REGIONS](IMGT_SCHEMES.md) definitions) 

- **fwr2**: Bases 115 to 165 (based on [IMGT_V_BY_REGIONS](IMGT_SCHEMES.md) definitions) 

- **cdr2**: Bases 166 to 195 (based on [IMGT_V_BY_REGIONS](IMGT_SCHEMES.md) definitions) 

- **fwr3**: Bases 196 to 312 (based on [IMGT_V_BY_REGIONS](IMGT_SCHEMES.md) definitions)

- **cdr3**: Bases 313 to (313 + `juncLength` - 6) - since junction sequnece 
includes (on the left) the last codon from fwr3, and (on the right)  
the first codon from fwr4.  

- **fwr4**: Bases (313 + `juncLength` - 6 + 1) to sequence_length. 

**For `regionDefinition`="IMGT_VDJ":**

- **fwr**: Bases	1 to 78 (based on [IMGT_V_BY_REGIONS](IMGT_SCHEMES.md) definitions)

Bases	115 to 165 (based on [IMGT_V_BY_REGIONS](IMGT_SCHEMES.md) definitions)
 
Bases	196 to 312 (based on [IMGT_V_BY_REGIONS](IMGT_SCHEMES.md) definitions)
 
Bases	(313 + `juncLength` - 6 + 1) to sequence_length.

- **cdr**: Bases	79 to 114 (based on [IMGT_V_BY_REGIONS](IMGT_SCHEMES.md) definitions)
 
Bases	166 to 195 (based on [IMGT_V_BY_REGIONS](IMGT_SCHEMES.md) definitions) 

Bases	313 to (313 + `juncLength` - 6) - since junction sequnece 
includes (on the left) the last codon from fwr3, and (on the right) 
the first codon from fwr4.  

Note: In case the `regionDefinition` argument is not one of the extended
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
regionDefinition = IMGT_VDJ_BY_REGIONS)
```








