**editBaseline** - *Edit the Baseline object*

Description
--------------------

`editBaseline` edits a field in a `Baseline` object.


Usage
--------------------
```
editBaseline(baseline, field, value)
```

Arguments
-------------------

baseline
:   `Baseline` object to be edited.

field
:   name of the field in the `Baseline` object to be edited.

value
:   value to set the `field`.




Value
-------------------

A `Baseline` object with the field of choice updated.



Examples
-------------------

```R
# Subset example data
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, c_call == "IGHG" & sample_id == "+7d")

# Make Baseline object
baseline <- calcBaseline(db, 
sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask", 
testStatistic="focused",
regionDefinition=IMGT_V,
targetingModel=HH_S5F,
nproc=1)

```

*calcBaseline will calculate observed and expected mutations for sequence_alignment using germline_alignment_d_mask as a reference.*
```
Calculating BASELINe probability density functions...

```


```R

# Edit the field "description"
baseline <- editBaseline(baseline, field="description", 
value="+7d IGHG")
```



See also
-------------------

See [Baseline](Baseline-class.md) for the input and return object.






