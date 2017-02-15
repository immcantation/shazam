





**editBaseline** - *Edit the Baseline object*

Description
--------------------

`editBaseline` edits a field in a `Baseline` object.


Usage
--------------------
```
editBaseline(baseline, field_name, value)
```

Arguments
-------------------

baseline
:   The `Baseline` S4 object to be edited.

field_name
:   Name of the field in the `Baseline` S4 object to be edited.

value
:   The value to set the `field_name`.




Value
-------------------

A `Baseline` object with the field of choice updated.



Examples
-------------------

```R
# Subset example data
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")

# Calculate BASELINe
baseline <- calcBaseline(db, 
sequenceColumn="SEQUENCE_IMGT",
germlineColumn="GERMLINE_IMGT_D_MASK", 
testStatistic="focused",
regionDefinition=IMGT_V,
targetingModel=HH_S5F,
nproc = 1)

```


```
Collapsing clonal sequences...
Calculating observed number of mutations...
Calculating the expected frequencies of mutations...
Calculating BASELINe probability density functions...

```


```R
# Edit the field "description"
baseline <- editBaseline(baseline, field_name = "description", 
value = "+7d IgA & IgG")
```



See also
-------------------

See [Baseline](Baseline-class.md) for the return object.



