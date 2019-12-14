**createRegionDefinition** - *Creates a RegionDefinition*

Description
--------------------

`createRegionDefinition` creates a `RegionDefinition`.


Usage
--------------------
```
createRegionDefinition(name = "", boundaries = factor(),
description = "", citation = "")
```

Arguments
-------------------

name
:   name of the region definition.

boundaries
:   `factor` defining the region boundaries of the sequence.
The levels and values of `boundaries` determine the 
number of regions (e.g. CDR and FWR).

description
:   description of the region definition and its source data.

citation
:   publication source.




Value
-------------------

A `RegionDefinition` object.



Examples
-------------------

```R
# Creates an empty RegionDefinition object
createRegionDefinition()
```


```
An object of class "RegionDefinition"
Slot "name":
[1] ""

Slot "description":
[1] ""

Slot "boundaries":
factor(0)
Levels: 

Slot "seqLength":
[1] 0

Slot "regions":
character(0)

Slot "labels":
character(0)

Slot "citation":
[1] ""


```



See also
-------------------

See [RegionDefinition](RegionDefinition-class.md) for the return object.






