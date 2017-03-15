





**createRegionDefinition** - *Creates a RegionDefinition*

Description
--------------------

`createRegionDefinition` creates a `RegionDefinition`.


Usage
--------------------
```
createRegionDefinition(name = "", boundaries = factor(), description = "",
citation = "")
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

**Error in eval(expr, envir, enclos)**: could not find function "createRegionDefinition"

See also
-------------------

See [RegionDefinition](RegionDefinition-class.md) for the return object.



