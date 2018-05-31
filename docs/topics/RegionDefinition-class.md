**RegionDefinition-class** - *S4 class defining a region definition*

Description
--------------------

`RegionDefinition` defines a common data structure for defining the region
boundaries of an Ig sequence.




Slots
-------------------



`name`
:   name of the RegionDefinition.

`description`
:   description of the model and its source.

`boundaries`
:   `factor` defining the region boundaries of the 
sequence. The levels and values of `boundaries` 
determine the number of regions.

`seqLength`
:   length of the sequence.

`regions`
:   levels of the boundaries; e.g, `c("CDR", "FWR")`.

`labels`
:   labels for the boundary and mutations combinations;
e.g., `c("CDR_R", "CDR_S", "FWR_R", "FWR_S")`.

`citation`
:   publication source.




See also
-------------------

See [IMGT_SCHEMES](IMGT_SCHEMES.md) for a set of predefined `RegionDefinition` objects.



