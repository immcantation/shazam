





**createMutationDefinition** - *Creates a MutationDefinition*

Description
--------------------

`createMutationDefinition` creates a `MutationDefinition`.


Usage
--------------------
```
createMutationDefinition(name, classes, description = "", citation = "")
```

Arguments
-------------------

name
:   name of the mutation definition.

classes
:   named character vectors with single-letter amino acid codes as names
and amino acid classes as values, with `NA` assigned to set of 
characters `c("X", "*", "-", ".")`. Replacement (R) is be 
defined as a change in amino acid class and silent (S) as no 
change in class.

description
:   description of the mutation definition and its source data.

citation
:   publication source.



Value
-------------------

A `MutationDefinition` object.



Examples
-------------------

```R
# Define hydropathy classes
library(alakazam)
hydropathy <- list(hydrophobic=c("A", "I", "L", "M", "F", "W", "V"),
hydrophilic=c("R", "N", "D", "C", "Q", "E", "K"),
neutral=c("G", "H", "P", "S", "T", "Y"))
chars <- unlist(hydropathy, use.names=FALSE)
classes <- setNames(translateStrings(chars, hydropathy), chars)

# Create hydropathy mutation definition
md <- createMutationDefinition("Hydropathy", classes)
```



See also
-------------------

See [MutationDefinition](MutationDefinition-class.md) for the return object.



