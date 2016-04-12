





**MutationDefinition-class** - *S4 class defining replacement and silent mutation definitions*

Description
--------------------

`MutationDefinition` defines a common data structure for defining the whether
a mutation is annotated as a replacement or silent mutation.



Slots
-------------------



`name`
:   name of the MutationDefinition.

`description`
:   description of the model and its source.

`classes`
:   named character vectors with single-letter amino acid codes as names
and amino acid classes as values, with `NA` assigned to set of 
characters `c("X", "*", "-", ".")`. Replacement (R) is be 
defined as a change in amino acid class and silent (S) as no 
change in class.

`codonTable`
:   matrix of codons (columns) and substitutions (rows).

`citation`
:   publication source.




See also
-------------------

See [MUTATION_SCHEMES](MUTATION_SCHEMES.md) for a set of predefined `MutationDefinition` objects.



