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

```

*

To cite the alakazam package in publications, please use:

  Gupta N, Vander Heiden J, Uduman M, Gadala-Maria D, Yaari G,
  Kleinstein S (2015). “Change-O: a toolkit for analyzing large-scale B
  cell immunoglobulin repertoire sequencing data.” _Bioinformatics_,
  1-3. doi:10.1093/bioinformatics/btv359
  <https://doi.org/10.1093/bioinformatics/btv359>.

To cite the Ig-specific lineage reconstruction and diversity methods,
please use:

  Stern J, Yaari G, Vander Heiden J, Church G, Donahue W, Hintzen R,
  Huttner A, Laman J, Nagra R, Nylander A, Pitt D, Ramanan S, Siddiqui
  B, Vigneault F, Kleinstein S, Hafler D, O'Connor K (2014). “B cells
  populating the multiple sclerosis brain mature in the draining
  cervical lymph nodes.” _Science Translational Medicine_, *6*(248),
  248ra107. doi:10.1126/scitranslmed.3008879
  <https://doi.org/10.1126/scitranslmed.3008879>.

To see these entries in BibTeX format, use 'format(<citation>,
bibtex=TRUE)', or 'toBibtex(.)'.*
```R
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






