





**TargetingModel-class** - *S4 class defining a targeting model*

Description
--------------------

`TargetingModel` defines a common data structure for mutability, substitution and
targeting of immunoglobulin (Ig) sequencing data in a 5-mer microsequence context.


Usage
--------------------
```
"plot"(x, y, ...)
```

Arguments
-------------------

y
:   ignored.

...
:   arguments to pass to [plotMutability](plotMutability.md).




Slots
-------------------



`name`
:   Name of the model.

`description`
:   Description of the model and its source data.

`species`
:   Genus and species of the source sequencing data.

`date`
:   Date the model was built.

`citation`
:   Publication source.

`substitution`
:   Normalized rates of the center nucleotide of a given 5-mer 
mutating to a different nucleotide. The substitution model 
is stored as a 5x3125 matrix of rates. Rows define
the mutated nucleotide at the center of each 5-mer, one of 
`c("A", "C", "G", "T", "N")`, and columns define the 
complete 5-mer of the unmutated nucleotide sequence.

`mutability`
:   Normalized rates of a given 5-mer being mutated. The 
mutability model is stored as a numeric vector of length 3125 
with mutability rates for each 5-mer. Note that "normalized" 
means that the mutability rates for the 1024 5-mers that 
contain no "N" at any position sums up to 1 (as opposed to 
the entire vector summing up to 1).

`targeting`
:   Rate matrix of a given mutation ocurring, defined as 
<code class = 'eq'>mutability * substitution</code>. The targeting model 
is stored as a 5x3125 matrix. Rows define
the mutated nucleotide at the center of each 5-mer, one of 
`c("A", "C", "G", "T", "N")`, and columns define the complete 5-mer 
of the unmutated nucleotide sequence.




See also
-------------------

See [createTargetingModel](createTargetingModel.md) building models from sequencing data.



