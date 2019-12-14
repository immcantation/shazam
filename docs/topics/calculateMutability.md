**calculateMutability** - *Calculate total mutability*

Description
--------------------

`calculateMutability` calculates the total (summed) mutability for a set of sequences 
based on a 5-mer nucleotide mutability model.


Usage
--------------------
```
calculateMutability(sequences, model = HH_S5F, progress = FALSE)
```

Arguments
-------------------

sequences
:   character vector of sequences.

model
:   [TargetingModel](TargetingModel-class.md) object with mutation likelihood information.

progress
:   if `TRUE` print a progress bar.




Value
-------------------

Numeric vector with a total mutability score for each sequence.



Examples
-------------------

```R
# Subset example data to one isotype and sample as a demo
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, ISOTYPE == "IgA" & SAMPLE == "-1h")

# Calculate mutability of germline sequences using \link{HH_S5F} model
mutability <- calculateMutability(sequences=db$GERMLINE_IMGT_D_MASK, model=HH_S5F)
```








