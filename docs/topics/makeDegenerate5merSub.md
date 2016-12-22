





**makeDegenerate5merSub** - *Make a degenerate 5-mer substitution model based on a 1-mer substitution model*

Description
--------------------

`makeDegenerate5merSub` populates substitution rates from a 1-mer substitution model
into 5-mers with corresponding central 1-mers.


Usage
--------------------
```
makeDegenerate5merSub(sub1mer, extended = FALSE)
```

Arguments
-------------------

sub1mer
:   a 4x4 matrix containing (normalized) substitution rates.
Row names should correspond to nucleotides to mutate from.
Column names should correspond to nucleotides to mutate into.
Nucleotides should include "A", "T", "G", and "C" 
(case-insensitive).

extended
:   whether to return the unextended (`extended=FALSE`) or 
extended (`extended=TRUE`) 5-mer substitution model. 
Default is `FALSE`.




Value
-------------------

For `extended=FALSE`, a 4x1024 matrix. For `extended=TRUE`, a 5x3125 
matrix.


Details
-------------------

As a concrete example, consider a 1-mer substitution model in which substitution
rates from "A" to "T", "G", and "C" are, respectively, 0.1, 0.6, and 0.3. In the 
resultant degenerate 5-mer substitution model, all the 5-mers (columns) that have 
an "A" as their central 1-mer would have substitution rates (rows) of 0.1, 0.6, and 
0.3 to "T", "G", and "C" respectively. 

When `extended=TRUE`, `extendSubstitutionMatrix` is called to extend
the 4x1024 substitution matrix.



Examples
-------------------

```R
# Make a degenerate 5-mer model (4x1024) based on HKL_S1F (4x4)
# Note: not to be confused with HKL_S5F@substitution, which is non-degenerate
degenerate5merSub <- makeDegenerate5merSub(sub1mer = HKL_S1F)

# Look at a few 5-mers
degenerate5merSub[, c("AAAAT", "AACAT", "AAGAT", "AATAT")]
```


```
  AAAAT AACAT AAGAT AATAT
A    NA  0.25  0.53  0.20
C  0.26    NA  0.31  0.57
G  0.50  0.30    NA  0.23
T  0.24  0.45  0.16    NA

```



See also
-------------------

See `[makeAverage1merSub](makeAverage1merSub.md)` for making a 1-mer substitution model by taking
the average of a 5-mer substitution model. See `[extendSubstitutionMatrix](extendSubstitutionMatrix.md)`
for extending the substitution matrix.



