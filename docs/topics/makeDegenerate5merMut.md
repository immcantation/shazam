**makeDegenerate5merMut** - *Make a degenerate 5-mer mutability model based on a 1-mer mutability model*

Description
--------------------

`makeDegenerate5merMut` populates mutability rates from a 1-mer mutability model
into 5-mers with corresponding central 1-mers.


Usage
--------------------
```
makeDegenerate5merMut(mut1mer, extended = FALSE)
```

Arguments
-------------------

mut1mer
:   a named vector of length 4 containing (normalized) 
mutability rates. Names should correspond to nucleotides, 
which should include "A", "T", "G", and "C" 
(case-insensitive).

extended
:   whether to return the unextended (`extended=FALSE`) or 
extended (`extended=TRUE`) 5-mer mutability model. 
Default is `FALSE`.




Value
-------------------

For `extended=FALSE`, a vector of length 1024. The vector returned is 
normalized. For `extended=TRUE`, a vector of length 3125.


Details
-------------------

As a concrete example, consider a 1-mer mutability model in which mutability
rates of "A", "T", "G", and "C" are, respectively, 0.14, 0.23, 0.31, and 0.32. 
In the resultant degenerate 5-mer mutability model, all the 5-mers that have 
an "A" as their central 1-mer would have mutability rate of 0.14/256, where 256 is
the number of such 5-mers. 

When `extended=TRUE`, `extendMutabilityMatrix` is called to extend the
mutability vector of length 1024 into a vector of length 3125.



Examples
-------------------

```R
# Make a degenerate 5-mer model (length of 1024) based on a 1-mer model
example1merMut <- c(A=0.2, T=0.1, C=0.4, G=0.3)
degenerate5merMut <- makeDegenerate5merMut(mut1mer = example1merMut)

# Look at a few 5-mers
degenerate5merMut[c("AAAAT", "AACAT", "AAGAT", "AATAT")]

```


```
      AAAAT       AACAT       AAGAT       AATAT 
0.000781250 0.001562500 0.001171875 0.000390625 

```


```R

# Normalized
sum(degenerate5merMut)
```


```
[1] 1

```



See also
-------------------

See [makeAverage1merMut](makeAverage1merMut.md) for making a 1-mer mutability model by 
taking the average of a 5-mer mutability model. See 
[extendMutabilityMatrix](extendMutabilityMatrix.md) for extending the mutability vector.



