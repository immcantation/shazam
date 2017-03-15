





**makeAverage1merMut** - *Make a 1-mer mutability model by averaging over a 5-mer mutability model*

Description
--------------------

`makeAverage1merMut` averages mutability rates in a 5-mer mutability model
to derive a 1-mer mutability model.


Usage
--------------------
```
makeAverage1merMut(mut5mer)
```

Arguments
-------------------

mut5mer
:   a named vector of length 1024 such as that returned by 
`createMutabilityMatrix` and that returned by
`makeDegenerate5merMut` with `extended=FALSE`.
Names should correspond to 5-mers made up of "A", "T", 
"G", and "C" (case-insensitive). `NA` values are 
allowed.




Value
-------------------

A named vector of length 4 containing normalized mutability rates.


Details
-------------------

For example, the mutability rate of "A" in the resultant 1-mer model
is derived by averaging the mutability rates of all the 5-mers that 
have an "A" as their central 1-mer, followed by normalization.



Examples
-------------------

```R
# Make a degenerate 5-mer model (length of 1024) based on a 1-mer model
example1merMut <- c(A=0.2, T=0.1, C=0.4, G=0.3)
degenerate5merMut <- makeDegenerate5merMut(mut1mer = example1merMut)

```

**Error in eval(expr, envir, enclos)**: could not find function "makeDegenerate5merMut"
```R
 
# Now make a 1-mer model by averaging over the degenerate 5-mer model
# Expected to get back example1merMut
makeAverage1merMut(mut5mer = degenerate5merMut)
```

**Error in eval(expr, envir, enclos)**: could not find function "makeAverage1merMut"

See also
-------------------

See [makeDegenerate5merMut](makeDegenerate5merMut.md) for making a degenerate 5-mer mutability
model based on a 1-mer mutability model.



