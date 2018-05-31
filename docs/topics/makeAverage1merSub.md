**makeAverage1merSub** - *Make a 1-mer substitution model by averaging over a 5-mer substitution model*

Description
--------------------

`makeAverage1merSub` averages substitution rates in a 5-mer substitution model
to derive a 1-mer substitution model.


Usage
--------------------
```
makeAverage1merSub(sub5mer)
```

Arguments
-------------------

sub5mer
:   a 4x1024 matrix such as that returned by 
`createSubstitutionMatrix` and that returned by
`makeDegenerate5merSub` with `extended=FALSE`.
Column names should correspond to 5-mers containing the 
central 1-mer to mutate from. Row names should correspond to 
nucleotides to mutate into. Nucleotides should include 
"A", "T", "G", and "C" (case-insensitive).




Value
-------------------

A 4x4 matrix with row names representing nucleotides to mutate from and column
names representing nucleotides to mutate into. Rates are normalized by row.


Details
-------------------

For example, the substitution rate from "A" to "T" in the resultant 1-mer model
is derived by averaging the substitution rates into a "T" of all the 5-mers that 
have an "A" as their central 1-mer.



Examples
-------------------

```R
# Make a degenerate 5-mer model (4x1024) based on HKL_S1F (4x4)
degenerate5merSub <- makeDegenerate5merSub(sub1mer = HKL_S1F)

# Now make a 1-mer model by averaging over the degenerate 5-mer model
# Expected to get back HKL_S1F
makeAverage1merSub(sub5mer = degenerate5merSub)
```


```
     A    C    G    T
A   NA 0.26 0.50 0.24
C 0.25   NA 0.30 0.45
G 0.53 0.31   NA 0.16
T 0.20 0.57 0.23   NA

```



See also
-------------------

See [makeDegenerate5merSub](makeDegenerate5merSub.md) for making a degenerate 5-mer substitution 
model based on a 1-mer substitution model.



