**MutabilityModel-class** - *S4 class defining a mutability model*

Description
--------------------

`MutabilityModel` defines a data structure for the 5-mer motif-based SHM targeting
mutability model.


Usage
--------------------
```
"print"(x)
```
```
"as.data.frame"(x)
```

Arguments
-------------------

x
:   `MutabilityModel` object.




Slots
-------------------



`.Data`
:   numeric vector containing 5-mer mutability estimates

`source`
:   character vector annotating whether the mutability was
inferred or directly measured.

`numMutS`
:   a number indicating the number of silent mutations used for 
estimating mutability

`numMutR`
:   a number indicating the number of replacement mutations used 
for estimating mutability









