**helperClonalConsensus** - *Helper function for calcClonalConsensus*

Description
--------------------

Helper function for calcClonalConsensus


Usage
--------------------
```
helperClonalConsensus(seqs, muFreqColumn = NULL, lenLimit = NULL,
mtd = c("mostCommon", "thresholdedFreq", "catchAll", "mostMutated",
"leastMutated"), minFreq = NULL, includeAmbiguous = FALSE,
breakTiesStochastic = FALSE, breakTiesByColumns = NULL, db = NULL)
```

Arguments
-------------------

seqs
:   a character vector of sequences.

muFreqColumn
:   `character` name of the column in db containing mutation
frequency. Applicable to and required for the `"mostMutated"`
and `"leastmutated"` methods. Default is `NULL`.

lenLimit
:   limit on consensus length.

mtd
:   method to calculate consensus sequence. One of
`"thresholdedFreq"`, `"mostCommon"`, `"catchAll"`,
`"mostMutated"`, or `"leastMutated"`. See "Methods" under
[collapseClones](collapseClones.md) for details.

minFreq
:   frequency threshold for calculating input consensus sequence.
Applicable to and required for the `"thresholdedFreq"` method.
A canonical choice is 0.6. Default is `NULL`.

includeAmbiguous
:   whether to use ambiguous characters to represent positions at
which there are multiple characters with frequencies that are at least
`minimumFrequency` or that are maximal (i.e. ties). Applicable to
and required for the `"thresholdedFreq"` and `"mostCommon"`
methods. Default is `FALSE`. See "Choosing ambiguous characters"
under [collapseClones](collapseClones.md) for rules on choosing ambiguous characters.

breakTiesStochastic
:   In case of ties, whether to randomly pick a sequence from sequences that
fulfill the criteria as consensus. Applicable to and required for all methods
except for `"catchAll"`. Default is `FALSE`. See "Methods"
under [collapseClones](collapseClones.md) for details.

breakTiesByColumns
:   A list of the form `list(c(col_1, col_2, ...), c(fun_1, fun_2, ...))`,
where `col_i` is a `character` name of a column in `db`,
and `fun_i` is a function to be applied on that column. Currently,
only `max` and `min` are supported. Note that the two `c()`'s
in `list()` are essential (i.e. if there is only 1 column, the list
should be of the form `list(c(col_1), c(func_1))`. Applicable to and
optional for the `"mostMutated"` and `"leastMutated"` methods.
If supplied, `fun_i`'s are applied on `col_i`'s to help break
ties. Default is `NULL`. See "Methods" under [collapseClones](collapseClones.md)
for details.

db
:   `data.frame` containing sequence data for a single clone.
Applicable to and required for the `"mostMutated"` and
`"leastmutated"` methods. Default is `NULL`.




Value
-------------------

A list containing `cons`, which is a character string that is the consensus sequence
for `seqs`; and `muFreq`, which is the maximal/minimal mutation frequency of
the consensus sequence for the `"mostMutated"` and `"leastMutated"` methods, or
`NULL` for all other methods.


Details
-------------------

See [collapseClones](collapseClones.md) for detailed documentation on methods and additional parameters.



Examples
-------------------

```R
# Subset example data
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, ISOTYPE %in% c("IgA", "IgG") & SAMPLE == "+7d")
# Data corresponding to a single clone
clone <- db[db$CLONE=="3192", ]
# Number of sequences in this clone
nrow(clone)

```


```
[1] 19

```


```R
# first compute mutation frequency for most/leastMutated methods
clone = observedMutations(db=clone, frequency=TRUE, combine=TRUE)
# manually create tie
clone = rbind(clone, clone[which.max(clone$MU_FREQ), ])
# Get consensus input sequence
# thresholdedFreq method, resolve ties deterministically without using ambiguous characters
consInput1 <- helperClonalConsensus(seqs=clone$SEQUENCE_IMGT,
muFreqColumn=NULL, lenLimit=NULL,
mtd="thresholdedFreq", minFreq=0.3,
includeAmbiguous=FALSE, 
breakTiesStochastic=FALSE,
breakTiesByColumns=NULL, db=NULL)$cons
```




