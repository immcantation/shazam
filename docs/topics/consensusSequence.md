**consensusSequence** - *Construct a consensus sequence*

Description
--------------------

Construct a consensus sequence


Usage
--------------------
```
consensusSequence(sequences, db = NULL, method = c("mostCommon",
"thresholdedFreq", "catchAll", "mostMutated", "leastMutated"),
minFreq = NULL, muFreqColumn = NULL, lenLimit = NULL,
includeAmbiguous = FALSE, breakTiesStochastic = FALSE,
breakTiesByColumns = NULL)
```

Arguments
-------------------

sequences
:   character vector of sequences.

db
:   `data.frame` containing sequence data for a single clone.
Applicable to and required for the `"mostMutated"` and
`"leastMutated"` methods. Default is `NULL`.

method
:   method to calculate consensus sequence. One of
`"thresholdedFreq"`, `"mostCommon"`, `"catchAll"`,
`"mostMutated"`, or `"leastMutated"`. See "Methods" under
[collapseClones](collapseClones.md) for details.

minFreq
:   frequency threshold for calculating input consensus sequence.
Applicable to and required for the `"thresholdedFreq"` method.
A canonical choice is 0.6. Default is `NULL`.

muFreqColumn
:   `character` name of the column in db containing mutation
frequency. Applicable to and required for the `"mostMutated"`
and `"leastMutated"` methods. Default is `NULL`.

lenLimit
:   limit on consensus length. if `NULL` then no length limit is set.

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




Value
-------------------

A list containing `cons`, which is a character string that is the consensus sequence
for `sequences`; and `muFreq`, which is the maximal/minimal mutation frequency of
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
clone <- subset(db, CLONE == "3192")

# First compute mutation frequency for most/leastMutated methods
clone <- observedMutations(clone, frequency=TRUE, combine=TRUE)

# Manually create a tie
clone <- rbind(clone, clone[which.max(clone$MU_FREQ), ])

# ThresholdedFreq method. 
# Resolve ties deterministically without using ambiguous characters
cons1 <- consensusSequence(clone$SEQUENCE_IMGT,
method="thresholdedFreq", minFreq=0.3,
includeAmbiguous=FALSE, 
breakTiesStochastic=FALSE)
cons1$cons
```


```
[1] "GAGGTGCAGCTGGTGGTCTCTGGGGGA...GGCTTGGTACAGCCAGGGCGGTCCCTAAGACTCTCCTGTACAGTTTCTGGATTCACCTTT............GGTGATTATGCTATGACGTGGATCCGCCAGGCTCCTGGGAAGGGGCTGGAGTGGGTCGGTTTCATTAGAAGCAAAACTTTTGGTGGGACAGCAGATTACGCCGCGTTTGTGAGA...GGCAGATTCACCATCTCAAGAGATGATTCCAAAAACATCGCCTATCTGCAATTGAACAGCCTGAAAACCGAGGACACAGGCGTCTATTACTGTGGTAGAGATCTCGCCGTAACTGACACAATAGGTGGTACTAACTGGTTCGACCCCTGGGGCCAGGGGACCCCGGTCACCGTCTCCTCAG"

```








