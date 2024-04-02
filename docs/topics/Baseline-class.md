**Baseline-class** - *S4 class defining a BASELINe (selection) object*

Description
--------------------

`Baseline` defines a common data structure the results of selection
analysis using the BASELINe method.


Usage
--------------------
```
"plot"(x, y, ...)
```
```
"summary"(object, nproc = 1)
```

Arguments
-------------------

x
:   `Baseline` object.

y
:   name of the column in the `db` slot of `baseline` 
containing primary identifiers.

...
:   arguments to pass to [plotBaselineDensity](plotBaselineDensity.md).

object
:   `Baseline` object.

nproc
:   number of cores to distribute the operation over.




Slots
-------------------



`description`
:   `character` providing general information regarding the 
sequences, selection analysis and/or object.

`db`
:   `data.frame` containing annotation information about 
the sequences and selection results.

`regionDefinition`
:   [RegionDefinition](RegionDefinition-class.md) object defining the regions
and boundaries of the Ig sequences.

`testStatistic`
:   `character` indicating the statistical framework 
used to test for selection. For example, `"local"` or 
`"focused"`.

`regions`
:   `character` vector defining the regions the BASELINe 
analysis was carried out on. For `"cdr"` and `"fwr"` 
or `"cdr1"`, `"cdr2"`, `"cdr3"`, etc.

`numbOfSeqs`
:   `matrix` of dimensions `r x c` containing the number of 
sequences or PDFs in each region, where:
`r` = number of rows = number of groups or sequences.
`c` = number of columns = number of regions.

`binomK`
:   `matrix` of dimensions `r x c` containing the number of 
successes in the binomial trials in each region, where:
`r` = number of rows = number of groups or sequences.
`c` = number of columns = number of regions.

`binomN`
:   `matrix` of dimensions `r x c` containing the total 
number of trials in the binomial in each region, where:
`r` = number of rows = number of groups or sequences.
`c` = number of columns = number of regions.

`binomP`
:   `matrix` of dimensions `r x c` containing the probability 
of success in one binomial trial in each region, where:
`r` = number of rows = number of groups or sequences.
`c` = number of columns = number of regions.

`pdfs`
:   `list` of matrices containing PDFs with one item for each 
defined region (e.g. `cdr` and `fwr`). Matrices have dimensions
`r x c` dimensions, where:
`r` = number of rows = number of sequences or groups. 
`c` = number of columns = length of the PDF (default 4001).

`stats`
:   `data.frame` of BASELINe statistics, 
including: mean selection strength (mean Sigma), 95% confidence 
intervals, and p-values with positive signs for the presence of 
positive selection and/or p-values with negative signs for the
presence of negative selection.




See also
-------------------

See [summarizeBaseline](summarizeBaseline.md) for more information on `@stats`.






