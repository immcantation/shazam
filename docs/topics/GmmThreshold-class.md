





**GmmThreshold-class** - *Output of the `gmm` method of findThreshold*

Description
--------------------

`GmmThreshold` contains output from the `gmm` method [findThreshold](findThreshold.md). 
It includes parameters of two Gaussian fits and threshold cut.


Usage
--------------------
```
"print"(x)
```
```
"plot"(x, y, ...)
```

Arguments
-------------------

x
:   GmmThreshold object

y
:   ignored.

...
:   arguments to pass to [plotGmmThreshold](plotGmmThreshold.md).




Slots
-------------------



`x`
:   input distance vector with NA or infinite values removed.

`model`
:   first-second fit functions.

`cutoff`
:   type of threshold cut.

`func1.0`
:   mixing weight of the first curve.

`func1.1`
:   second parameter of the first curve (either mean of a Normal distribution or shape of a Gamma distribution).

`func1.2`
:   third parameter of the first curve (either standard deviation of a Normal distribution or scale of a Gamma distribution).

`func2.0`
:   mixing weight of the second curve.

`func2.1`
:   second parameter of the second curve (either mean of a Normal distribution or shape of a Gamma distribution)..

`func2.2`
:   third parameter of the second curve (either standard deviation of a Normal distribution or scale of a Gamma distribution).

`loglk`
:   the fit log-likelihood.

`threshold`
:   the threshold cut.

`sensitivity`
:   the sensitivity.

`specificity`
:   the specificity.

`pvalue`
:   the p-value from Hartigans' dip statistic (HDS) test, with values less than 0.05 indicating significant bimodality.




See also
-------------------

[findThreshold](findThreshold.md)



