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

`a1`
:   mixing weight of the first curve.

`b1`
:   second parameter of the first curve. Either the mean of a Normal 
distribution or shape of a Gamma distribution.

`c1`
:   third parameter of the first curve. Either the standard deviation of a 
Normal distribution or scale of a Gamma distribution.

`a2`
:   mixing weight of the second curve.

`b2`
:   second parameter of the second curve. Either the mean of a Normal 
distribution or shape of a Gamma distribution.

`c2`
:   third parameter of the second curve. Either the standard deviation 
of a Normal distribution or scale of a Gamma distribution.

`loglk`
:   log-likelihood of the fit.

`threshold`
:   threshold.

`sensitivity`
:   sensitivity.

`specificity`
:   specificity.

`pvalue`
:   p-value from Hartigans' dip statistic (HDS) test. 
Values less than 0.05 indicate significant bimodality.




See also
-------------------

[findThreshold](findThreshold.md)






