





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

`omega1`
:   first Gaussain mixing proportion.

`omega2`
:   second Gaussain mixing proportion.

`mu1`
:   first Gaussian mean.

`mu2`
:   second Gaussain mean.

`sigma1`
:   first Gaussain standard deviation.

`sigma2`
:   second Gaussain standard deviation.

`threshold`
:   optimum threshold cut.




See also
-------------------

[findThreshold](findThreshold.md)



