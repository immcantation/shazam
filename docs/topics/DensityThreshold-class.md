





**DensityThreshold-class** - *Output of the `dens` method of findThreshold*

Description
--------------------

`DensityThreshold` contains output from the `dens` method [findThreshold](findThreshold.md).


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
:   DensityThreshold object

y
:   ignored.

...
:   arguments to pass to [plotDensityThreshold](plotDensityThreshold.md).




Slots
-------------------



`x`
:   input distance vector with NA or infinite values removed.

`bandwidth`
:   bandwidth value fit during density estimation.

`xdens`
:   x-axis (distance value) vector for smoothed density estimate.

`ydens`
:   y-axis (density) vector for smoothed density estimate.

`threshold`
:   distance threshold that separates two modes of the input distribution.




See also
-------------------

[findThreshold](findThreshold.md)



