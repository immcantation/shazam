





**findThreshold** - *Find distance threshold*

Description
--------------------

`findThreshold` automtically determines and optimal threshold for clonal assignment of
Ig sequences using a vector of nearest neighbor distances. It provides two alternative methods 
using either a Guassian Mixture Model fit (`method="gmm"`) or kernel density 
fit (`method="density"`).


Usage
--------------------
```
findThreshold(data, method = c("gmm", "density"), cutEdge = 0.9,
cross = NULL, subsample = NULL)
```

Arguments
-------------------

data
:   numeric vector containing nearest neighbor distances.

method
:   string defining the method to use for determining the optimal threshold.
One of `"gmm"` or `"density"`. See Details for methodological
descriptions.

cutEdge
:   upper range (a fraction of the data density) to rule initialization of 
Gaussian fit parameters. Default value is equal to 90
Applies only to the `"gmm"` method.

cross
:   supplementary nearest neighbor distance vector output from [distToNearest](distToNearest.md) 
for initialization of the Gaussian fit parameters. 
Applies only to the `"gmm"` method.

subsample
:   number of distances to subsample for speeding up bandwidth inference.
Applies only to the `"density"` method. If `NULL` no subsampling
is performed. As bandwith inferrence is computationally expensive, subsampling
is recommended for large data sets.




Value
-------------------


+  `"gmm"` method:      Returns a [GmmThreshold](GmmThreshold-class.md) object including the optimum 
`threshold` and the Gaussian fit parameters.
+  `"density"` method:  Returns a [DensityThreshold](DensityThreshold-class.md) object including the optimum 
`threshold` and the density fit parameters.



Details
-------------------


+  `"gmm"`:     Performs a Gaussian Mixture Model (GMM) procedure, 
including the Expectation Maximization (EM) algorithm, for learning 
the parameters  of two univariate Gaussians which fit the bimodal 
distribution entries. Retrieving the fit parameters, it then calculates
the optimum threshold, where the average of the sensitivity plus 
specificity reaches its maximum.
+  `"density"`: Fits a binned approximation to the ordinary kernel density estimate
to the nearest neighbor distances after determining the optimal
bandwidth for the density estimate via least-squares cross-validation of 
the 4th derivative of the kernel density estimator. The optimal threshold
is set as the minimum value in the valley in the density estimate
between the two modes of the distribution.




Examples
-------------------

```R
# Subset example data to one sample as a demo
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, SAMPLE == "-1h")

# Use nucleotide Hamming distance and normalize by junction length
db <- distToNearest(db, model="ham", normalize="len", nproc=1)

```

**Error in eval(expr, envir, enclos)**: could not find function "distToNearest"
```R

# Find threshold using the "gmm" method
output <- findThreshold(db$DIST_NEAREST, method="gmm")

```

**Error in eval(expr, envir, enclos)**: could not find function "findThreshold"
```R
print(output)

```

**Error in print(output)**: object 'output' not found
```R
# Plot "gmm" method results
plot(output, binwidth=0.02)

```

**Error in plot(output, binwidth = 0.02)**: object 'output' not found
```R

# Find threshold using the "density" method 
output <- findThreshold(db$DIST_NEAREST, method="density")

```

**Error in eval(expr, envir, enclos)**: could not find function "findThreshold"
```R
print(output)

```

**Error in print(output)**: object 'output' not found
```R
# Plot "density" method results
plot(output)
```

**Error in plot(output)**: object 'output' not found

See also
-------------------

See [distToNearest](distToNearest.md) for generating the nearest neighbor distance vectors.
See [plotGmmThreshold](plotGmmThreshold.md) and [plotDensityThreshold](plotDensityThreshold.md) for plotting output.



