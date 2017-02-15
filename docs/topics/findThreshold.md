





**findThreshold** - *Find distance threshold*

Description
--------------------

Switch between `"gmm"` or `"dens"` method to infer value of the threshold between the two modes in a bimodal distribution.


Usage
--------------------
```
findThreshold(data, method = c("gmm", "dens"), cross = NULL,
cutEdge = 0.9, subsample = NULL)
```

Arguments
-------------------

data
:   input data (a numeric vector) containing distance distribution.

method
:   either one of the `"gmm"` or `"dens"` techniques.

cross
:   a supplementary info (numeric vector) invoked from [distToNearest](distToNearest.md) 
function, to support initialization of the Gaussian fit parameters.

cutEdge
:   upper range (a fraction of the data density) to rule initialization of 
Gaussian fit parameters. Default value is equal to <code class = 'eq'>90</code>% of the entries.

subsample
:   number of distances to subsample for speeding up bandwidth inference.




Value
-------------------


+  Calling `"gmm"` method, it returns an object including optimum "`threshold`" 
cut and the Gaussian fit parameters, such as mixing proportion ("`omega1`" and "`omega2`"), 
mean ("`mu1`" and "`mu2`"), and standard deviation ("`sigma1`" and "`sigma2`").
Returns "`NULL`" if no fit has found. See also [GmmResults](GmmResults-class.md) class.
+  Calling `"dens"` method, it returns distance threshold that separates two modes 
of the input distribution. Returns "`NULL`" if no cut has found. 
See also [DensResults](DensResults-class.md) class.



Details
-------------------


+  `"gmm"`: This function follows a Gaussian Mixture Model (GMM) procedure, 
including the Expectation Maximization (EM) algorithm, for learning the parameters  
of two univariate Gaussians which fit the bimodal distribution entries. 
Retrieving the fit parameters, it then calculates, analytically, the optimum threshold, 
where the average of the Sensitivity plus Specificity reaches its maximum. This threshold 
can be then invoked for assigning Ig sequences to clonal groups.

+  `"dens"`: The distance to nearest neighbor can be used to estimate a threshold for assigning Ig
sequences to clonal groups. A histogram of the resulting vector is often bimodal, 
with the ideal threshold being a value that separates the two modes. This function takes 
as input a vector of such distances and infers the ideal threshold.




Examples
-------------------

```R
# Subset example data to one sample as a demo
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, SAMPLE == "-1h")

# Use nucleotide Hamming distance and normalize by junction length
db <- distToNearest(db, model="ham", first=FALSE, normalize="length", nproc=1)

# To find the Threshold cut use findThreshold-switch for "gmm" method 
output <- findThreshold(db$DIST_NEAREST, method="gmm", cutEdge=0.9)

```


```
[1] "The number of non-NA entries= 958"
[1] "The 'gmm' would be done in 51 iterations"
##################################################

```


```R

# Retrieve outputs:
omega <- c(output@omega1, output@omega2)
mu <-    c(output@mu1,    output@mu2) 
sigma <- c(output@sigma1, output@sigma2) 
threshold <- output@threshold

# To check the quality of the fit performance and corresponding 
# threshold location use "plotGmmFit" function in shazam. 

# using findThreshold switch for "dens" method
output <- findThreshold(db$DIST_NEAREST, method="dens")

# Retrieve outputs:
threshold <- output@threshold

# Plot histogram of non-NA distances
p1 <- ggplot(data=subset(db, !is.na(DIST_NEAREST))) + theme_bw() + 
ggtitle("Distance to nearest: hs1f") + xlab("distance") +
geom_histogram(aes(x=DIST_NEAREST), binwidth=0.025, 
fill="steelblue", color="white") + 
geom_vline(xintercept=threshold, linetype="dashed")
plot(p1)
```

![4](findThreshold-4.png)



