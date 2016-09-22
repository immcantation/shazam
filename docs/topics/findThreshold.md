





**findThreshold** - *Find distance threshold*

Description
--------------------

Infer value of the minimum between the two modes in a bimodal distribution.


Usage
--------------------
```
findThreshold(distances, subsample = NULL)
```

Arguments
-------------------

distances
:   numeric vector of distances.

subsample
:   number of distances to subsample for speeding up bandwidth inference.



Value
-------------------

Returns distance threshold that separates two modes of the input distribution.

Details
-------------------

The distance to nearest neighbor can be used to estimate a threshold for assigning Ig
sequences to clonal groups. A histogram of the resulting vector is often bimodal, 
with the ideal threshold being a value that separates the two modes. This function takes 
as input a vector of such distances and infers the ideal threshold.



Examples
-------------------

```R
# Subset example data to one sample as a demo
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, SAMPLE == "-1h")

# Use genotyped V assignments, HS1F model, and normalize by junction length
dist_hs1f <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", 
model="hs1f", first=FALSE, normalize="length")
threshold <- findThreshold(dist_hs1f$DIST_NEAREST)

# Plot histogram of non-NA distances
p1 <- ggplot(data=subset(dist_hs1f, !is.na(DIST_NEAREST))) + theme_bw() + 
ggtitle("Distance to nearest: hs1f") + xlab("distance") +
geom_histogram(aes(x=DIST_NEAREST), binwidth=0.025, 
fill="steelblue", color="white") + 
geom_vline(xintercept=threshold, linetype="dashed")
plot(p1)
```

![2](findThreshold-2.png)


See also
-------------------

See [distToNearest](distToNearest.md) for details on generating the input distance vector.



