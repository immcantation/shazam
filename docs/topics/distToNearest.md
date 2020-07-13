**distToNearest** - *Distance to nearest neighbor*

Description
--------------------

Get non-zero distance of every heavy chain (`IGH`) sequence (as defined by 
`sequenceColumn`) to its nearest sequence in a partition of heavy chains sharing the same 
V gene, J gene, and junction length (VJL), or in a partition of single cells with heavy chains
sharing the same heavy chain VJL combination, or of single cells with heavy and light chains 
sharing the same heavy chain VJL and light chain VJL combinations.


Usage
--------------------
```
distToNearest(
db,
sequenceColumn = "junction",
vCallColumn = "v_call",
jCallColumn = "j_call",
model = c("ham", "aa", "hh_s1f", "hh_s5f", "mk_rs1nf", "mk_rs5nf", "m1n_compat",
"hs1f_compat"),
normalize = c("len", "none"),
symmetry = c("avg", "min"),
first = TRUE,
VJthenLen = TRUE,
nproc = 1,
fields = NULL,
cross = NULL,
mst = FALSE,
subsample = NULL,
progress = FALSE,
cellIdColumn = NULL,
locusColumn = NULL,
onlyHeavy = TRUE,
keepVJLgroup = TRUE
)
```

Arguments
-------------------

db
:   data.frame containing sequence data.

sequenceColumn
:   name of the column containing the junction for grouping and for calculating
nearest neighbot distances. Note that while both heavy and light chain junctions
may be used for VJL grouping, only the heavy chain junction is used to calculate 
distances.

vCallColumn
:   name of the column containing the V-segment allele calls.

jCallColumn
:   name of the column containing the J-segment allele calls.

model
:   underlying SHM model, which must be one of 
`c("ham", "aa", "hh_s1f", "hh_s5f", "mk_rs1nf", "hs1f_compat", "m1n_compat")`.
See Details for further information.

normalize
:   method of normalization. The default is `"len"`, which 
divides the distance by the length of the sequence group. If 
`"none"` then no normalization if performed.

symmetry
:   if model is hs5f, distance between seq1 and seq2 is either the
average (avg) of seq1->seq2 and seq2->seq1 or the minimum (min).

first
:   if `TRUE` only the first call of the gene assignments 
is used. if `FALSE` the union of ambiguous gene 
assignments is used to group all sequences with any 
overlapping gene calls.

VJthenLen
:   a Boolean value specifying whether to perform partitioning as a 2-stage
process. If `TRUE`, partitions are made first based on V and J
annotations, and then further split based on junction lengths corresponding 
to `sequenceColumn`. If `FALSE`, perform partition as a 1-stage 
process during which V annotation, J annotation, and junction length are used 
to create partitions simultaneously. Defaults to `TRUE`.

nproc
:   number of cores to distribute the function over.

fields
:   additional fields to use for grouping.

cross
:   character vector of column names to use for grouping to calculate 
distances across groups. Meaning the columns that define self versus others.

mst
:   if `TRUE`, return comma-separated branch lengths from minimum 
spanning tree.

subsample
:   number of sequences to subsample for speeding up pairwise-distance-matrix calculation. 
Subsampling is performed without replacement in each VJL group of heavy chain sequences. 
If `subsample` is larger than the unique number of heavy chain sequences in each 
VJL group, then the subsampling process is ignored for that group. For each heavy chain
sequence in `db`, the reported `dist_nearest` is the distance to the closest
heavy chain sequence in the subsampled set for the VJL group. If `NULL` no 
subsampling is performed.

progress
:   if `TRUE` print a progress bar.

cellIdColumn
:   name of the character column containing cell IDs. Only applicable and required for 
single-cell mode.

locusColumn
:   name of the character column containing locus information. Only applicable and 
required for single-cell mode.

onlyHeavy
:   use only `IGH` (for BCR data) or `TRB/TRD` (for TCR data) 
sequences for grouping. Only applicable and required for single-cell mode. 
Default is `TRUE`. Also see [groupGenes](http://www.rdocumentation.org/packages/alakazam/topics/groupGenes).

keepVJLgroup
:   a Boolean value specifying whether to keep in the output the the column 
column indicating grouping based on VJL combinations. Only applicable for
1-stage partitioning (i.e. `VJthenLen=FALSE`). Also see 
[groupGenes](http://www.rdocumentation.org/packages/alakazam/topics/groupGenes).




Value
-------------------

Returns a modified `db` data.frame with nearest neighbor distances between heavy chain
sequences in the `dist_nearest` column if `cross=NULL`. If `cross` was 
specified, distances will be added as the `cross_dist_nearest` column. 

Note that distances between light chain sequences are not calculated, even if light chains 
were used for VJL grouping via `onlyHeavy=FALSE`. Light chain sequences, if any,
will have `NA` in the `dist_nearest` field.

Note that the output `vCallColumn` and `jCallColumn` columns will be converted to 
`character` if they were `factor` in the input `db`.


Details
-------------------

To invoke single-cell mode, both `cellIdColumn` and `locusColumn` must be supplied. 
Otherwise, the function will run under non-single-cell mode.

Under single-cell mode, only heavy chain sequences will be used for calculating nearest neighbor
distances. Under non-single-cell mode, all input sequences will be used for calculating nearest
neighbor distances, regardless of the values in the `locusColumn` field (if present).

Values in the `locusColumn` column must be one of `c("IGH", "IGI", "IGK", "IGL"` for BCR 
data or `"TRA", "TRB", "TRD", "TRG")` for TCR data. Otherwise, the function returns an 
error message and stops.

For single-cell mode, the input format is the same as that for [groupGenes](http://www.rdocumentation.org/packages/alakazam/topics/groupGenes). 
Namely, each row represents a sequence/chain. Sequences/chains from the same cell are linked
by a cell ID in the `cellIdColumn` field. Under this mode, there is a choice of whether 
grouping should be done using `IGH` for BCR data or `TRB/TRD` for TCR data 
sequences only, or using both `IGH and IGK/IGL` for BCR data or 
`TRB/TRD and TRA/TRG` for TCR data sequences. This is governed by `onlyHeavy`.

Note that for `distToNearest`, a cell with multiple `IGH` (for BCR data) or 
multiple `TRB/TRD` (for TCR data) is not allowed.

The distance to nearest (heavy chain) neighbor can be used to estimate a threshold for assigning 
Ig sequences to clonal groups. A histogram of the resulting vector is often bimodal, with the 
ideal threshold being a value that separates the two modes.

The following distance measures are accepted by the `model` parameter.


+  `"ham"`:          Single nucleotide Hamming distance matrix from [getDNAMatrix](http://www.rdocumentation.org/packages/alakazam/topics/getDNAMatrix) 
with gaps assigned zero distance.
+  `"aa"`:           Single amino acid Hamming distance matrix from [getAAMatrix](http://www.rdocumentation.org/packages/alakazam/topics/getAAMatrix).
+  `"hh_s1f"`:       Human single nucleotide distance matrix derived from [HH_S1F](HH_S1F.md) with 
[calcTargetingDistance](calcTargetingDistance.md).
+  `"hh_s5f"`:       Human 5-mer nucleotide context distance matix derived from [HH_S5F](HH_S5F.md) with 
[calcTargetingDistance](calcTargetingDistance.md).
+  `"mk_rs1nf"`:     Mouse single nucleotide distance matrix derived from [MK_RS1NF](MK_RS1NF.md) with 
[calcTargetingDistance](calcTargetingDistance.md).
+  `"mk_rs5nf"`:     Mouse 5-mer nucleotide context distance matrix derived from [MK_RS1NF](MK_RS1NF.md) with 
[calcTargetingDistance](calcTargetingDistance.md).
+  `"hs1f_compat"`:  Backwards compatible human single nucleotide distance matrix used in 
SHazaM v0.1.4 and Change-O v0.3.3.
+  `"m1n_compat"`:   Backwards compatibley mouse single nucleotide distance matrix used in 
SHazaM v0.1.4 and Change-O v0.3.3.


Note on `NA`s: if, for a given combination of V gene, J gene, and sequence length,
there is only 1  heavy chain sequence (as defined by `sequenceColumn`), `NA` is 
returned instead of a distance (since it has no heavy chain neighbor). If for a given combination 
there are multiple heavy chain sequences but only 1 unique one, (in which case every heavy cahin 
sequence in this group is the de facto nearest neighbor to each other, thus giving rise to distances 
of 0), `NA`s are returned instead of zero-distances.

Note on `subsample`: Subsampling is performed independently in each VJL group for heavy chain
sequences. If `subsample` is larger than number of heavy chain sequences in the group, it is 
ignored. In other words, subsampling is performed only on groups in which the number of heavy chain 
sequences is equal to or greater than `subsample`. `dist_nearest` has values calculated 
using all heavy chain sequences in the group for groups with fewer than `subsample` heavy chain
sequences, and values calculated using a subset of heavy chain sequences for the larger groups. 
To select a value of `subsample`, it can be useful to explore the group sizes in `db` 
(and the number of heavy chain sequences in those groups).


References
-------------------


1. Smith DS, et al. Di- and trinucleotide target preferences of somatic 
mutagenesis in normal and autoreactive B cells. 
J Immunol. 1996 156:2642-52. 
1. Glanville J, Kuo TC, von Budingen H-C, et al. 
Naive antibody gene-segment frequencies are heritable and unaltered by 
chronic lymphocyte ablation. 
Proc Natl Acad Sci USA. 2011 108(50):20066-71.
1. Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
on synonymous mutations from high-throughput immunoglobulin sequencing data. 
Front Immunol. 2013 4:358.
 



Examples
-------------------

```R
# Subset example data to one sample as a demo
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, sample_id == "-1h")

# Use genotyped V assignments, Hamming distance, and normalize by junction length
# First partition based on V and J assignments, then by junction length
# Take into consideration ambiguous V and J annotations
dist <- distToNearest(db, sequenceColumn="junction", 
vCallColumn="v_call_genotyped", jCallColumn="j_call",
model="ham", first=FALSE, VJthenLen=TRUE, normalize="len")

# Plot histogram of non-NA distances
p1 <- ggplot(data=subset(dist, !is.na(dist_nearest))) + 
theme_bw() + 
ggtitle("Distance to nearest: Hamming") + 
xlab("distance") +
geom_histogram(aes(x=dist_nearest), binwidth=0.025, 
fill="steelblue", color="white")
plot(p1)
```

![2](distToNearest-2.png)


See also
-------------------

See [calcTargetingDistance](calcTargetingDistance.md) for generating nucleotide distance matrices 
from a [TargetingModel](TargetingModel-class.md) object. See [HH_S5F](HH_S5F.md), [HH_S1F](HH_S1F.md), 
[MK_RS1NF](MK_RS1NF.md), [getDNAMatrix](http://www.rdocumentation.org/packages/alakazam/topics/getDNAMatrix), and [getAAMatrix](http://www.rdocumentation.org/packages/alakazam/topics/getAAMatrix)
for individual model details.






