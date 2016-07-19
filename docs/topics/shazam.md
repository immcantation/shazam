





**shazam** - *The shazam package*

Description
--------------------

Provides tools for advanced anaylisis of immunoglobulin (Ig) somatic hypermutation 
(SHM), including BASELINe, a novel method for quantifying antigen-driven selection in 
high-throughput Ig sequencing data.



Details
-------------------

Dramatic improvements in high-throughput sequencing technologies now enable 
large-scale characterization of Ig repertoires, defined as the collection of transmembrane 
antigen-receptor proteins located on the surface of T and B lymphocytes.

The `shazam` package provides tools for advanced analysis of Ig sequences following 
germline segment assignment. Namely, the analysis of SHM. 
Which includes:
 
+ Statistical analysis of SHM patterns 
Computational models and analyses of SHM have separated the process 
into two independent components: 

1. A mutability model that defines where mutations occur.
1. A nucleotide substitution model that defines the resulting mutation.

Collectively these are what form the targeting model of SHM. `shazam` 
provides tools to build these mutability and substitution (i.e. targeting) 
models.

+ BASELINe 
Bayesian Estimation of Antigen-driven Selection in Ig Sequences is a 
novel method for quantifying antigen-driven selection in high-throughput
Ig sequence data. The targeting model created using `shazam` is used 
to estimate the null distribution of expected mutation frequencies in 
BASELINe.

+ Distance calculations 
Based on the underlying SHM targeting (calculated using `shazam`) one 
can compute evolutionary distances between sequences or groups of 
sequences. This information is particularly useful in understanding and 
defining clonal relationships.

+ SHMulate 
`shazam` also provides tools for simulating immunoglobulin (Ig) somatic hypermutation.
 

Below are the functions in `shazam` broken down by the main tasks described
above:

Targeting models
-------------------



+ [createTargetingModel](createTargetingModel.md):     Build a 5-mer targeting model.
+ [plotMutability](plotMutability.md):           Plot 5-mer mutability rates.


Mutational profiling
-------------------



+ [collapseClones](collapseClones.md):           Build clonal consensus sequence.
+ [observedMutations](observedMutations.md):        Compute observed mutation counts and frequencies.
+ [expectedMutations](expectedMutations.md):        Compute expected mutation frequencies.
+ [shmulateSeq](shmulateSeq.md):              Simulate mutations in a single sequence.
+ [shmulateTree](shmulateTree.md):             Simulate mutations over a lineage tree.


Selection analysis
-------------------



+ [calcBaseline](calcBaseline.md):             Calculate the BASELINe probability
density functions (PDFs).
+ [groupBaseline](groupBaseline.md):            Combine PDFs from sequences grouped
by biological or experimental relevance.
+ [summarizeBaseline](summarizeBaseline.md):        Compute summary statistics from BASELINe PDFs.
+ [testBaseline](testBaseline.md):             Perform significance testing for the difference
between BASELINe PDFs.
+ [plotBaselineDensity](plotBaselineDensity.md):      Plot the probability density functions
resulting from selection analysis.
+ [plotBaselineSummary](plotBaselineSummary.md):      Plot summary stastistics resulting from 
selection analysis.


Distance profiling
-------------------



+ [distToNearest](distToNearest.md):            Tune clonal assignment thresholds by calculating 
distances to nearest-neighbors.
+ [calcTargetingDistance](calcTargetingDistance.md):    Construct a nucleotide distance matrix from a 
5-mer targeting model.


Simulation
-------------------



+  [shmulateSeq](shmulateSeq.md):                Simulate mutations in a single sequence.
 +  [shmulateTree](shmulateTree.md):               Simulate sequences to populate a tree.


References
-------------------


1. Hershberg U, et al. Improved methods for detecting selection by mutation 
analysis of Ig V region sequences. 
Int Immunol. 2008 20(5):683-94.
1. Uduman M, et al. Detecting selection in immunoglobulin sequences. 
Nucleic Acids Res. 2011 39(Web Server issue):W499-504.
1. Yaari G, et al. Quantifying selection in high-throughput immunoglobulin 
sequencing data sets. 
Nucleic Acids Res. 2012 40(17):e134.
1. Yaari G, et al. Models of somatic hypermutation targeting and substitution based 
on synonymous mutations from high-throughput immunoglobulin sequencing data. 
Front Immunol. 2013 4:358.
 





