





**shazam** - *The shazam package*

Description
--------------------

Dramatic improvements in high-throughput sequencing technologies now enable 
large-scale characterization of Ig repertoires, defined as the collection of transmembrane 
antigen-receptor proteins located on the surface of T and B lymphocytes. The `shazam`
package provides tools for advanced analysis of somatic hypermutation (SHM) in
immunoglobulin (Ig) sequences. The key functions in `shazam`, broken down topic, are 
described below.



Mutational profiling
-------------------


`shazam` provides tools to quantify the extent and nature of SHM within
full length V(D)J sequences as well as sub-regions (eg, FWR and CDR).
Quantification of expected mutational loaded, under specific SHM targeting 
models, can also be performed along with model driven simulations of SHM.


+ [collapseClones](collapseClones.md):           Build clonal consensus sequence.
+ [observedMutations](observedMutations.md):        Compute observed mutation counts and frequencies.
+ [expectedMutations](expectedMutations.md):        Compute expected mutation frequencies.
+ [shmulateSeq](shmulateSeq.md):              Simulate mutations in a single sequence.
+ [shmulateTree](shmulateTree.md):             Simulate mutations over a lineage tree.


SHM targeting models
-------------------


Computational models and analyses of SHM have separated the process 
into two independent components: 

1. A mutability model that defines where mutations occur.
1. A nucleotide substitution model that defines the resulting mutation.

Collectively these are what form the targeting model of SHM. `shazam` 
provides empirically derived targeting models for both humans and mice,
along with tools to build these mutability and substitution models from data.


+ [createTargetingModel](createTargetingModel.md):     Build a 5-mer targeting model.
+ [plotMutability](plotMutability.md):           Plot 5-mer mutability rates.
+ [HS5FModel](HS5FModel.md):                Human 5-mer SHM targeting model.
+ [MRS5NFModel](MRS5NFModel.md):              Mouse 5-mer SHM targeting model.


Quantification of selection pressure
-------------------


Bayesian Estimation of Antigen-driven Selection in Ig Sequences is a 
novel method for quantifying antigen-driven selection in high-throughput
Ig sequence data. Targeting models created using `shazam` can be used 
to estimate the null distribution of expected mutation frequencies used
by BASELINe, providing measures of selection pressure informed by known 
AID targeting biases.


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


Mutational distance calculation
-------------------


`shazam` provides tools to compute evolutionary distances between 
sequences or groups of sequences, which can leverage SHM targeting 
models. This information is particularly useful in understanding and 
defining clonal relationships.


+ [findThreshold](findThreshold.md):            Identify clonal assignment threshold based on 
distances to nearest neighbors.
+ [distToNearest](distToNearest.md):            Tune clonal assignment thresholds by calculating 
distances to nearest neighbors.
+ [calcTargetingDistance](calcTargetingDistance.md):    Construct a nucleotide distance matrix from a 
5-mer targeting model.


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
1. Cui A, et al. A model of somatic hypermutation targeting in mice based on 
high-throughput immunoglobulin sequencing data. Under review.
 





