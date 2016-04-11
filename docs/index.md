SHazaM
-------------------------------------------------------------------------------

SHazaM provides tools for advanced analysis of somatic hypermutation (SHM) in 
immunoglobulin (Ig) sequences. Shazam focuses on the following four 
analysis topics:

1. **Statistical analysis of SHM patterns**  
   Shazam provides tools to build SHM targeting models. Models of SHM may be 
   divided into two independent components: (a) a mutability model that defines 
   where mutations occur and (b) a nucleotide substitution model that defines 
   the resulting mutation. Collectively these two components define an SHM 
   targeting model.
2. **Quantification of mutations**  
   Shazam includes methods for determine the rate of observed and expected 
   mutations under various criteria. Mutational profiling criteria include 
   rates under SHM targeting models, mutations specific to CDR and FWR 
   regions, and physicochemical property dependent substitution rates.
3. **BASELINe**  
   Bayesian Estimation of Antigen-driven Selection in Ig sequences is a 
   method for quantifying selection pressure in high-throughput Ig 
   sequencing data. Targeting models created using Shazam may be used 
   to estimate the null distribution of expected mutation frequencies in 
   BASELINe.
4. **Model-dependent distance calculations**  
   Based on the underlying SHM targeting model one can compute evolutionary 
   distances between sequences or groups of sequences. This information is 
   particularly useful in understanding and defining clonal relationships.

Related Projects
-------------------------------------------------------------------------------

* [pRESTO](http://presto.readthedocs.org) - 
  Raw read assembly, quality control and UMI processing 
* [Change-O](http://changeo.readthedocs.org) - 
  V(D)J alignment standardization and clonal clustering
* [Alakazam](http://alakazam.readthedocs.org) - 
  Lineage reconstruction, V(D)J gene usage, repertoire diversity and 
  physicochemical property analysis
* [TIgGER](http://clip.med.yale.edu/tigger) - 
  Polymorphism detection and genotyping

Dependencies
---------------

**Depends:** alakazam, ggplot2  
**Imports:** ape, data.table, doParallel, dplyr, foreach, graphics, grid, iterators, methods, lazyeval, parallel, SDMTools, scales, seqinr, stats, stringi, tidyr, utils  
**Suggests:** knitr, rmarkdown, testthat

Authors
---------------

[Mohamed Uduman](mailto:mohamed.uduman@yale.edu) (aut)  
[Gur Yaari](mailto:gur.yaari@biu.ac.il) (aut)  
[Namita Gupta](mailto:namita.gupta@yale.edu) (aut)  
[Jason Vander Heiden](mailto:jason.vanderheiden@yale.edu) (aut, cre)  
[Ang Cui](mailto:angcui@mit.edu) (ctb)  
[Susanna Marquez](mailto:susanna.marquez@yale.edu) (ctb)  
[Steven Kleinstein](mailto:steven.kleinstein@yale.edu) (aut, cph)

Citing
---------------


To cite the SHazaM package in publications, please use:

Gupta N, Vander Heiden J, Uduman M, Gadala-Maria D, Yaari G and Kleinstein S (2015). “Change-O: a toolkit for analyzing
large-scale B cell immunoglobulin repertoire sequencing data.” _Bioinformatics_, pp. 1-3. [http://doi.org/10.1093/bioinformatics/btv359](http://doi.org/10.1093/bioinformatics/btv359).

To cite the selection analysis methods, please use:

Yaari G, Uduman M and Kleinstein S (2012). “Quantifying selection in high-throughput Immunoglobulin sequencing data
sets.” _Nucleic acids research_, *40*(17), pp. e134. [http://doi.org/10.1093/nar/gks457](http://doi.org/10.1093/nar/gks457).

To cite the HS5F model and the targeting model generation methods, please use:

Yaari G, Vander Heiden J, Uduman M, Gadala-Maria D, Gupta N, Stern J, O'Connor K, Hafler D, Lasserson U, Vigneault F and
Kleinstein S (2013). “Models of somatic hypermutation targeting and substitution based on synonymous mutations from
high-throughput immunoglobulin sequencing data.” _Frontiers in Immunology_, *4*(358), pp. 1-11. [http://doi.org/10.3389/fimmu.2013.00358](http://doi.org/10.3389/fimmu.2013.00358).

