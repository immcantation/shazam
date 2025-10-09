[![CRAN status](https://www.r-pkg.org/badges/version/shazam)](https://cran.r-project.org/package=shazam)
[![](http://cranlogs.r-pkg.org/badges/grand-total/shazam)](https://www.r-pkg.org/pkg/shazam)
[![](https://cranlogs.r-pkg.org/badges/shazam)](https://www.r-pkg.org/pkg/shazam)
[![](https://img.shields.io/static/v1?label=AIRR-C%20sw-tools%20v1&message=compliant&color=008AFF&labelColor=000000&style=plastic)](https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html)

SHazaM
-------------------------------------------------------------------------------

SHazaM is part of the [Immcantation](http://immcantation.readthedocs.io) 
analysis framework for Adaptive Immune Receptor Repertoire sequencing 
(AIRR-seq) and provides tools for advanced analysis of somatic hypermutation 
(SHM) in immunoglobulin (Ig) sequences. Shazam focuses on the following  
analysis topics:

1. **Quantification of mutational load**  
   SHazaM includes methods for determine the rate of observed and expected 
   mutations under various criteria. Mutational profiling criteria include 
   rates under SHM targeting models, mutations specific to CDR and FWR 
   regions, and physicochemical property dependent substitution rates.
2. **Statistical models of SHM targeting patterns**  
   Models of SHM may be divided into two independent components: 
   (a) a mutability model that defines where mutations occur and (b) a 
   nucleotide substitution model that defines the resulting mutation. 
   Collectively these two components define an SHM targeting model.
   SHazaM provides empirically derived SHM 5-mer context mutation models 
   for both humans and mice, as well tools to build SHM targeting models
   from data. 
3. **Analysis of selection pressure using BASELINe**  
   The Bayesian Estimation of Antigen-driven Selection in Ig Sequences 
   (BASELINe) method is a novel method for quantifying antigen-driven 
   selection in high-throughput Ig sequence data. BASELINe uses SHM 
   targeting models can be used to estimate the null distribution of 
   expected mutation frequencies, and provide measures of selection 
   pressure informed by known AID targeting biases.
4. **Model-dependent distance calculations**  
   SHazaM provides methods to compute evolutionary distances between 
   sequences or set of sequences based on SHM targeting models. This 
   information is particularly useful in understanding and defining 
   clonal relationships.

Contact
-------------------------------------------------------------------------------

If you need help or have any questions, please contact the [Immcantation Group](mailto:immcantation@googlegroups.com).

If you have discovered a bug or have a feature request, you can open an issue using the [issue tracker](https://github.com/immcantation/shazam/issues).

To receive alerts about Immcantation releases, news, events, and tutorials, join the [Immcantation News](https://groups.google.com/g/immcantation-news) Google Group. [Membership settings](https://groups.google.com/g/immcantation-news/membership) can be adjusted to change the frequency of email updates.

