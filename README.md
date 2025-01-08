[![](http://cranlogs.r-pkg.org/badges/grand-total/shazam)](https://www.r-pkg.org/pkg/shazam)
[![](https://cranlogs.r-pkg.org/badges/shazam)](https://www.r-pkg.org/pkg/shazam)
[![](https://img.shields.io/static/v1?label=AIRR-C%20sw-tools%20v1&message=compliant&color=008AFF&labelColor=000000&style=plastic)](https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html)

!!! important "2025 Immcantation Users Group Meeting"
    *Are you an Immcantation user and/or interested in adaptive immune receptor repertoire analysis?*
    
    Register now for the upcoming Immcantation Users Group Meeting!
    It will be held virtually on **January 30th, 2025, from 10 to 1:30pm (ET)**.
    All talks will be from user-submitted abstracts.

    Full information here: [https://immcantation.github.io/users-meeting/](https://immcantation.github.io/users-meeting/)

**IMPORTANT!** 
SHazaM has moved to https://github.com/immcantation/shazam

To update Git configuration settings use:

```
   git config user.email "your-gh-user@email.com"
   git config user.name "your-gh-user-name"
   git remote set-url origin git@github.com:immcantation/shazam.git
```

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

For help and questions please contact the [Immcantation Group](mailto:immcantation@googlegroups.com)
or use the [issue tracker](https://github.com/immcantation/shazam/issues?q=is%3Aissue+is%3Aopen+).
