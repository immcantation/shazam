# [![](http://cranlogs.r-pkg.org/badges/grand-total/shazam)](https://www.r-pkg.org/pkg/shazam) [![](https://cranlogs.r-pkg.org/badges/shazam)](https://www.r-pkg.org/pkg/shazam) [![](https://img.shields.io/static/v1?label=AIRR-C%20sw-tools%20v1&message=compliant&color=008AFF&labelColor=000000&style=plastic)](https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html)


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



## Dependencies

**Depends:** ggplot2  
**Imports:** alakazam, ape, diptest, doParallel, dplyr, foreach, graphics, grid, igraph, iterators, KernSmooth, lazyeval, MASS, methods, parallel, progress, rlang, scales, seqinr, stats, stringi, tidyr, tidyselect, utils  
**Suggests:** knitr, rmarkdown, testthat


## Authors

[Mohamed Uduman](mailto:mohamed.uduman@yale.edu) (aut)  
[Namita Gupta](mailto:namita.gupta@yale.edu) (aut)  
[Susanna Marquez](mailto:susanna.marquez@yale.edu) (aut, cre)  
[Julian Zhou](mailto:julian.zhou@yale.edu) (aut)  
[Nima Nouri](mailto:nima.nouri@yale.edu) (aut)  
[Noah Yann Lee](mailto:noah.yann.lee@yale.edu) (aut)  
[Ang Cui](mailto:angcui@mit.edu) (ctb)  
[Jason Vander Heiden](mailto:jason.vanderheiden@gmail.com) (aut)  
[Gur Yaari](mailto:gur.yaari@biu.ac.il) (aut)  
[Steven Kleinstein](mailto:steven.kleinstein@yale.edu) (aut, cph)


## Citing

To cite the SHazaM package in publications, please use:

  Gupta N, Vander Heiden J, Uduman M, Gadala-Maria D, Yaari G,
  Kleinstein S (2015). “Change-O: a toolkit for analyzing large-scale B
  cell immunoglobulin repertoire sequencing data.” _Bioinformatics_,
  1-3. doi:10.1093/bioinformatics/btv359
  <https://doi.org/10.1093/bioinformatics/btv359>.

To cite the selection analysis methods, please use:

  Yaari G, Uduman M, Kleinstein S (2012). “Quantifying selection in
  high-throughput Immunoglobulin sequencing data sets.” _Nucleic acids
  research_, *40*(17), e134. doi:10.1093/nar/gks457
  <https://doi.org/10.1093/nar/gks457>.

To cite the HH_S5F model and the targeting model generation methods,
please use:

  Yaari G, Vander Heiden J, Uduman M, Gadala-Maria D, Gupta N, Stern J,
  O'Connor K, Hafler D, Lasserson U, Vigneault F, Kleinstein S (2013).
  “Models of somatic hypermutation targeting and substitution based on
  synonymous mutations from high-throughput immunoglobulin sequencing
  data.” _Frontiers in Immunology_, *4*(358), 1-11.
  doi:10.3389/fimmu.2013.00358
  <https://doi.org/10.3389/fimmu.2013.00358>.

To cite the HKL_S1F, HKL_S5F, MK_RS1NF, and MK_RS5NF models, please
use:

  Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T,
  O'Connor K, Vigneault F, Shlomchik M, Kleinstein S (2016). “A Model
  of Somatic Hypermutation Targeting in Mice Based on High-Throughput
  Ig Sequencing Data.” _The Journal of Immunology_, *197*(9),
  3566-3574. doi:10.4049/jimmunol.1502263
  <https://doi.org/10.4049/jimmunol.1502263>.

To cite the mixture model for threshold determination, please use:

  Nouri N, Kleinstein S (2018). “Optimized Threshold Inference for
  Partitioning of Clones From High-Throughput B Cell Repertoire
  Sequencing Data.” _Frontiers in Immunology_, *9*.
  doi:10.3389/fimmu.2018.01687
  <https://doi.org/10.3389/fimmu.2018.01687>.



## License

AGPL-3
