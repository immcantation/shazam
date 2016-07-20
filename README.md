SHazaM
-------------------------------------------------------------------------------

SHazaM is part of the [Immcantation](http://immcantation.readthedocs.io) 
analysis framework for Adaptive Immune Receptor Repertoire sequencing 
(AIRR-seq) and provides tools for advanced analysis of somatic hypermutation 
(SHM) in immunoglobulin (Ig) sequences. Shazam focuses on the following four 
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
