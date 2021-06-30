**IMGT_SCHEMES** - *IMGT unique numbering schemes*

Description
--------------------

Sequence region definitions according to the IMGT unique numbering scheme.






Format
-------------------

A [RegionDefinition](RegionDefinition-class.md) object defining:

+ `IMGT_V`:               The IMGT numbered V segment up to position nucleotide 312.
This definition combines the CDR1 and CDR2 into a single CDR region,
and FWR1, FWR2 and FWR3 into a single FWR region. CDR3 and FWR4 are
excluded as they are downstream of nucleotide 312.
+ `IMGT_V_BY_CODONS`:     The IMGT numbered V segment up to position nucleotide 312.
This definition treats each codon, from codon 1 to codon 104, as a 
distinct region.
+ `IMGT_V_BY_REGIONS`:    The IMGT numbered V segment up to position nucleotide 312.
This defines separate regions for each of CDR1, CDR2,
FWR1, FWR2 and FWR3. CDR3 and FWR4 are
excluded as they are downstream of nucleotide 312.
+ `IMGT_V_BY_SEGMENTS`:   The IMGT numbered V segment up to position nucleotide 312.
This definition has no subdivisons and treats the entire V segment
as a single region.
+ `IMGT_VDJ`:             IMGT numbered regions for CDR1-3 and FWR1-4 with combined CDR and FWR 
definitions spanning CDR1-3 and FWR1-4, respectively.
Note, unless the definition object has been updated using [setRegionBoundaries](setRegionBoundaries.md) 
this schema will have a value of `0` for the `seqLength` slot and
the `boundaries` slot will be empty. This is because
these slots depend on the junction length which is unknown in the template 
scheme. After [setRegionBoundaries](setRegionBoundaries.md) has been run, these slots will be populated
with the appropriate values for the specied sequence and junction length.
+ `IMGT_VDJ_BY_REGIONS`:  The IMGT numbered regions for FWR1-4 and CDR1-3 with separate region boundaries
for each of CDR1, CDR2, CDR3, FWR1, FWR2, FWR3 and FWR4. 
Note, unless the definition object has been updated using [setRegionBoundaries](setRegionBoundaries.md) 
this schema will have a value of `0` for the `seqLength` slot and
the `boundaries` slot will be empty. This is because
these slots depend on the junction length which is unknown in the template 
scheme. After [setRegionBoundaries](setRegionBoundaries.md) has been run, these slots will be populated
with the appropriate values for the specied sequence and junction length.



References
-------------------


1. Lefranc MP, et al. IMGT unique numbering for immunoglobulin and T cell 
receptor variable domains and Ig superfamily V-like domains. 
Developmental and comparative immunology. 2003 27:55-77.










