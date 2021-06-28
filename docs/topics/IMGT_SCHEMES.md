**IMGT_SCHEMES** - *IMGT unique numbering schemes*

Description
--------------------

Sequence region definitions according to the IMGT unique numbering scheme.






Format
-------------------

A [RegionDefinition](RegionDefinition-class.md) object defining:

+ `IMGT_V`:              The IMGT numbered V segment up to position nucleotide 312.
This definition combines the CDR1 and CDR2 into a single CDR region,
and FWR1, FWR2 and FWR3 into a single FWR region. CDR3 and FWR4 are
excluded as they are downstream of nucleotide 312.
+ `IMGT_V_BY_CODONS`:    The IMGT numbered V segment up to position nucleotide 312.
This definition treats each codon, from codon 1 to codon 104, as a 
distinct region.
+ `IMGT_V_BY_REGIONS`:   The IMGT numbered V segment up to position nucleotide 312.
This defines separate regions for each of CDR1, CDR2,
FWR1, FWR2 and FWR3. CDR3 and FWR4 are
excluded as they are downstream of nucleotide 312.
+ `IMGT_V_BY_SEGMENTS`:  The IMGT numbered V segment up to position nucleotide 312.
This definition has no subdivisons and treats the entire V segment
as a single region.
+ `IMGT_VDJ`:            The IMGT numbered segments of FWR1/2/3/4 and CDR1/2/3.
This definition combines regions of CDR1, CDR2, CDR3 into a single CDR region, 
and FWR1, FWR2, FWR3 FWR4 into a single FWR region.
Note that until function [makeRegion](makeRegion.md) will be applied
- this `IMGT_SCHEMES` will have the slot `seqLength`
with value 0, and the `boundaries` slot will be empty. This is since
these slots depend on the junction length which is unknown yet.
After [makeRegion](makeRegion.md) is applied - these slots get specific values
per the specific sequence and junction length.
+ `IMGT_VDJ_BY_REGIONS`:    The IMGT numbered segments of FWR1/2/3/4 and CDR1/2/3.
This defines separate regions for each of CDR1, CDR2, CDR3, 
FWR1, FWR2, FWR3 and FWR4. 
Note that until function [makeRegion](makeRegion.md) will be applied
- this `IMGT_SCHEMES` will have the slot `seqLength`
with value 0, and the `boundaries` slot will be empty. This is since
these slots depend on the junction length which is unknown yet.
After [makeRegion](makeRegion.md) is applied - these slots get specific values
per the specific sequence and junction length.



References
-------------------


1. Lefranc MP, et al. IMGT unique numbering for immunoglobulin and T cell 
receptor variable domains and Ig superfamily V-like domains. 
Developmental and comparative immunology. 2003 27:55-77.










