# ClonalEvolution
Clonal Evolution analyses of primary and lymph node metastatic patients. 
This code assumes you already performed: 
	1) Mutaiton calling (e.g. Mutect2, VarScan, etc)
	2) Performed some sort of purity/ploidy caller (ASCAT, FACETS, etc)

Here we used Mutect2 and facets. 
Using these outputs we first remove likely FFPE artefacts using microsec. 
We filter the VCF files of the FFPE (Primary samples only)
We run vcf2maf script (Primary and lymph node)
we make Conipher annotattion tables 
Run conipher for each primary tumor sample





