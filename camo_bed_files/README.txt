# Description of Camo .bed Files

To run our camo variant rescue pipeline we provide three bed files for each genome build.
Following the description below, these three files can be used to call camo variants in 
CDS camo regions of the genome

All 3 .bed files contain the same column descriptions. 

The first 3 columns specify genomic coordinates

The fourth column contains a semi-colon delimited list of region_IDs. The first ID in 
the list is the region specified by that line's coordinates, all of the other IDs are 
regions that are in the same camo group and are all camouflaged by each other

The fifth column shows repeat number for that camo group (i.e., the number of times this 
camouflaged region is repeated and the number regions listed in column four)


## align_to.bed

Should be used in masking the genome to create a camo mask genome. Shows the one region in each
camo group where all reads for that group should align. This bed is expanded by 50 bp, complemented,
and then passed into bedtools maskFasta to mask the reference

## extraction.bed (previously the 'realign.bed').

Shows which regions should be extracted for realignment. Should be passed as
input into 'samtools view' to extract low quality MAPQ reads from camo regions

## GATK.bed

This bed specifies the interval GATK will use to call variants once the reads are realigned to 
the camo-masked genome. It is similar to the align_to.bed. The only difference is the GATK.bed only
has CDS regions listed and is restricted to the actual camo boundaries and not the boundaries of the 
gene body element that contain the camo region. 
