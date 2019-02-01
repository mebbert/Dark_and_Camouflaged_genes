# Systematic analysis of dark and camouflaged genes
We anticipate having all code and .bed files posted here by the end of January 2019.

**IN PROGRESS**
### Update 1 February 2019
Filled in description of how to use our camouflaged bed files. Still working on writing rest of
README
### Update 29 January 2019
We just uploaded scripts and .bed files. We are working on a write up for how to run the
camouflaged gene analysis.

## What are dark genes?

## Using our camouflaged .bed files

We provide two .bed files to define and call variants in camouflaged genes. The .bed files were
created from 10 whole genomes generated with Illumina paired-end read Sequencing with a 100bp read
length. The realign .bed file defines all gene-body elements (i.e. exon, intron, UTR, etc.) that
contain a camo region and which other gene-body element(s) share high sequence identity which causes
the regions to camouflage  eachother. The fourth column of the realign .bed file contains a semi-colon delimited
list of gene-body element IDs. The first ID corresponds to the given genomic position, while all
subsequent region IDs are regions with high sequence identity that are camouflaged together. The
realign file is used as an argument to samtools view to extract all low-MAPQ reads that fall within 
camouflaged genes. The entire gene-body element is used in the realign file as opposed to just the 
camouflaged region in order to potentially capture more variants in flanking sequences (which tend 
to be highly homologous yet failed to meet our strict cutoffs).

The alignto .bed file is similar to the realign .bed file, but it only contains the coordinates for
one gene-body element per camoulog group. For example, CR1 exons 10, 21, and 26 are nearly identical
and are all camouflaged to eachother. The realign .bed will contain the coordinates for all three
exons, while the alignto .bed will only contain coordinates for exon 10. Thus, the alignto file
defines which repeated region will be used to align the extracted reads, and which repeated regions
will be masked. The complement of the alignto .bed is used to mask the genome, and then the alignto
.bed itself is used as the interval .bed in GATK to call variants. When masking the genome, each
region in the alignto .bed is expanded by 50bp so that reads extracted near the boundary of a camo
region can accurately align. 

The final column of the alignto and realign .bed lists the repeat number of the
camoulog group--the number of times this sequence appears in the genome with high homology. We
recommend only calling variants in regions with ≤ 5 repeat number. Depending on this repeat number we
also call variants using GATK with different ploidies. Variants within a region where the camoulog 
group repeat number of 2, are called with a ploidy of 4. Whereas in the CR1 example, we would use a
HaplotypeCaller ploidy of 6.

In order to use these .bed files to call camouflaged variants in your own data set, follow the scripts found
in the **06\_CREATE\_BED\_FILE** and **07\_RUN\_ADSP** directories. A brief outline of the workflow is as
follows:

For each repeat number to be tested, do the following:

1. Expand the alignto .bed file by 50 and use it to mask the genome for this ploidy

2. Use the realign .bed file to extract reads from the correct repeat number regions

3. Align these reads to the new camo masked reference with bwa mem

4. Use GATK HaplotypeCaller to call variants specifying the correct ploidy and setting the interval
   to the alignto .bed

5. Combine and Genotype gVCF files across the whole cohort to create the ploidy specific VCF

6. Filter out false positive variants with QD and InbreedingCoeff cutoffs 

7. Convert polyploid VCFs to diploid and use plink2 to test associations

Any variants found using these methods should be investigated independentally and then experimentally 
validated to ensure they are not false positives and then to determine in which specific camouflaged region 
the variant actually lies. 

## Rerunning our analysis

In order to replicate our results or to apply these same methods to other datasets, we provide all
of our scripts organized into 11 distinct steps. Each step is encompassed into a directory with all
necesary scripts to run that step. Each directory contains a submit.sh shell script to automate the
running of each step. At the top each submit.sh script is a list of entry points where the
user should fill in the correct paths for their data before running the script.  
