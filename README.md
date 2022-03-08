# Systematic analysis of dark and camouflaged genes

In our original paper [Systematic analysis of dark and camouflaged genes reveals disease-relevant
genes hiding in plain sight](https://doi.org/10.1186/s13059-019-1707-2),
we presented a list of dark and camouflaged genes and .bed files specifying their genomic coordinates. 
We originally provided the exact scripts to replicate our work, but those scripts were not well adapted
to using this technique on other data sets. Here, we are providing a [Nextflow](https://www.nextflow.io/)
workflow that we hope is much more user friendly. 
these regions in short-read data. 

---

## What are dark genes?

Essentially, dark genes/regions are genomic regions that cannot be adequately assembled or aligned using
standard short-read sequncing technologies, preventing researchers from calling mutations in these regions.
We define two subclasses of dark regions (using very conservative cutoffs): `dark by depth` (where ≤ 5
reads map to the region) and `dark by MAPQ` (where reads align but ≥ 90% of reads have MAPQ < 10). A specific
subset of dark by MAPQ genes are `camouflaged genes` (or regions), where camouflaged genes/regions are the result
of duplication events in the genome. The cutoffs we used are extremely conservative, which means the problem is actually
larger than what is reported in our paper.

## Example Dark and Camouflaged genes

![Example of dark and camouflaged regions](./imgs/dark_camo_example.png)

The IGV pile-ups above show examples of these dark regions: **(a)** _HLA-DRB5_ is `dark by depth`, because
no reads align to this region; **(b)** _HSPA1A_ is `dark by MAPQ` because reads align, but the aligned reads
have a poor MAPQ (MAPQ == 0) and would be completely overlooked by variant callers; and **(c)** _HSPA1A_ is `dark
by MAPQ` because of a gene duplication event: it is camouflaged by HSPA1B. The two genes are nearly
identical so aligners can't determine from which gene a read originated. Thus, the aligner randomly
assigns the reads to one of the two locations and gives it a mapping quality of 0. **NOTE:** the
reason these reads have a MAPQ of 0 is _not_ because the reads do not align well&mdash;they
align perfectly&mdash;but they align to _multiple places_ perfectly.

Dark regions were not a new concept before our paper, but we discovered dark regions are more prevelant
across the human genome than most realized. Based on standard
whole-genome Illumina sequencing data, we identified 36794 dark regions in 6054 gene bodies (3804
protein-coding) from pathways important to human health, development, and reproduction. Of the 6054
gene bodies, 527 (8.7%) were 100% dark (117 protein-coding) and 2128 (35.2%) were ≥5% dark (592
protein-coding). We found 2855 dark regions were in protein-coding exons (CDS) across 748 genes.
Many of these genes are important to human health and disease. 

We have also developed a conceptually simple algorithm to resolve most of these camouflaged regions in short-read
data. As a proof of concept, we originally our algorithm to the Alzhiemer's Disease Sequencing Project (ADSP)
whole genome/exome case-control study and uploaded those results [NIAGADS](https://www.niagads.org/home). We discovered
an interesting loss-of-function frameshift variant in _CR1_, a top-five Alzheimer's disease, that was found only in 5 Alzhiemers
disease cases (0 controls). While we do not currently have enough statistical power to formally
assess whether the frameshift variant we discovered is driving Alzheimer's disease, it serves as an
important proof of concept that important variants are being overlooked. 

---

## Running our analysis

In our original release, we divided up our scripts into 11 distinct steps to replicate the results
from our paper. We have now broken the workflow into two separate workflows: **(1)** a workflow to identify
dark and camouflaged regions for a given genome version and to mask that genome for use in the
second workflow; and **(2)** to rescue variants from camouflaged regions for a given data set (e.g., as we
did for the ADSP).

> IMPORTANT: While we typically think of the human reference genome in major versions (e.g.,
> hg19/b37 and GRCh38), the camouflaged regions can vary dramatically depending on which contigs are
> included in the specific version used for alignment. e.g., the 1000 Genome consortium included all
> unplaced contigs and all HLA variants. In one sense, this is ideal to help reads align to their
> proper place. This approach dramatically increases the number of regions that will be camouflaged,
> and thus will be overlooked. i.e., _the more complete and accurate the reference genome becomes,
> the more regions will be camouflaged when using short-read technologies._ **The real solution to
> resolving camouflaged regions is to using long-read sequencing platorms (e.g., [Oxford Nanopore
> Technologies](https://nanoporetech.com/) and [PacBio](https://www.pacb.com/)).**

---

## Using our camouflaged .bed files

We encourage researchers to use our methods and the provided .bed files to call camo variants within
their own data sets.

We outline our short-read data pipeline to rescue camo variants below:

![Rescue Camo Variants Image](./imgs/rescue_pipeline.png)

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
will be masked. The complement of the alignto .bed is used to mask the genome. When masking the genome, each
region in the alignto .bed is expanded by 50bp so that reads extracted near the boundary of a camo
region can accurately align. 

The final column of the alignto and realign .bed lists the repeat number of the
camoulog group--the number of times this sequence appears in the genome with high homology. We
recommend only calling variants in regions with ≤ 5 repeat number. Depending on this repeat number we
also call variants using GATK with different ploidies. Variants within a region where the camoulog 
group repeat number of 2, are called with a ploidy of 4. Whereas in the CR1 example, we would use a
HaplotypeCaller ploidy of 6.

The GATK.bed provided is essentially the same as the align\_to.bed with the exception that it is 
restricted to CDS regions and only lists the interval of that actual camo region and not the whole 
gene body. The GATK.bed is passed into GATK to define the interval over which to call variants. We
only called variants strictly in camo regions and only in CDS regions to avoid calling too many 
false positives.

In order to use these .bed files to call camouflaged variants in your own data set, follow the scripts found
in the steps **06\_CREATE\_BED\_FILE** to **10\_VARIANT\_FILTERING** directories. A brief outline of the workflow is as
follows:

For each repeat number to be tested (we recommend ≤ 5), do the following:

1. Expand the alignto .bed file by 50 and use it to mask the genome for this ploidy

2. Use the realign .bed file to extract reads from the correct repeat number regions

3. Align these reads to the new camo masked reference with bwa mem

4. Use GATK HaplotypeCaller to call variants specifying the correct ploidy and setting the interval
   to the GATK .bed

5. Combine and Genotype gVCF files across the whole cohort to create the ploidy specific VCF

6. Filter out false positive variants with QD and InbreedingCoeff cutoffs and remove variants found
   in reference-based artifact positions

7. Convert polyploid VCFs to diploid and use plink2 to test associations

Any variants found using these methods should be investigated independentally and then experimentally 
validated to ensure they are not false positives and then to determine in which specific camouflaged region 
the variant actually lies. 
