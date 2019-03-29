import json
import random
import gzip
from sets import Set
from collections import defaultdict
import io

json_vcf = gzip.open('HGMD_PRO.vcf.tsv.bgz', 'r')

dark_cds_rl100 = open('../hg38/illuminaRL100/illuminaRL100.hg38.percent_dark_genes.txt', 'r')

gene_assocs_rl100_CDS_only = io.open('CDS_dark_gene_disease_associations-illumina_RL_100.txt', 'w', encoding = 'utf-8')
pheno_out = io.open('phenotypes_gene_names-illumina_RL_100.txt', 'w', encoding = 'utf-8')

# Get protein-coding genes with >= 5% exonic regions (CDS and UTR) that are dark
dark_genes_rl100_CDS_only = {}
for line in dark_cds_rl100:
    if line.startswith('gene_name'):
        continue
    line = line.strip()
    toks = line.split('\t')
    gene = toks[0]
    biotype = toks[1]
    perc_cds_dark = float(toks[2])
    perc_exons_dark = float(toks[4])

    if biotype == 'protein coding' and perc_cds_dark >= 5:
        dark_genes_rl100_CDS_only[gene] = gene


genes_rl100_CDS_only = Set()
hits_rl100_CDS_only = {}
all_NGMD_genes = Set()
all_NGMD_phenotypes = Set()
pheno_dist = defaultdict(int)
gene_pheno_tuples = Set()
for line in json_vcf:
    line = line.strip()
    toks = line.split('\t')
    chrom = toks[0]
    start = toks[1]
    end = toks[2]
    json_data = json.loads(toks[3])

    gene = json_data['GENE']
    all_NGMD_genes.add(gene)
    pheno = ','.join(json_data['PHEN'])
    all_NGMD_phenotypes.add(pheno)
    gene_pheno_tuples.add((gene, pheno))

    key = gene + ':' + pheno
    if gene in dark_genes_rl100_CDS_only:
        genes_rl100_CDS_only.add(gene)
        if key not in hits_rl100_CDS_only:
            pheno_dist[pheno] += 1
            hits_rl100_CDS_only[key] = key
            pheno_out.write(u'%s\t%s\n' % (pheno, gene))

pheno_out.close()

print 'Found ' + str(len(all_NGMD_phenotypes)) + ' unique phenotypes in NGMD'

print 'Found ' + str(len(genes_rl100_CDS_only)) + ' genes across ' + str(len(hits_rl100_CDS_only)) + ' disease associations (RL 100; CDS only).'

gene_assocs_rl100_CDS_only.write(u'num\t%s' % '\t'.join(all_NGMD_phenotypes))
gene_assocs_rl100_CDS_only.write(u'\n1')
for pheno in all_NGMD_phenotypes:
    gene_assocs_rl100_CDS_only.write(u'\t%d' % pheno_dist[pheno])

exit()

all_NGMD_genes = list(all_NGMD_genes)
all_NGMD_phenotypes = list(all_NGMD_phenotypes)
gene_pheno_tuples = list(gene_pheno_tuples)

null_distribution = []
for i in range(10000):
    gene_set = random.sample(all_NGMD_genes, len(genes_rl100_CDS_only))
    hits = {}
    null_pheno_dist = defaultdict(int)
    for gene,pheno in gene_pheno_tuples:
        key = gene + ':' + pheno
        if gene in gene_set:
            if key not in hits:
                null_pheno_dist[pheno] += 1
                hits[key] = key
    null_distribution.append(null_pheno_dist)

null_matrix = io.open('null_matrix.txt', 'w', encoding='utf-8')
null_matrix.write(u'bootstrap_num\t%s' % '\t'.join(all_NGMD_phenotypes))
for i,dist in enumerate(null_distribution):
    null_matrix.write(u'\n%d' % i)
    for pheno in all_NGMD_phenotypes:
        null_matrix.write(u'\t%d' % dist[pheno])
null_matrix.close()

    

