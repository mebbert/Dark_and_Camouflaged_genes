import numpy as np
import sys
import math

def extractAnnos(tok):
	annos = {}
	annos_toks = tok.split(';')
	for annos_tok in annos_toks:
		key,value = annos_tok.split('=')
		annos[key] = value
	return annos

def choose(n,r):
	f = math.factorial
	return f(n) / f(r) / f(n-r)

def getInbreeding(genotypes):
	ploidy = genotypes[0].count('/') + 1
	ogeno = [0] * (ploidy + 1) #observed genotype counts
	samples = 0
	ref_total = 0
	for geno in genotypes:
		if geno[0] == '.': continue
		samples += 1
		ref_count = geno.count('0')
		ref_total += ref_count
		ogeno[ploidy - ref_count] += 1
	if sum(ogeno[1:]) == 1 and ogeno[-1] == 1: return float('inf')
	elif sum(ogeno[1:]) < 10: return 0.0
	p = ref_total / float(ploidy * samples)
	q = 1 - p
	egeno = [0] * (ploidy + 1) #expected genotype counts from hwe
	for n in range(ploidy + 1):
		egeno[n] = samples * choose(ploidy,n) * (p**(ploidy - n)) * (q**n)

	inbreedingCoeff = round(1 - ((sum(ogeno[1:ploidy]) + 1) / float(sum(egeno[1:ploidy]) + 1)), 5)
	return inbreedingCoeff

def updateVcf(vcf):
	vcf = open(vcf, 'r')
	for line in vcf:
		if line[0] == '#': 
			sys.stdout.write(line)
			continue
		toks = line.strip().split('\t')
		#Skip Sex Chromosomes
		if toks[0] == "X" or toks[0] == "Y":
			return
		if toks[6] != "PASS" and toks[6] != ".": 
			sys.stdout.write(line)
			continue
		genos = [ x.split(':')[0] for x in toks[9:] ]
		inbreeding = getInbreeding(genos)
		if inbreeding == float('inf'):
			toks[6] = "homo_singleton"
		else: 
			toks[6] = "PASS"
			toks[7] += ";InbreedingCoeff=%f" % inbreeding
		sys.stdout.write('%s\n' % '\t'.join(toks))


if __name__ == "__main__":
	updateVcf(sys.argv[1])
