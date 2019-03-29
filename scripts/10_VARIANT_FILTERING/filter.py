import sys
from collections import defaultdict

def loadFalsePositives(fp_bed):
	fp_bed = open(fp_bed, 'r')
	fps = set()
	for line in fp_bed:
		toks = line.strip().split("\t")
		fps.add((toks[0], int(toks[1])))

	return fps

def extractAnnos(tok):
	annos = defaultdict(int)
	annos_toks = tok.split(';')
	for annos_tok in annos_toks:
		key,value = annos_tok.split('=')
		annos[key] = value
	return annos

def updateVcf(vcf, fp_bed):
	fps = loadFalsePositives(fp_bed)
	vcf = open(vcf, 'r')
	for line in vcf:
		if line[0] == '#': 
			sys.stdout.write(line)
			continue
		toks = line.strip().split('\t')
		#Skip Sex Chromosomes
		if toks[0] == "X" or toks[0] == "Y":
			return
		if toks[6] != "PASS": 
			sys.stdout.write(line)
			continue
			
		annos = extractAnnos(toks[7])
		QD = float(annos['QD'])
		inbreedingCoeff = float(annos['InbreedingCoeff'])
		if (toks[0], int(toks[1])) in fps:
			toks[6] = "reference-based artifact"
		elif QD < 2.0:
			toks[6] = "qd_filter"
		#elif inbreedingCoeff  > .5 or inbreedingCoeff < -.5:
		#			toks[6] = "inbreeding_filter"
		sys.stdout.write('%s\n' % '\t'.join(toks))


if __name__ == "__main__":
	updateVcf(sys.argv[1], sys.argv[2])
