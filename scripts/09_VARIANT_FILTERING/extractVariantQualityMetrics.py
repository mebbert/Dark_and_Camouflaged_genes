#!/usr/bin/python
import sys
from collections import defaultdict

def extractInfo(info):
	info_toks = info.split(';')
	info_dict = defaultdict(lambda: "NA")
	for info_tok in info_toks:
		key,value = info_tok.split('=')
		info_dict[key] = value
	return info_dict

def extractVariantMetrics(vcf):
	vcf = open(vcf, 'r')
	#A,G are purines T,C are pyrimidines
	#dictionary to be used to determine TiTv ratio
	base_type = {'A':'pur', 'G':'pur', 'T':'pyr', 'C':'pyr'}
	print("chrom\tstart\tend\ttype\tsubtype\tPASS\tset\tAF\tQD\tFS\tSOR\tMQ\tMQRankSum\tReadPosRankSum\tBaseQRankSum\tClippingRankSum\tInbreedingCoeff")
	for line in vcf:
		if line[0] == "#": continue
		toks = line.strip().split('\t')
		chrom = toks[0]
		pos = int(toks[1])
		ref = toks[3]
		alt = toks[4]
		PASS = toks[6]

		if ',' in alt: continue

		variant_type = "."
		subtype = "."
		if len(ref) == 1 and len(alt) == 1:
			variant_type = "SNP"
			if base_type[ref] == base_type[alt]: subtype = "Ti" #transition
			else: subtype = "Tv" #transversion
		elif len(ref) != len(alt):
			variant_type = "INDEL"
			if len(ref) > len(alt): subtype = "DEL"
			else: subtype = "INS"

		info = extractInfo(toks[7])
		print('\t'.join([chrom, str(pos - 1), str(pos), variant_type, subtype, PASS, info['set'], info['AF'], info['QD'], info['FS'], info['SOR'], info['MQ'],
			info['MQRankSum'], info['ReadPosRankSum'], info['BaseQRankSum'],
			info['ClippingRankSum'], info['InbreedingCoeff']))

if __name__ == "__main__":
	extractVariantMetrics(sys.argv[1])

		
		
