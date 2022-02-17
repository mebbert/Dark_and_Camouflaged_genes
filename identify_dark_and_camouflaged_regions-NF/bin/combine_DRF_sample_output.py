#!/usr/bin/env python3
import os
import sys
import shutil
import subprocess
import gzip
# import numpy

def get_contig_order(hg_ref):
    hg_ref_faidx = open(hg_ref + ".fai", "r")

    contig_order = {}
    for i, line in enumerate(hg_ref_faidx):
        toks = line.strip().split('\t')

        contig = toks[0]
        contig_order[contig] = i

    return contig_order


def concatenate_DRF_files(file_genomic_regions, combined_out_file):

    # Loop over all of the files in order and concatenate

    # Create empty file (or delete contents, if exists) and write header
    f = open(combined_out_file, 'w')
    f.write("chrom\tstart\tend\tnMapQBelowThreshold\tdepth\tpercMapQBelowThreshold\n")
    f.close()
    for i, row in enumerate(file_genomic_regions):
        in_filename = row[0]

        print("Concatenating " + in_filename + " to " + combined_out_file)

        # Throw out the header for all but the first file
        if i != 0:
            # command = ' '.join(["gunzip -c", in_filename, "| tail -n +2 - | gzip - >>", combined_out_file])
            command = ' '.join(["cat", in_filename, ">>", combined_out_file])
        else:
            command = ' '.join(['cp', in_filename, combined_out_file])

        print("Command: " + command)
        subprocess.call(command, shell=True)


def main(hg_ref, sample_low_mapq_bed_dir, combined_out_file):

    contig_order = get_contig_order(hg_ref)

    # Loop over all files in the sample dir and read
    # the second line (first is a header) to determine
    # which genomic region it contains
    file_genomic_regions = []
    for filename in os.listdir(sample_low_mapq_bed_dir):
        abs_path = os.path.join(sample_low_mapq_bed_dir, filename)

        # ignore non-gz files
        if not filename.endswith(".gz"):
            continue

        print("Opening file as gzip: " + abs_path)
        f = gzip.open(abs_path, "rt")

        # read two lines
        line = f.readline()
        line = f.readline()

        # Some files will be blank (except for header) from incomplete genomic regions
        # that have 'N' or 'n' nucleotides in the reference (e.g., centromeres,
        # telomeres, etc.). Skip these files.
        if line == "":
            continue

        toks = line.strip().split('\t')
        contig = toks[0]
        start = int(toks[1])

        # Store:
        #  1. file name
        #  2. order the contig/region is listed in the .fai file
        #  3. start position of the first entry
        region = [abs_path, contig_order.get(contig), start]
        file_genomic_regions.append(region)

    # Sort file_genomic_regions by the contig order in the .fai and then the
    # start of the first entry. This will determine the order to
    # concatenate the files.
    file_genomic_regions.sort(key = lambda row: (row[1], row[2]))

#    for region in file_genomic_regions:
#        print(region)

    concatenate_DRF_files(file_genomic_regions, combined_out_file)
    

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
