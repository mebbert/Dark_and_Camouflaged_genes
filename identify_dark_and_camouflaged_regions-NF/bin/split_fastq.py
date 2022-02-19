#!/usr/bin/env python3

import sys
import gzip
from Bio import SeqIO


def main(sample_name, in_fq_gz_file, n_reads_per_batch):

    n_reads_per_batch = int(n_reads_per_batch)
    with gzip.open(in_fq_gz_file, "rt") as in_handle:
        record_iter = SeqIO.parse(in_handle, "fastq")

        i = 0
        while True:

            out_filename = sample_name + ".interleaved_R1_R2.split_%05d.fastq.gz" % i
            i += 1
            with gzip.open(out_filename, "wt") as out_handle:

                n_reads = 0
                while n_reads < n_reads_per_batch:
                    n_reads += 1
                    try:
                        read = next(record_iter)
                    except StopIteration:
                        read = None
                    if read is None:
                        # End of file
                        return
                    SeqIO.write(read, out_handle, "fastq")
                
                print("Wrote %i reads to %s." % (n_reads, out_filename))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
