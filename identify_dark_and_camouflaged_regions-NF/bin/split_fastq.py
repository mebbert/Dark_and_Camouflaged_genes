#!/usr/bin/env python3

import sys
import gzip
from Bio import SeqIO

# Code obtained from: https://biopython.org/wiki/Split_large_file


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


def main(sample_name, in_fq_gz_file, n_reads_per_batch):

    with gzip.open(in_fq_gz_file, "rt") as in_handle:
        record_iter = SeqIO.parse(in_handle, "fastq")

        for i, batch in enumerate(batch_iterator(record_iter, int(n_reads_per_batch))):
            out_filename = sample_name + ".interleaved_R1_R2.split_%05d.fastq.gz" % i
            with gzip.open(out_filename, "wt") as out_handle:
                count = SeqIO.write(batch, out_handle, "fastq")
            print("Wrote %i records to %s" % (count, out_filename))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
