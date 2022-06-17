#!/usr/bin/python3
import sys

def main():
    blat_bed = open('./blat.results.bed', 'r')
    print("we are in the script")
    small_blat_bed = open('./small.blat.results.bed', 'w')
    large_blat_bed = open('./large.blat.results.bed', 'w')
    for line in blat_bed:
        if line.startswith("#"):
            continue
        toks = line.strip().split()
        if int(toks[2]) - int(toks[1]) > 45000:
            print("writing to large file")
            print(int(toks[2]), int(toks[1]), int(toks[2])- int(toks[1]))
            large_blat_bed.write(line)
        else:
            small_blat_bed.write(line)

    blat_bed.close()
    small_blat_bed.close()
    large_blat_bed.close()
                                    

if __name__ == "__main__":
    main()


