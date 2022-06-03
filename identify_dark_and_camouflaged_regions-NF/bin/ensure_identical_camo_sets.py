#!/usr/bin/env python3

import sys


def format_camo_sets(lines):
    bad_lines = []
    camo_set = set()
    for line in lines:
        camo_set.update(line[3].split(';'))

    for line in lines:
        #line[3] = list(camo_set)
        if line[3].split(';').sort() != list(camo_set).sort():
            bad_lines.append(line)

    return bad_lines


def main(mask_bed_file):
    try:
        bad_lines = []
        with open(mask_bed_file, 'r') as mask_bed:
            previous_lines = []
            for line in mask_bed:
                current_line = line.strip().split('\t')
                if previous_lines != [] and (current_line[0] != previous_lines[-1][0] or current_line[1] > previous_lines[-1][2]):
                    if len(previous_lines) > 1:
                        bad_lines = format_camo_sets(previous_lines) + bad_lines
                        previous_lines = []
                    
                else:
                    previous_lines.append(current_line)
        print(bad_lines)
    except ValueError as e:
        print(e)

if __name__ == "__main__":
    main(sys.argv[1])

