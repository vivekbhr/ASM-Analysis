#!/usr/bin/env python
# author: Fidel Ramirez, Sarah Diehl
# date: 9.11.2013

# version 3: improved removal of multiple mapping reads by recommendations from Friederike:
# - early versions of bowtie2 (< 2.2.0) contained a bug that led to irregular and often faulty setting of the XS tag for paired-end data
# - the mere presence of the XS tag is not enough to indicate a multi-read, instead, one should additionally check the AS tag that contains the alignment score of the read alignment that is reported - if the score of the XS tag == AS tag, then the read is a truly multi-mapped read because the reported alignment score is not better than the score of the additional alignments that were found
# - moreover, tophat uses the XS tag to indicate the strand of a read, thus, filtering for the presence of the XS tag without checking the score will lead to over-deletion of reads
# - a third tag, NH, should also be checked - it should be set to 1 if only one alignment was found

# version 2: added regional filtering with BED file

import sys
import argparse

# bx python
from bx.intervals.intersection import IntervalTree, Interval

def BED_to_interval_tree(BED_file):
    """
    Creates an index of intervals for each BED entri

    :param BED_file: file handler of a BED file
    """
    from bx.intervals.intersection import IntervalTree, Interval

    bed_interval_tree = {}
    for line in BED_file:
        if line[0] == "#": continue
        fields = line.strip().split()
        chrom, start_bed, end_bed, = fields[0], int(fields[1]), int(fields[2])

        if chrom not in bed_interval_tree:
            bed_interval_tree[chrom] = IntervalTree()

        # skip if a region overlaps with a region already seen
        """
        if len(bed_interval_tree[chrom].find(start_bed, start_bed + 1)) > 0:
            continue
        """

        bed_interval_tree[chrom].add_interval(Interval(start_bed, end_bed))

#    print "finish proccessing bed file"
    return bed_interval_tree


def main():
    p = argparse.ArgumentParser(
        description = "This script filters sam-files generated by "
        "BOWTIE2 (!) according to user-specified options. It was "
        "initially written by Fidel and modified by Sarah & Friederike. "
        "At least one option has to be chosen, otherwise there will be no "
        "output. It should be used in a pipeline as follows\n /package/samtools/samtools view -h - | %(prog)s",
        usage = "/package/samtools/samtools view -h - | %(prog)s")

    p.add_argument('--remove_blklist',
                 help="BED file containing regions to filter out. Usually a black list "
                   "regions want to be filter out.",
                 type=argparse.FileType('r'))
    p.add_argument('--random', action="store_true",
                 help="if set, all reads mapping to 'random' chromosome are removed")
    p.add_argument('--chrM', action="store_true",
                 help = "if set, all reads mapping to the mitochondrial chromosome are removed")
    p.add_argument('--multiple', action="store_true",
                 help = "if set, all reads where more than one alignment "
                 "was reported by bowtie2 are removed")

    args = p.parse_args()

    if args.remove_blklist:
        filter_out = BED_to_interval_tree(args.filter_out_from_BED)

    for line in sys.stdin:

        field = line.split("\t")

        # filter header - skips random and/or chrM depending on the chosen option
        if line[0] == '@':
    	    if (args.random and field[1][-7:] == '_random') or (args.chrM and field[1] == 'SN:chrM'):
    	        continue
    	    else:
    	        print line,
    	        continue

	    # skip random and mitochondrial chromosomes (reads)
        if (args.random and field[2][-7:] == '_random') or (args.chrM and field[2] == 'chrM'):
    	    continue

        # filtering based on a bed file
        if args.remove_blklist and field[2] in filter_out:
            pos = int(field[3])
            match = filter_out[field[2]].find(pos-300, pos + 300)
            if len(match)>0:
                continue


        # remove reads with more than one reported alignment (indicated by the presence of <XS:i:> in the sam file)
        if args.multiple and len(field) >= 12:
            alignScore = None
            additionalAlignScore = None
            alignNumbers = None
            for subfield in field[12:]:
                (name,ftype,value) = subfield.split(":")
                if name=='AS':
                    alignScore = int(value)
                if name=='XS':
                    additionalAlignScore = int(value)
                if name=='NH':
                    alignNumbers = int(value)
            if additionalAlignScore != None and alignScore == additionalAlignScore:
                continue
            elif alignNumbers > 1:
                continue

            line = "\t".join(field)

        print line,

if __name__ == '__main__':
    main()
