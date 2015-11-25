#!/usr/bin/env python

'''
This program filters and splits BAM files processed by suspenders:
Originally the split script was written by @friedue which I modified extensively.
Reads are written into different output alignment files (.bam) depending on
their genomic origin (which is noted in the po:i tag). If a third output file
is indicated, reads which map equally well at the same position in both genomes
will also be reported (po:i:3). Currently, this script ignores reads that were
assigned by suspender's "RANDOM" filter'.
Possible improvements:
let user choose which reads should be ignored (not just
RANDOM reads, perhaps);
allow SAM format for in- and output;
multiprocessing
'''

import sys
import argparse
import getopt
import pysam
import os
import doctest
import re
# bx python (bed file handler for filtering out blklist )
from bx.intervals.intersection import IntervalTree, Interval

def get_args():
    parser=argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve',
        description='Filter or Split suspenders-merged BAM files')

    subparsers = parser.add_subparsers(
        title="commands",
        dest='command',
        metavar='')

    ## I am splitting this into two scripts : FILTER and SPLIT
    filter_mode = subparsers.add_parser(
        'filter',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Filter the BAM file with multiple criteria.",
        usage='%(prog)s filter '
              '--remove_blklist /path/to/blklist-regions.bed '
              '--random '
              '-out filtered.bam \n')
    ## For Filtering
    filter_mode.add_argument('--remove_blklist',
                 help="BED file containing regions to filter out. Usually a black list "
                   "regions want to be filter out.",
                 type=argparse.FileType('r'))
    filter_mode.add_argument('--random', action="store_true",
                 help="if set, all reads mapping to 'random' chromosome are removed")
    filter_mode.add_argument('--chrM', action="store_true",
                 help = "if set, all reads mapping to the mitochondrial chromosome are removed")
    filter_mode.add_argument('--multiple', action="store_true",
                 help = "if set, all reads where more than one alignment "
                 "was reported by bowtie2 are removed")

    #args = p.parse_args()

    split_mode = subparsers.add_parser(
        'split',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Splits the BAM file by the reads mapping to maternal/paternal allele (or equally mapped)",
        usage='%(prog)s '
              '--BAMfile filtered.bam '
              '--outfile1 filtered_maternal.bam'
              '--outfile2 filtered_paternal.bam'
              '--outfile3 filtered_equalmapped.bam')

    ## For splitting
    split_mode.add_argument('--BAMfile', '-in', type=str, required=True, help="suspender-generated BAM file")
    split_mode.add_argument('--outfile1', '-alt1', type=str, required=True, help = 'Prefix for output file with alternative genome 1')
    split_mode.add_argument('--outfile2', '-alt2', type=str, required=True, help = 'Prefix for output file with alternative genome 2')
    split_mode.add_argument('--outfile3', '-neither', type=str, help = 'Prefix for output file with reads that mapped equally well to both alternative genomes (optional).')
    split_mode.add_argument('--filterDuplicates', '-dup', type=str, default='yes', help = 'remove duplicates; default: yes')
    split_mode.add_argument('--filterForMappingQuality', '-mapQ', type=int, help = 'indicate a minimum mapping quality that each read (pair) should have, default: 0')
    split_mode.add_argument('--removeMultiMapped', '-vMulti', action='store_true', default=False, help = 'exclude reads that have more than one alignment based on the NH and/or XS tag; default: not set')
    split_mode.add_argument('--coordinateSorting', '-coordSort', action='store_true', default=True, help = 'sort the resulting bam files according to read coordinates and index them; default: not set')

    args=parser.parse_args()
    return args

##### ALRIGHT, All arguments parsed.. NOW WE WRITE THE FUNCTIONS

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

    return bed_interval_tree


def prepare_header(BamFile, Sorting):
    '''replaces the entry that indicates the sorting of a BAM file with a PICARD-consistent entry'''
    replacements = {'query' : 'queryname', 'unsort' : 'unsorted' , 'coordinate' : 'coordinate'}
    headerList = ['query', 'coordinate','unsort']
    oH = BamFile.header

    if Sorting:
        oH['HD']['SO'] = 'coordinate'
    else:
        for k in headerList:
            if re.search(k, oH['HD']['SO']):
                oH['HD']['SO'] = replacements[k]

    return(oH)


def remove_multiAlignments(Read, KeepInfo):
    '''Checks for NH and XS tag to remove reads that aligned more than once,
    returns a boolean whether the read should be kept or not
    '''
    keepRead = KeepInfo

    # the highest possible alignment score in end-to-end align.bamment is 0.
    # For each mismatch or gap open/extension,
    # that score is decremented by some amount --> negative values
    try:
        alignScore = [item for item in Read.tags if item[0] == 'AS'][0][1]
    except IndexError:
        alignScore = None

    # XS --> tophat uses it to indicate which strand the read originated from (+/-)
    # the XS flag is set by bowtie2, when a second valid alignment for the same read is found
    # if the alignment score of the reported read AS==XS,
    # then that read is multimapping
    # caution: bowtie2 < version 2.0. had a bug in setting XS tags for paired-end data
    try:
        additionalAlignScore = [item for item in Read.tags if item[0] == 'XS'][0][1]
    except IndexError:
        additionalAlignScore = None

    # NH:i tag should indicate how many alignments were found
    # it is not set by bowtie2 or bwa
    # not sure if this is a good filter
    try:
        alignNumbers = [item for item in Read.tags if item[0] == 'NH'][0][1]
    except IndexError:
        alignNumbers = None


    if isinstance(additionalAlignScore, (int, long, float, complex)) and alignScore == additionalAlignScore:
        keepRead=False
    elif alignNumbers > 1:
        keepRead=False

    return(keepRead)

def check_mappingQuality(Read, KeepInfo, MinMapQ):
    '''Returns boolean based on the mapping quality filter indicated by the user'''
    keepRead = KeepInfo
    if Read.mapq < MinMapQ:
        keepRead = False

    return(keepRead)


def save_reads(Read, KeepInfo, outFile1, outFile2, outFile3, c):
    '''Checks the origin of a read and saves it to the respective file'''
    ReadOrigin = [item for item in Read.tags if item[0] == 'po'][0][1]

    if(ReadOrigin == 1 and KeepInfo):
            outFile1.write(Read)
            c[0] += 1
    elif(ReadOrigin == 2 and KeepInfo):
            outFile2.write(Read)
            c[1] += 1
    elif KeepInfo:
            c[2] += 1
            if outFile3 is not None:
                outFile3.write(Read)
    return(c)


def sort_output(outPrefix):
    '''Sorts the output file by read coordinate'''
    pysam.sort(outPrefix+'.originalSort.bam', outPrefix + '.coordSort')
    ## Build the bam index for output
    pysam.index(outPrefix + '.coordSort.bam')
    ## Remove original sorted file
    os.remove(outPrefix+'.originalSort.bam')

def main():

    args = get_args()
    ## IF ARGS==FILTER THEN DO THE FILTERING STUFF

    if args.command == 'filter':
        # adding samfilter functions first
        if args.remove_blklist:
            filter_out = BED_to_interval_tree(args.remove_blklist)

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


    ## ELSE DO THE SPLITTING STUFF
    if args.command == 'split':
        infile = pysam.Samfile(args.BAMfile, "rb")
        newHeader = prepare_header(infile, args.coordinateSorting)

        out1 = args.outfile1
        out2 = args.outfile2
        out1_noSort = pysam.Samfile(out1 + '.originalSort.bam', "wb", header = newHeader ) # template will cause the original header to be used, one could modify the header
        # by generating an appropriate dictionary structure, newHeader, and then use that: header = newHeader
        out2_noSort = pysam.Samfile(out2 + '.originalSort.bam' , "wb", header = newHeader )
        if args.outfile3:
            out3 = args.outfile3
            out3_noSort = pysam.Samfile(out3 + '.originalSort.bam' , "wb", header = newHeader )

        else:
            out3 = None

        keep = True
        origin = 0
        Counts = [0,0,0] # keeping track of reads per origin

        # seen is a dictionary of read names
        # the value associated is a boolean
        # that when true, it means to keep the read
        seen = {}
        seenRead = {}

        for DNAread in infile:
            #check if read should be kept
            if DNAread.flag == 4: # excluding unmapped reads
                keep = False

            if args.filterDuplicates == 'yes' and DNAread.is_duplicate:
                keep=False

            if keep and [item for item in DNAread.tags if item[0] == 'ct'][0][1] == 'R':
                keep = False

            if keep and args.removeMultiMapped:
                keep = remove_multiAlignments(DNAread, True)

            if keep and args.filterForMappingQuality:
                keep = check_mappingQuality(DNAread, True, args.filterForMappingQuality)

            # if the current read is the second mate,
            # save both reads (if they passed)
            # and delete the entries from the dictionaries
            if DNAread.qname in seen:
                keep1 = seen[DNAread.qname]
                keep2 = keep
            # if either of the reads can be kept, keep both
                if keep1 or keep2:
                    Counts = save_reads(DNAread, True, outFile1 = out1_noSort, outFile2 = out2_noSort, outFile3 = out3_noSort, c = Counts)
                    Counts = save_reads(seenRead[DNAread.qname], True, outFile1 = out1_noSort, outFile2 = out2_noSort, outFile3 = out3_noSort, c = Counts)
                del(seen[DNAread.qname])
                del(seenRead[DNAread.qname])

            # if the current read is not within the seen dictionary
            # first check if the read is unpaired
            # if unpaired, save and forget
            # if paired, save and wait for mate
            elif DNAread.mate_is_unmapped or not DNAread.is_paired:
                Counts = save_reads(DNAread, keep, outFile1 = out1_noSort, outFile2 = out2_noSort, outFile3 = out3_noSort, c = Counts)
            elif DNAread.is_proper_pair:
                seen[DNAread.qname] = keep
                seenRead[DNAread.qname] = DNAread

            #origin = 0
            keep = True

        for entry in seen:
            Counts = save_reads(seenRead[entry], seen[entry], outFile1 = out1_noSort, outFile2 = out2_noSort, outFile3 = out3_noSort, c = Counts)
            keep = True

        out1_noSort.close()
        out2_noSort.close()
        if args.outfile3:
            out3_noSort.close()

        CountOut='''
        reads mapped to alternative genome 1 ({1}): {0}\n
        reads mapped to alternative genome 2 ({3}): {2}\n
        undistinguishable reads {4}'''.format(Counts[0], args.outfile1, Counts[1], args.outfile2, Counts[2])

        print CountOut

        if args.coordinateSorting:
            sort_output(out1)
            sort_output(out2)
            if args.outfile3:
                sort_output(out3)

if __name__ == '__main__':
    main()
