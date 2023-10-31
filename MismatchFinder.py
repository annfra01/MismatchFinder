import argparse
import os
import pandas
import pysam
from GenomicRange import GenomicRange
from Mismatches import Mismatches

# Press Umschalt+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

global args

""" get arguments from command line
"""


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", dest='BAM')
    parser.add_argument("-i", "--index", dest='INDEX')
    parser.add_argument("-r", "--reference", dest='REF')
    global args
    args = parser.parse_args()


def read_bam(file_name, index_file, mismatches):
    count_mismatches = 0
    bam_file = pysam.AlignmentFile(file_name, "rb", index_filename=index_file)
    print(bam_file.references)
    print(bam_file.lengths)
    end_pos = 0
    for read in bam_file.fetch():
        genomic_range = GenomicRange()
        no_mismatches = read.get_tag(tag="nM", with_value_type=True)[0]
        mismatch_type = ""
        cigar_string = read.cigarstring
        cigar_tuples = read.cigartuples
        genomic_range.strand = get_strand(read)
        read_start = read.pos
        genomic_range.chr_name = read.reference_name
        if no_mismatches == 0:
            continue
        # if no_mismatches == 1:
        # print("start-pos before: ", read_start)
        # count_mismatches += 1
        # if 'I' in cigar_string:
        #    mismatch_type = "Insertion"
        #    end_pos, read_start = get_in_del_pos(cigar_tuples, read_start)
        # elif 'D' in cigar_string:
        #    mismatch_type = "Deletion"
        #    end_pos, read_start = get_in_del_pos(cigar_tuples, read_start)
        # elif 'M' in cigar_string:

        # elif no_mismatches == 2:
        # count_mismatches += 2
        # end_pos = read_start + 2
        # if 'I' in cigar_string:
        #    mismatch_type = "Insertion"
        #    # print(cigar_tuples[0])
        # elif 'D' in cigar_string:
        #    mismatch_type = "Deletion"
        #    # print(cigar_tuples[0])

        for mismatch in range(no_mismatches):
            count_mismatches += 1
            if 'I' in cigar_string:
                genomic_range.mismatch_type = "Insertion"
                get_in_del_pos(cigar_tuples, read_start, genomic_range)
            elif 'D' in cigar_string:
                genomic_range.mismatch_type = "Deletion"
                get_in_del_pos(cigar_tuples, read_start, genomic_range)
            mismatches.add_mismatch(read.query_name, genomic_range)
            # print(genomic_range.chr_name,
            #      # read_start,
            #      genomic_range.mismatch_start,
            #      # end_pos,
            #      genomic_range.mismatch_end,
            #      # no_mismatches,
            #      genomic_range.mismatch_type,
            #      genomic_range.name,
            #      # read.cigarstring,
            #      # read.cigartuples,
            #      genomic_range.strand,
            #      sep='\t')

    print("total mismatches:", count_mismatches)  # different number to grep Befehl
    bam_file.close()


""" get position of Insertion or Deletion 
    Args: cigar tuples, start position of read
    Returns: start and end postion of Insertion or Deletion
"""


def get_in_del_pos(cigar_tuples, read_start, genomic_range):
    for x in range(len(cigar_tuples)):
        type_mismatch_cigar = cigar_tuples[x][0]
        length_mismatch_match = cigar_tuples[x][1]
        if type_mismatch_cigar != 1 and type_mismatch_cigar != 2:
            genomic_range.mismatch_start = read_start + length_mismatch_match
        elif type_mismatch_cigar == 1 or type_mismatch_cigar == 2:
            break
    genomic_range.mismatch_end = genomic_range.mismatch_start + 1


""" get direction of read; if flag in SAM file is 0, then forward strand
                           if flag in SAM file is 16, then reverse strand
    Args: current read
    Returns: read direction (+ or -)
"""


def get_strand(read):
    strand_flag = read.flag  # is_forward oder is_reverse gibts auch
    if strand_flag == 0:
        strand = "+"
    else:
        strand = "-"
    return strand


def mismatch_finder():
    mismatches = Mismatches()
    parse_args()
    global args
    print("args ", args)
    bam_file = args.BAM
    index_file = args.INDEX
    if not os.path.exists(bam_file):
        raise OSError("Could not find {}.".format(bam_file))  # doesnt work yet

    # with pysam.AlignmentFile(bam_file, "r") as bam:
    read_bam(bam_file, index_file, mismatches)
    mismatches.create_csv()


if __name__ == '__main__':
    mismatch_finder()
