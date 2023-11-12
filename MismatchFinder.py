import argparse
import os
import re
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
        no_mismatches = read.get_tag(tag="nM", with_value_type=True)[0]
        mismatch_type = ""
        cigar_string = read.cigarstring
        cigar_tuples = read.cigartuples
        current_md_cigar_pos = 0
        current_pos = read.pos
        if no_mismatches == 0:
            continue
        for mismatch in range(no_mismatches):
            genomic_range = GenomicRange()
            genomic_range.chr_name = read.reference_name
            genomic_range.strand = get_strand(read)
            if no_mismatches == 2:
                print(no_mismatches)
            count_mismatches += 1
            if 'I' in cigar_string:
                genomic_range.mismatch_type = "Insertion"
                # length_mismatch_match = cigar_tuples[x][1]
                current_pos, current_md_cigar_pos = get_in_del_pos(cigar_tuples, current_pos, genomic_range, current_md_cigar_pos)
            elif 'D' in cigar_string:
                genomic_range.mismatch_type = "Deletion"
                current_pos, current_md_cigar_pos = get_in_del_pos(cigar_tuples, current_pos, genomic_range, current_md_cigar_pos)
            elif 'M' in cigar_string:
                md_tag = read.get_tag(tag="MD", with_value_type=True)[0]
                md_sub = re.sub(r'([\\^]*[ACGT]+)[0]*', ' \\1 ', md_tag)
                md_split = re.split('[ ]+', md_sub)
                genomic_range.mismatch_type = "Mismatch"
                current_pos, current_md_cigar_pos = get_mismatch_pos(md_split, current_pos, genomic_range, current_md_cigar_pos)
            mismatches.add_mismatch(read.query_name, genomic_range)
            #print(genomic_range.chr_name,
            #      read_start,
            #      genomic_range.mismatch_start,
            #      # end_pos,
            #      genomic_range.mismatch_end,
            #      no_mismatches,
            #      genomic_range.mismatch_type,
            #      genomic_range.name,
            #      read.cigarstring,
            #      # read.cigartuples,
            #      genomic_range.strand,
            #      sep='\t')
    print("total mismatches:", count_mismatches)  # different number to grep Befehl
    bam_file.close()


def get_mismatch_pos(md_split, current_pos, genomic_range, current_md_pos):
    genomic_range.mismatch_start = current_pos
    for x in range(current_md_pos, len(md_split)):
        current_md_pos += 1
        # print(md_split[x])
        current_element = md_split[x]
        if not current_element.isdigit():
            genomic_range.mismatch_start += 1
            break
        else:
            genomic_range.mismatch_start += int(current_element)
    genomic_range.mismatch_end = genomic_range.mismatch_start + 1
    return genomic_range.mismatch_start, current_md_pos


""" get position of Insertion or Deletion 
    Args: cigar tuples, start position of read
    Returns: start and end postion of Insertion or Deletion
"""


def get_in_del_pos(cigar_tuples, read_start, genomic_range, current_cigar_pos):
    for x in range(current_cigar_pos, len(cigar_tuples)):
        current_cigar_pos += 1
        type_mismatch_cigar = cigar_tuples[x][0]
        length_mismatch_match = cigar_tuples[x][1]
        if type_mismatch_cigar != 1 and type_mismatch_cigar != 2:
            genomic_range.mismatch_start = read_start + length_mismatch_match
        elif type_mismatch_cigar == 1 or type_mismatch_cigar == 2:
            genomic_range.mismatch_start += 1
            break
    genomic_range.mismatch_end = genomic_range.mismatch_start + 1
    return genomic_range.mismatch_start, current_cigar_pos


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
