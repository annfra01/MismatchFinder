import argparse
import os
import re
import pysam
from GenomicRange import GenomicRange
from Mismatches import Mismatches

global args

""" gets arguments (bam file and index file) from command line
"""


def parse_args():
    parser = argparse.ArgumentParser(prog="MismatchFinder",
                                     description="This tool finds mismatches in genomic sequences.")
    parser.add_argument("-f", "--file", dest='BAM', help='This file needs to be a .bam file.')
    parser.add_argument("-i", "--index", dest='INDEX', help='This file needs to be a .bai file.')
    global args
    args = parser.parse_args()


""" reads bam-file and finds all mismatches 
    :arg    bam-file, index-file, object Mismatches
"""


def find_mismatches(file_name, index_file, mismatches):
    count_mismatches = 0
    bam_file = pysam.AlignmentFile(file_name, "rb", index_filename=index_file)
    count_reads = 0
    count_reads_with_mismatch = 0
    for read in bam_file.fetch():
        count_reads += 1
        no_all_mismatches = read.get_tag(tag="NM", with_value_type=True)[0]
        no_reg_mismatches = read.get_tag(tag="nM", with_value_type=True)[0]
        cigar_string = read.cigarstring
        cigar_tuples = read.cigartuples
        current_md_cigar_pos = 0
        current_pos = read.reference_start  # instead of read.pos to skip soft clipping
        has_deletion = 'D' in cigar_string
        has_insertion = 'I' in cigar_string
        if no_all_mismatches == 0:
            continue
        count_reads_with_mismatch += 1
        if has_insertion:
            count_mismatches += 1
            genomic_range = create_genomic_range(read)
            genomic_range.mismatch_type = "Insertion"
            # length_mismatch_match = cigar_tuples[x][1]
            get_in_del_pos(cigar_tuples, current_pos, genomic_range,
                           current_md_cigar_pos, has_deletion)
            mismatches.add_mismatch(read.query_name, genomic_range)
            # print(genomic_range.chr_name,
            #      read.query_name,
            #      #      read_start,
            #      genomic_range.mismatch_start,
            #      #      # end_pos,
            #      genomic_range.mismatch_end,
            #      #      no_mismatches,
            #      genomic_range.mismatch_type,
            #      #      genomic_range.name,
            #      read.cigarstring,
            #      #      # read.cigartuples,
            #      #      genomic_range.strand,
            #      sep='\t')
        elif has_deletion:
            count_mismatches += 1
            genomic_range = create_genomic_range(read)
            genomic_range.mismatch_type = "Deletion"
            get_in_del_pos(cigar_tuples, current_pos, genomic_range, current_md_cigar_pos, has_deletion)
            mismatches.add_mismatch(read.query_name, genomic_range)
            # print(genomic_range.chr_name,
            #      read.query_name,
            #      #      read_start,
            #      genomic_range.mismatch_start,
            #      #      # end_pos,
            #      genomic_range.mismatch_end,
            #      #      no_mismatches,
            #      genomic_range.mismatch_type,
            #      #      genomic_range.name,
            #      read.cigarstring,
            #     #      # read.cigartuples,
            #      #      genomic_range.strand,
            #      sep='\t')
        for mismatch in range(no_reg_mismatches):
            genomic_range = create_genomic_range(read)
            genomic_range.mismatch_type = "Mismatch"
            count_mismatches += 1
            md_tag = read.get_tag(tag="MD", with_value_type=True)[0]
            md_sub = re.sub(r'([ACGT]+)', ' \\1 ', md_tag)
            md_split = re.split(' ', md_sub)
            current_pos, current_md_cigar_pos = get_mismatch_pos(md_split, current_pos, genomic_range,
                                                                 current_md_cigar_pos, cigar_tuples)
            mismatches.add_mismatch(read.query_name, genomic_range)
            # print(genomic_range.chr_name,
            #      read.query_name,
            #      #      read_start,
            #      genomic_range.mismatch_start,
            #      #      # end_pos,
            #      genomic_range.mismatch_end,
            #      #      no_mismatches,
            #      genomic_range.mismatch_type,
            #      #      genomic_range.name,
            #      read.cigarstring,
            #      #      # read.cigartuples,
            #      #      genomic_range.strand,
            #      sep='\t')
    print("total mapped reads:", count_reads)
    print("total mapped reads with mismatches:", count_reads_with_mismatch)
    print("total mismatches:", count_mismatches)  # different number to grep Befehl
    bam_file.close()


""" creates new GenomicRange object
"""


def create_genomic_range(read):
    genomic_range = GenomicRange()
    genomic_range.chr_name = read.reference_name
    genomic_range.strand = get_strand(read)
    return genomic_range


""" gets regular mismatch position
    :arg    MD-tag splitted, current position in reference, genomic range object, current position
            in MD-tag, cigar tuple
    :returns mismatch start as new current position in reference, current position in MD-tag
"""


def get_mismatch_pos(md_split, current_pos, genomic_range, current_md_pos, cigar_tuples):
    current_cigar_pos = current_pos
    current_cigar_pos += cigar_tuples[0][1]
    current_cigar_pos_count = 0
    # if cigar_tuples[0][0] == 4:
    #    genomic_range.mismatch_start = current_pos
    genomic_range.mismatch_start = current_pos
    for x in range(current_md_pos, len(md_split)):
        current_md_pos += 1
        # print(md_split[x])
        current_element = md_split[x]
        current_element = current_element.split("^")[0]
        #
        # print(current_element)
        if not current_element.isdigit() and len(current_element) == 1:
            genomic_range.mismatch_start += 1
            break
        elif not current_element.isdigit() and len(current_element) > 1:
            genomic_range.mismatch_start += len(current_element)
        else:
            genomic_range.mismatch_start += int(current_element)
            if current_cigar_pos < genomic_range.mismatch_start:
                current_cigar_pos_count += 1
                if cigar_tuples[current_cigar_pos_count][0] == 3:
                    genomic_range.mismatch_start += cigar_tuples[current_cigar_pos_count][1]
    genomic_range.mismatch_end = genomic_range.mismatch_start + 1
    return genomic_range.mismatch_start, current_md_pos


""" gets position of Insertion or Deletion 
    :arg cigar tuples, start position of read
    :returns start and end postion of Insertion or Deletion
"""


def get_in_del_pos(cigar_tuples, read_start, genomic_range, current_cigar_pos, has_deletion):
    for x in range(current_cigar_pos, len(cigar_tuples)):
        current_cigar_pos += 1
        type_mismatch_cigar = cigar_tuples[x][0]
        length_mismatch_match = cigar_tuples[x][1]
        if type_mismatch_cigar != 1 and type_mismatch_cigar != 2:
            genomic_range.mismatch_start = read_start + length_mismatch_match
        elif type_mismatch_cigar == 1 or type_mismatch_cigar == 2:
            if has_deletion:
                genomic_range.mismatch_start += 1
            break
    genomic_range.mismatch_end = genomic_range.mismatch_start + length_mismatch_match
    # return genomic_range.mismatch_start, current_cigar_pos


""" get direction of read; if flag in SAM file is 0, then forward strand
                           if flag in SAM file is 16, then reverse strand
    :arg current read
    :returns read direction (+ or -)
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
    bam_file = args.BAM
    index_file = args.INDEX
    if not os.path.exists(bam_file):
        raise OSError("Could not find {}.".format(bam_file))  # doesnt work yet
    find_mismatches(bam_file, index_file, mismatches)
    mismatches.create_csv()


if __name__ == '__main__':
    mismatch_finder()
