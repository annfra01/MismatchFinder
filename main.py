import argparse
import os

import pysam


# Press Umschalt+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", dest='BAM')
    args = parser.parse_args()

    return args


def read_bam(file_name):
    count_mismatches = 0
    bam_file = pysam.AlignmentFile(file_name, "rb")
    # print(bam_file)

    print(bam_file.references)
    print(bam_file.lengths)
    end_pos = 0
    for read in bam_file.fetch():
        no_mismatches = read.get_tag(tag="nM", with_value_type=True)[0]
        mismatch_type = ""
        cigar_string = read.cigarstring
        cigar_tuples = read.cigartuples
        strand = get_strand(read)
        start_pos = read.pos
        chromosome = read.reference_name
        if no_mismatches == 0:
            continue
        if no_mismatches == 1:
            # print("start-pos before: ", start_pos)
            count_mismatches += 1
            if 'I' in cigar_string:
                mismatch_type = "Insertion"
                end_pos, start_pos = get_in_del_pos(cigar_tuples, start_pos)
            elif 'D' in cigar_string:
                mismatch_type = "Deletion"
                end_pos, start_pos = get_in_del_pos(cigar_tuples, start_pos)
        elif no_mismatches == 2:
            count_mismatches += 2
            end_pos = start_pos + 2
            if 'I' in cigar_string:
                mismatch_type = "Insertion"
                # print(cigar_tuples[0])
            elif 'D' in cigar_string:
                mismatch_type = "Deletion"
                # print(cigar_tuples[0])
        # print("start-pos after: ", start_pos)
        # if 'I' in cigar_string:
        print(chromosome,
              start_pos,
              end_pos,
              # no_mismatches,
              mismatch_type,
              # read.cigarstring,
              # read.cigartuples,
              strand,
              sep='\t')
        # going through cigartuple if mismatch is deletion or insertion
        # print(start_pos)
        # for x in range(len(cigar_tuples)):
        #    type_mismatch = cigar_tuples[x][0]
    #     length_mismatch_match = cigar_tuples[x][1]
    #    if type_mismatch == 0:
    #        start_pos += length_mismatch_match
    #    elif type_mismatch == 1 or type_mismatch == 2:
    #        break
    # print(start_pos)

    # print(read.cigarstring)
    # print(read.get_tags(with_value_type=True))
    # print(0 in read.get_tag(tag="nM", with_value_type=True))
    print("total mismatches:", count_mismatches) # different number to grep Befehl
    # bam_file.close()


""" get position of Insertion or Deletion 
    Args: cigar tuples, start position of read
    Returns: start and end postion of Insertion or Deletion
"""


def get_in_del_pos(cigar_tuples, start_pos):
    for x in range(len(cigar_tuples)):
        type_mismatch_cigar = cigar_tuples[x][0]
        length_mismatch_match = cigar_tuples[x][1]
        if type_mismatch_cigar != 1 and type_mismatch_cigar != 2:
            start_pos += length_mismatch_match
        elif type_mismatch_cigar == 1 or type_mismatch_cigar == 2:
            break
    end_pos = start_pos + 1
    return end_pos, start_pos


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


def main():
    print(parse_args())
    bam_file = parse_args().BAM
    if not os.path.exists(bam_file):
        raise OSError("Could not find {}.".format(bam_file))  # doesnt work yet

    # with pysam.AlignmentFile(bam_file, "r") as bam:
    read_bam(bam_file)


if __name__ == '__main__':
    main()
