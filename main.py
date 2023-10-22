# This is a sample Python script.
import argparse
import pysam


# Press Umschalt+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file")
    args = parser.parse_args()
    return args.file


def read_bam(file_name):
    bam_file = pysam.AlignmentFile(file_name)
    print(bam_file.header)
    print(bam_file)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    read_bam(get_arguments())

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
