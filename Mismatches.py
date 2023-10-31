import pandas


class Mismatches:
    def __init__(self):
        self.all_mismatches = {}

    def add_mismatch(self, read_id, genomic_range):
        self.all_mismatches[read_id] = [genomic_range.chr_name, genomic_range.mismatch_start,
                                        genomic_range.mismatch_end, genomic_range.mismatch_type,
                                        genomic_range.name, genomic_range.strand]

    def create_csv(self):
        mismatches_data_frame = pandas.DataFrame.from_dict(self.all_mismatches, orient='index')
        mismatches_data_frame.to_csv('mismatches.bed', index=False, sep='\t', header=False)
