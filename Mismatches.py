import pandas


class Mismatches:
    def __init__(self):
        self.all_mismatches = {}  # includes a list of lists in case there is more than one mismatch per read

    """ adds one mismatch to dictionary which stores all mismatches
    """

    def add_mismatch(self, read_id, genomic_range):
        if read_id in self.all_mismatches:
            self.all_mismatches[read_id].append([genomic_range.chr_name, genomic_range.mismatch_start,
                                                 genomic_range.mismatch_end, genomic_range.mismatch_type,
                                                 genomic_range.name, genomic_range.strand])
        else:
            self.all_mismatches[read_id] = [[genomic_range.chr_name, genomic_range.mismatch_start,
                                             genomic_range.mismatch_end, genomic_range.mismatch_type,
                                             genomic_range.name, genomic_range.strand]]

    def dict_to_list(self):
        mismatches_list = []
        for read in self.all_mismatches:
            for mismatch in (self.all_mismatches[read]):
                mismatches_list.append(mismatch)
        return mismatches_list

    """ creates csv from dict with all mismatches
    """

    def create_csv(self):
        mismatches_data_frame = pandas.DataFrame.from_dict(self.all_mismatches, orient='index')
        mismatches_data_frame.to_csv('example.bed', index=False, sep='\t', header=False)

    def create_csv_from_list(self):
        mismatches_data_frame = pandas.DataFrame(self.dict_to_list())
        mismatches_data_frame.to_csv('mismatches.bed', index=False, sep='\t', header=False)
