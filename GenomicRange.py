class GenomicRange:
    def __init__(self, chr_name="", mismatch_start=0, mismatch_end=0, mismatch_type="", strand=""):
        self.chr_name = chr_name
        self.mismatch_start = mismatch_start
        self.mismatch_end = mismatch_end
        self.mismatch_type = mismatch_type
        self.name = "."
        self.strand = strand
