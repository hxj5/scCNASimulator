# config.py - global configuration


import sys
from ..blib.region import format_chrom


class Config:
    def __ini__(self):
        self.bam_fn = None
        self.out_dir = None

        self.cell_anno_fn = None
        self.cnv_profile_fn = None
        self.umi_dir = None

        self.cell_tag = None
        self.umi_tag = None

        self.cell_anno = None
        self.cnv_profile = None

    def show(self, fp = sys.stdout):
        pass

