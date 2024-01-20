# config.py - global configuration


import sys
from ..blib.region import format_chrom


class Config:
    def __ini__(self):
        self.sam_fn = None
        self.out_dir = None

        self.cell_anno_fn = None
        self.cnv_profile_fn = None
        self.umi_dir = None

        self.cell_tag = None
        self.umi_tag = None

        self.merged_cnv_profile_fn = None

        self.umi_fn_prefix = "allele_umi"
        self.umi_fn_suffix = ".tsv"

        self.out_sam_fn = None


    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stderr

        s =  "%s\n" % prefix
        s += "%ssam_file = %s\n" % (prefix, self.sam_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%scell_anno_fn = %s\n" % (prefix, self.cell_anno_fn)
        s += "%scnv_profile_fn = %s\n" % (prefix, self.cnv_profile_fn)
        s += "%sumi_dir = %s\n" % (prefix, self.umi_dir)
        s += "%scell_tag = %s\n" % (prefix, self.cell_tag)
        s += "%sumi_tag = %s\n" % (prefix, self.umi_tag)
        s += "%s\n" % prefix

        s += "%smerged_cnv_profile_fn = %s\n" % (prefix, self.merged_cnv_profile_fn)
        s += "%sumi_fn_prefix = %s\n" % (prefix, self.umi_fn_prefix)
        s += "%sumi_fn_suffix = %s\n" % (prefix, self.umi_fn_suffix)
        s += "%sout_sam_file = %s\n" % (prefix, self.out_sam_fn)
        s += "%s\n" % prefix

        s += "%s\n" % prefix

        fp.write(s)

