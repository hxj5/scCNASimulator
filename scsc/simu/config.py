# config.py - global configuration


import sys
from ..utils.region import format_chrom


class Config:
    def __init__(self):
        self.defaults = DefaultConfig()

        self.sam_fn = None
        self.out_dir = None

        self.cell_anno_fn = None
        self.ref_cell_types_str = None
        self.ref_cell_types = None

        self.cnv_profile_fn = None
        self.feature_fn = None

        self.baf_dir = None
        self.umi_dir = None

        self.cell_tag = self.defaults.CELL_TAG
        self.umi_tag = self.defaults.UMI_TAG

        self.debug = self.defaults.DEBUG

        self.merged_cnv_profile_fn = None

        self.baf_fn_prefix = "scsc"
        self.umi_fn_prefix = "allele_umi"
        self.umi_fn_suffix = ".tsv"

        self.out_sam_fn = None
        self.out_cnv_fn = None
        self.out_feature_fn = None
        self.out_cell_anno_fn = None
        self.out_umi_stat_fn = None


    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stderr

        s =  "%s\n" % prefix
        s += "%ssam_file = %s\n" % (prefix, self.sam_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)

        s += "%scell_anno_fn = %s\n" % (prefix, self.cell_anno_fn)
        s += "%sref_cell_types = %s\n" % (prefix, self.ref_cell_types_str)

        s += "%scnv_profile_fn = %s\n" % (prefix, self.cnv_profile_fn)
        s += "%sfeature_fn = %s\n" % (prefix, self.feature_fn)

        s += "%sbaf_dir = %s\n" % (prefix, self.baf_dir)
        s += "%sumi_dir = %s\n" % (prefix, self.umi_dir)

        s += "%scell_tag = %s\n" % (prefix, self.cell_tag)
        s += "%sumi_tag = %s\n" % (prefix, self.umi_tag)
        s += "%sdebug = %d\n" % (prefix, self.debug)
        s += "%s\n" % prefix

        s += "%sbaf_fn_prefix = %s\n" % (prefix, self.baf_fn_prefix)
        s += "%sumi_fn_prefix = %s\n" % (prefix, self.umi_fn_prefix)
        s += "%sumi_fn_suffix = %s\n" % (prefix, self.umi_fn_suffix)

        s += "%sout_sam_file = %s\n" % (prefix, self.out_sam_fn)
        s += "%sout_cnv_file = %s\n" % (prefix, self.out_cnv_fn)
        s += "%sout_feature_file = %s\n" % (prefix, self.out_feature_fn)
        s += "%sout_cell_anno_file = %s\n" % (prefix, self.out_cell_anno_fn)
        s += "%sout_umi_stat_file = %s\n" % (prefix, self.out_umi_stat_fn)
        s += "%s\n" % prefix

        s += "%s\n" % prefix

        fp.write(s)


class DefaultConfig:
    def __init__(self):
        self.CELL_TAG = "CB"
        self.UMI_TAG = "UB"
        self.UMI_TAG_BC = "UB"    # the default umi tag for 10x data.
        self.DEBUG = 0
