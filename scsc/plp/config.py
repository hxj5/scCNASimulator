# Global Configure

import sys
from ..app import APP


class Config:
    def __init__(self):
        self.defaults = DefaultConfig()

        self.sam_fn = None
        self.out_dir = None
        self.region_fn = None
        self.snp_fn = None
        self.barcode_fn = None
        self.debug = self.defaults.PLP_DEBUG

        self.cell_tag = self.defaults.PLP_CELL_TAG
        self.umi_tag = self.defaults.PLP_UMI_TAG
        self.nproc = self.defaults.PLP_NPROC
        self.min_count = self.defaults.PLP_MIN_COUNT
        self.min_maf = self.defaults.PLP_MIN_MAF
        self.output_all_reg = self.defaults.PLP_OUTPUT_ALL_REG
        self.no_dup_hap = self.defaults.PLP_NO_DUP_HAP

        self.min_mapq = self.defaults.PLP_MIN_MAPQ
        self.min_len = self.defaults.PLP_MIN_LEN
        self.incl_flag = self.defaults.PLP_INCL_FLAG
        self.excl_flag = -1
        self.no_orphan = self.defaults.PLP_NO_ORPHAN

        self.sam = None          # a pysam::AlignmentFile object.
        self.barcodes = None     # list of barcode strings.
        self.reg_list = None     # list of gene/block objects.
        self.snp_set = None      # set of SNPs.

        self.out_prefix = APP + "."
        self.out_region_fn = None
        self.out_sample_fn = None
        self.out_ad_fn = None
        self.out_dp_fn = None
        self.out_oth_fn = None   


    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stderr

        s =  "%s\n" % prefix
        s += "%ssam_file = %s\n" % (prefix, self.sam_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%sregion_file = %s\n" % (prefix, self.region_fn)
        s += "%ssnp_file = %s\n" % (prefix, self.snp_fn)
        s += "%sbarcode_file = %s\n" % (prefix, self.barcode_fn)
        s += "%sdebug = %d\n" % (prefix, self.debug)
        s += "%s\n" % prefix

        s += "%scell_tag = %s\n" % (prefix, self.cell_tag)
        s += "%sumi_tag = %s\n" % (prefix, self.umi_tag)
        s += "%snumber_of_processes = %d\n" % (prefix, self.nproc)
        s += "%smin_count = %d\n" % (prefix, self.min_count)
        s += "%smin_maf = %f\n" % (prefix, self.min_maf)
        s += "%soutput_all_reg = %s\n" % (prefix, self.output_all_reg)
        s += "%sno_dup_hap = %s\n" % (prefix, self.no_dup_hap)
        s += "%s\n" % prefix

        s += "%smin_mapq = %d\n" % (prefix, self.min_mapq)
        s += "%smin_len = %d\n" % (prefix, self.min_len)
        s += "%sinclude_flag = %d\n" % (prefix, self.incl_flag)
        s += "%sexclude_flag = %d\n" % (prefix, self.excl_flag)
        s += "%sno_orphan = %s\n" % (prefix, self.no_orphan)
        s += "%s\n" % prefix

        s += "%snumber_of_barcodes = %d\n" % (prefix, len(self.barcodes) if \
                self.barcodes is not None else -1)
        s += "%snumber_of_regions = %d\n" % (prefix, len(self.reg_list) if \
                self.reg_list is not None else -1)
        s += "%snumber_of_snps = %d\n" % (prefix, self.snp_set.get_n() if \
                self.snp_set is not None else -1)
        s += "%s\n" % prefix

        s += "%soutput_region_file = %s\n" % (prefix, self.out_region_fn)
        s += "%soutput_sample_file = %s\n" % (prefix, self.out_sample_fn)
        s += "%soutput_ad_file = %s\n" % (prefix, self.out_ad_fn)
        s += "%soutput_dp_file = %s\n" % (prefix, self.out_dp_fn)
        s += "%soutput_oth_file = %s\n" % (prefix, self.out_oth_fn)
        s += "%s\n" % prefix

        fp.write(s)


    def use_umi(self):
        return self.umi_tag is not None


class DefaultConfig:

    def __init__(self):
        self.PLP_DEBUG = 0
        self.PLP_CELL_TAG = "CB"
        self.PLP_UMI_TAG = "UB"
        self.PLP_UMI_TAG_BC = "UB"    # the default umi tag for 10x data.
        self.PLP_NPROC = 1
        self.PLP_MIN_COUNT = 1 
        self.PLP_MIN_MAF = 0
        self.PLP_OUTPUT_ALL_REG = False
        self.PLP_NO_DUP_HAP = True

        self.PLP_MIN_MAPQ = 20
        self.PLP_MIN_LEN = 30
        self.PLP_INCL_FLAG = 0
        self.PLP_EXCL_FLAG_UMI = 772
        self.PLP_EXCL_FLAG_XUMI = 1796
        self.PLP_NO_ORPHAN = True


if __name__ == "__main__":
    conf = Config()
    conf.show()

