# pipeline.py - pipeline running all steps.


import getopt
import os
import sys
from sys import stdout, stderr

from .app import APP, VERSION
from .blib.w_assert import assert_e, assert_n
from .pileup import pileup_main
from .plp.config import Config as PlpConfig
from .simu.cnv import merge_cnv_profile
from .simu.config import Config as SimuConfig
from .simulation import simu_main


class Config:
    def __init__(self):
        self.sam_fn = None
        self.out_dir = None
        self.barcode_fn = None

        self.cell_anno_fn = None
        self.ref_cell_types_str = None

        self.cnv_profile_fn = None
        self.feature_fn = None
        self.snp_fn = None

        self.cell_tag = None
        self.umi_tag = None
        self.nproc = None
        self.min_count = None
        self.min_maf = None
        self.output_all_reg = True
        self.no_dup_hap = True
        self.debug = None

        self.min_mapq = None
        self.min_len = None
        self.incl_flag = None
        self.excl_flag = None
        self.no_orphan = None

        self.merged_cnv_profile_fn = None


def prepare_args(conf):
    assert_e(conf.sam_fn)
    assert_e(conf.barcode_fn)

    assert_n(conf.out_dir)
    if not os.path.exists(conf.out_dir):
        os.mkdir(conf.out_dir)

    assert_e(conf.cell_anno_fn)
    assert_n(conf.ref_cell_types_str)

    assert_e(conf.cnv_profile_fn)
    assert_e(conf.feature_fn)
    assert_e(conf.snp_fn)

    conf.merged_cnv_profile_fn = os.path.join(conf.out_dir, "merged.cnv_profile.tsv")


def conf2plp_argv(conf):
    args = [APP, "pileup"]

    if conf.sam_fn is not None:
        args.extend(["--sam", conf.sam_fn])
    if conf.out_dir is not None:
        plp_dir = os.path.join(conf.out_dir, "pileup")
        args.extend(["--outdir", plp_dir])
    if conf.barcode_fn is not None:
        args.extend(["--barcode", conf.barcode_fn])

    if conf.feature_fn is not None:
        args.extend(["--region", conf.feature_fn])
    if conf.snp_fn is not None:
        args.extend(["--phasedSNP", conf.snp_fn])

    if conf.cell_tag is not None:
        args.extend(["--cellTAG", conf.cell_tag])
    if conf.umi_tag is not None:
        args.extend(["--UMItag", conf.umi_tag])
    if conf.nproc is not None:
        args.extend(["--nproc", conf.nproc])
    if conf.min_count is not None:
        args.extend(["--minCOUNT", conf.min_count])
    if conf.min_maf is not None:
        args.extend(["--minMAF", conf.min_maf])
    if conf.output_all_reg is not None:
        if conf.output_all_reg:
            args.append("--outputAllReg")
    if conf.no_dup_hap is not None:
        if not conf.no_dup_hap:
            args.append("--countDupHap")
    if conf.debug is not None:
        args.extend(["--debug", conf.debug])

    if conf.min_mapq is not None:
        args.extend(["--minMAPQ", conf.min_mapq])
    if conf.min_len is not None:
        args.extend(["--minLEN", conf.min_len])
    if conf.incl_flag is not None:
        args.extend(["--inclFLAG", conf.incl_flag])
    if conf.excl_flag is not None:
        args.extend(["--exclFLAG", conf.excl_flag])
    if conf.no_orphan is not None:
        if not conf.no_orphan:
            args.append("--countORPHAN")

    args = [str(x) for x in args]
    return(args)


def conf2simu_argv(conf, plp_conf):
    args = [APP, "simu"]

    if conf.sam_fn is not None:
        args.extend(["--sam", conf.sam_fn])
    if conf.out_dir is not None:
        simu_dir = os.path.join(conf.out_dir, "simu")
        args.extend(["--outdir", simu_dir])

    if conf.cell_anno_fn is not None:
        args.extend(["--cellAnno", conf.cell_anno_fn])
    if conf.ref_cell_types_str is not None:
        args.extend(["--refCellTypes", conf.ref_cell_types_str])

    if conf.merged_cnv_profile_fn is not None:
        args.extend(["--cnvProfile", conf.merged_cnv_profile_fn])
    if conf.feature_fn is not None:
        args.extend(["--feature", conf.feature_fn])

    assert_e(plp_conf.out_dir)
    args.extend(["--BAFdir", plp_conf.out_dir])
    assert_e(plp_conf.umi_dir)
    args.extend(["--UMIdir", plp_conf.umi_dir])

    args.extend(["--cellTAG", plp_conf.cell_tag])
    args.extend(["--UMItag", plp_conf.umi_tag])
    if conf.debug is not None:
        args.extend(["--debug", conf.debug])

    args = [str(x) for x in args]
    return(args)


def usage(fp = stderr, plp_conf = None, simu_conf = None):
    s =  "\n"
    s += "Version: %s\n" % VERSION
    s += "Usage:   %s %s <options>\n" % (APP, COMMAND)
    s += "\n" 
    s += "Options:\n"
    s += "  -s, --sam FILE          Indexed sam/bam/cram file.\n"
    s += "  -O, --outdir DIR        Output directory for sparse matrices.\n"
    s += "  -b, --barcode FILE      A plain file listing all effective cell barcode.\n"
    s += "      --cellAnno FILE     Cell annotation file, 2 columns <cell> <clone_id>.\n"
    s += "      --refCellTypes STR  Reference cell types, comma separated.\n"
    s += "      --cnvProfile FILE   CNV profile file, 7 columns.\n"
    s += "      --feature FILE      Feature annotation file, typically for genes; 4 columns.\n"
    s += "  -P, --phasedSNP FILE    A TSV or VCF file listing phased SNPs (i.e., containing phased GT).\n"
    s += "  -h, --help              Print this message and exit.\n"
    s += "\n"
    s += "Optional arguments:\n"
    s += "  -p, --nproc INT         Number of processes [%d]\n" % plp_conf.NPROC
    s += "      --cellTAG STR       Tag for cell barcodes [%s]\n" % plp_conf.CELL_TAG
    s += "      --UMItag STR        Tag for UMI, set to None when reads only [%s]\n" % plp_conf.UMI_TAG
    s += "      --minCOUNT INT      Mininum aggragated count for SNP [%d]\n" % plp_conf.MIN_COUNT
    s += "      --minMAF FLOAT      Mininum minor allele fraction for SNP [%f]\n" % plp_conf.MIN_MAF
    s += "  -D, --debug INT         Used by developer for debugging [%d]\n" % plp_conf.DEBUG
    s += "\n"
    s += "Read filtering:\n"
    s += "  --inclFLAG INT    Required flags: skip reads with all mask bits unset [%d]\n" % plp_conf.INCL_FLAG
    s += "  --exclFLAG INT    Filter flags: skip reads with any mask bits set [%d\n" % plp_conf.EXCL_FLAG_UMI
    s += "                    (when use UMI) or %d (otherwise)]\n" % plp_conf.EXCL_FLAG_XUMI
    s += "  --minLEN INT      Minimum mapped length for read filtering [%d]\n" % plp_conf.MIN_LEN
    s += "  --minMAPQ INT     Minimum MAPQ for read filtering [%d]\n" % plp_conf.MIN_MAPQ
    s += "  --countORPHAN     If use, do not skip anomalous read pairs.\n"
    s += "\n"
    s += "Note:\n"
    s += "For file format details of '--cellAnno', '--cnvProfile', '--feature', and\n"
    s += "'--phasedSNP', please see the full manual at\n"
    s += "https://github.com/hxj5/scCNVSimulator/blob/master/docs/manual.rst\n"
    s += "\n"

    fp.write(s)


def pipeline_main(argv):
    """Command-Line interface.
    @param argv   A list of cmdline parameters [list]
    @return       0 if success, -1 otherwise [int]
    """
    func = "pipeline_main"

    conf = Config()
    plp_conf = PlpConfig()
    simu_conf = SimuConfig()

    if len(argv) <= 2:
        usage(stderr, plp_conf.defaults, simu_conf.defaults)
        sys.exit(1)

    stdout.write("[I::%s] start ...\n" % (func, ))

    opts, args = getopt.getopt(argv[2:], "-s:-O:-b:-P:-h-p:-D:", [
                    "sam=", "outdir=", "barcode=",
                    "cellAnno=", "refCellTypes=", 
                    "cnvProfile=", "feature=",
                    "phasedSNP=",
                    "help",

                    "nproc=", 
                    "cellTAG=", "UMItag=", 
                    "minCOUNT=", "minMAF=", 
                    "debug=",
                    
                    "inclFLAG=", "exclFLAG=", "minLEN=", "minMAPQ=", "countORPHAN"
                ])

    for op, val in opts:
        if len(op) > 2:
            op = op.lower()
        if op in   ("-s", "--sam"): conf.sam_fn = val
        elif op in ("-O", "--outdir"): conf.out_dir = val
        elif op in ("-b", "--barcode"): conf.barcode_fn = val
        elif op in (      "--cellanno"): conf.cell_anno_fn = val
        elif op in (      "--refcelltypes"): conf.ref_cell_types_str = val
        elif op in (      "--cnvprofile"): conf.cnv_profile_fn = val
        elif op in (      "--feature"): conf.feature_fn = val
        elif op in ("-P", "--phasedsnp"): conf.snp_fn = val
        elif op in ("-h", "--help"): usage(stderr, plp_conf.defaults, simu_conf.defaults); sys.exit(1)

        elif op in ("-p", "--nproc"): conf.nproc = int(val)
        elif op in ("--celltag"): conf.cell_tag = val
        elif op in ("--umitag"): conf.umi_tag = val
        elif op in ("--mincount"): conf.min_count = int(val)
        elif op in ("--minmaf"): conf.min_maf = float(val)
        elif op in ("-D", "--debug"): conf.debug = int(val)

        elif op in ("--inclflag"): conf.incl_flag = int(val)
        elif op in ("--exclflag"): conf.excl_flag = int(val)
        elif op in ("--minlen"): conf.min_len = int(val)
        elif op in ("--minmapq"): conf.min_mapq = float(val)
        elif op in ("--countorphan"): conf.no_orphan = False

        else:
            stderr.write("[E::%s] invalid option: '%s'.\n" % (func, op))
            return(-1)

    stdout.write("[I::%s] check parameters ...\n" % (func, ))

    prepare_args(conf)

    stdout.write("[I::%s] merge clone-specific CNV profile ...\n" % (func, ))

    merge_cnv_profile(conf.cnv_profile_fn, conf.merged_cnv_profile_fn, 
                    max_gap = 1, verbose = True)

    stdout.write("[I::%s] pileup allele-specific UMIs ...\n" % (func, ))

    plp_argv = conf2plp_argv(conf)
    ret = pileup_main(plp_argv, plp_conf)
    if ret < 0:
        stderr.write("[E::%s] pileup failed; errcode %d.\n" % (func, ret))
        return(ret)

    stdout.write("[I::%s] simulate CNVs ...\n" % (func, ))

    simu_argv = conf2simu_argv(conf, plp_conf)
    ret = simu_main(simu_argv, simu_conf)
    if ret < 0:
        stderr.write("[E::%s] simulation failed; errcode %d.\n" % (func, ret))
        return(ret)

    stdout.write("[I::%s] All Done!\n" % (func, ))
    return(0)
    

COMMAND = "pipeline"


if __name__ == "__main__":
    pipeline_main(sys.argv)
