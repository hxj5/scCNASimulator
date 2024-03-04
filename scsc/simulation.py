# simulation.py - CNV simulation


import getopt
import numpy as np
import os
import pysam
import sys
from sys import stdout, stderr
import time

from .app import APP, VERSION
from .blib.w_assert import assert_e, assert_n
from .blib.region import format_chrom
from .simu.core import simu_cnv
from .simu.allele import load_allele_umi
from .simu.cnv import load_cnv_profile, merge_cnv_profile
from .simu.config import Config
from .simu.utils import load_cell_anno


def prepare_args(conf):
    assert_e(conf.sam_fn)

    if not os.path.exists(conf.out_dir):
        os.mkdir(conf.out_dir)

    assert_e(conf.cell_anno_fn)
    assert_e(conf.cnv_profile_fn)
    assert_e(conf.umi_dir)

    assert_n(conf.cell_tag)
    assert_n(conf.umi_tag)

    conf.merged_cnv_profile_fn = os.path.join(conf.out_dir, "merged.cnv_profile.tsv")
    conf.out_sam_fn = os.path.join(conf.out_dir, "out.bam")
    conf.out_umi_stat_fn = os.path.join(conf.out_dir, "cnv_umi_stat.tsv")


def simu_core(argv, conf):
    func = "simu_core"
    ret = -1
    cmdline = "[%s] unexpected command line" % func

    start_time = time.time()

    try:
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time))
        stdout.write("[I::%s] start time: %s.\n" % (func, time_str))

        cmdline = " ".join(argv)
        stdout.write("[I::%s] CMD: %s\n" % (func, cmdline))
    
        prepare_args(conf)
        conf.show(stderr)
    
        # load clone annotation.
        stdout.write("[I::%s] load clone annotation.\n" % func)
        cell_anno = load_cell_anno(conf.cell_anno_fn)
    
        # merge CNV profile.
        stdout.write("[I::%s] merge CNV profile.\n" % func)
        merge_cnv_profile(conf.cnv_profile_fn, conf.merged_cnv_profile_fn, 
                                    max_gap = 1, verbose = True)
        cnv_profile = load_cnv_profile(conf.merged_cnv_profile_fn, sep = "\t",
                                    verbose = True)
    
        # load allele-specific UMIs.
        stdout.write("[I::%s] load allele-specific UMIs.\n" % func)
        fn_list = []
        for fn in os.listdir(conf.umi_dir):
            if fn.startswith(conf.umi_fn_prefix) and fn.endswith(conf.umi_fn_suffix):
                fn_list.append(os.path.join(conf.umi_dir, fn))
        allele_umi = load_allele_umi(fn_list, verbose = True)
    
        # simulate copy number variations.
        stdout.write("[I::%s] simulate copy number variations.\n" % func)
        in_sam = pysam.AlignmentFile(conf.sam_fn, "rb")
        out_sam = pysam.AlignmentFile(conf.out_sam_fn, "wb", template = in_sam)
        umi_stat = simu_cnv(
            in_sam = in_sam,
            out_sam = out_sam,
            cell_anno = cell_anno,
            cnv_profile = cnv_profile,
            allele_umi = allele_umi,
            cell_tag = conf.cell_tag,
            umi_tag = conf.umi_tag,
        )
        in_sam.close()
        out_sam.close()

        # index the BAM file
        stdout.write("[I::%s] index the BAM file.\n" % func)
        pysam.index(conf.out_sam_fn)

        # output UMI stat
        fp = open(conf.out_umi_stat_fn, "w")
        s_keys = ("A0", "A0_amp", "A0_del", "A1", "A1_amp", "A1_del", 
                "AMB", "AMB_amp", "AMB_del", "invalid")
        s = "\t".join(s_keys) + "\n"
        fp.write(s)
        for cell, c_dat in umi_stat.items():
            s = "\t".join([str(c_dat[k]) for k in s_keys]) + "\n"
            fp.write(s)
        fp.close()

    except ValueError as e:
        stderr.write("[E::%s] '%s'\n" % (func, str(e)))
        stdout.write("[E::%s] Running program failed.\n" % func)
        stdout.write("[E::%s] Quiting ...\n" % func)
        ret = -1

    else:
        stdout.write("[I::%s] All Done!\n" % func)
        ret = 0

    finally:
        stdout.write("[I::%s] CMD: %s\n" % (func, cmdline))
        end_time = time.time()
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
        stdout.write("[I::%s] end time: %s\n" % (func, time_str))
        stdout.write("[I::%s] time spent: %.2fs\n" % (func, end_time - start_time))

    return(ret)
            

def usage(fp = stderr):
    s =  "\n" 
    s += "Version: %s\n" % (VERSION, )
    s += "Usage: %s %s <options>\n" % (APP, COMMAND)
    s += "\n" 
    s += "Options:\n"
    s += "  --sam FILE             Indexed BAM/SAM/CRAM file.\n"
    s += "  --outdir DIR           Output dir.\n"
    s += "  --cellAnno FILE        Cell annotation file, 2 columns.\n"
    s += "  --cnvProfile FILE      CNV profile file, 7 columns.\n"
    s += "  --UMIdir DIR           Dir storing gene-specific UMI files.\n"
    s += "  --cellTAG STR          Cell barcode tag.\n"
    s += "  --UMItag STR           UMI tag.\n"
    s += "  --version              Print version and exit.\n"
    s += "  --help                 Print this message and exit.\n"
    s += "\n"

    fp.write(s)


def simu_main(argv, conf = None):
    func = "simu_main"

    if len(argv) <= 2:
        usage(stderr)
        sys.exit(1)

    if conf is None:
        conf = Config()

    opts, args = getopt.getopt(argv[2:], "", [
        "sam=",
        "outdir=",
        "cellAnno=", "cnvProfile=",
        "UMIdir=",
        "cellTAG=", "UMItag=",
        "version", "help"
    ])

    for op, val in opts:
        if len(op) > 2:
            op = op.lower()
        if op in ("--sam"): conf.sam_fn = val
        elif op in ("--outdir"): conf.out_dir = val
        elif op in ("--cellanno"): conf.cell_anno_fn = val
        elif op in ("--cnvprofile"): conf.cnv_profile_fn = val
        elif op in ("--umidir"): conf.umi_dir = val
        elif op in ("--celltag"): conf.cell_tag = val
        elif op in ("--umitag"): conf.umi_tag = val
        elif op in ("--version"): stderr.write("%s\n" % VERSION); sys.exit(1)
        elif op in ("--help"): usage(); sys.exit(1)
        else:
            stderr.write("[E::%s] invalid option: '%s'.\n" % (func, op))
            return(-1)    

    ret = simu_core(argv, conf)
    return(ret)


COMMAND = "simu"
CONF_SEED = 123


if __name__ == "__main__":
    simu_main(sys.argv)

