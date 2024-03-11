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
from .simu.baf import BAFIO, BAFCellReg
from .simu.core import simu_cnv
from .simu.allele import load_allele_umi
from .simu.cnv import load_cnv_profile, merge_cnv_profile
from .simu.config import Config
from .simu.utils import load_cell_anno, save_cell_anno, \
                        load_features, save_features


def prepare_args(conf):
    assert_e(conf.sam_fn)

    if not os.path.exists(conf.out_dir):
        os.mkdir(conf.out_dir)

    assert_e(conf.cell_anno_fn)
    assert_n(conf.ref_cell_types_str)
    conf.ref_cell_types = [s.strip().strip('"') for s in \
        conf.ref_cell_types_str.split(",")]

    assert_e(conf.cnv_profile_fn)
    assert_e(conf.feature_fn)
    assert_e(conf.baf_dir)
    assert_e(conf.umi_dir)

    assert_n(conf.cell_tag)
    assert_n(conf.umi_tag)

    conf.out_sam_fn = os.path.join(conf.out_dir, "out.bam")
    conf.out_cnv_fn = os.path.join(conf.out_dir, "merged_cnv_profile.tsv")
    conf.out_feature_fn = os.path.join(conf.out_dir, "features.tsv")
    conf.out_cell_anno_fn = os.path.join(conf.out_dir, "cell_anno.tsv")
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

        # load allele-specific UMIs.
        stdout.write("[I::%s] load allele-specific UMIs.\n" % func)
        fn_list = []
        for fn in os.listdir(conf.umi_dir):
            if fn.startswith(conf.umi_fn_prefix) and fn.endswith(conf.umi_fn_suffix):
                fn_list.append(os.path.join(conf.umi_dir, fn))
        allele_umi = load_allele_umi(fn_list, verbose = True)
    
        # load clone annotation.
        stdout.write("[I::%s] load clone annotation.\n" % func)
        cell_anno = load_cell_anno(conf.cell_anno_fn)
        save_cell_anno(cell_anno, conf.out_cell_anno_fn)

        # load cell-region BAF matrix
        stdout.write("[I::%s] load cell-region BAF matrix.\n" % func)
        baf_io = BAFIO(conf.baf_dir, conf.baf_fn_prefix)
        baf_adata = baf_io.load_data()
        baf_adata = baf_adata.transpose()

        if conf.debug:
            stderr.write("[D::%s] baf_adata is:\n" % func)
            stderr.write("%s\n" % str(baf_adata))
        
        # calc cell-region BAF (allelic imbalance information)
        stdout.write("[I::%s] calc cell-region BAF (allelic imbalance).\n" % func)
        cellreg_baf = BAFCellReg(baf_adata, cell_anno,   \
            cell_type_key = "cell_type",
            ref_cell_types = conf.ref_cell_types,
            theo = 0.5)
        
        if conf.debug:
            stderr.write("[D::%s] cellreg_baf is:\n" % func)
            stderr.write("%s\n" % str(cellreg_baf.baf))

        # merge CNV profile.
        stdout.write("[I::%s] merge CNV profile.\n" % func)
        merge_cnv_profile(conf.cnv_profile_fn, conf.out_cnv_fn, 
                                    max_gap = 1, verbose = True)
        cnv_profile = load_cnv_profile(conf.out_cnv_fn, sep = "\t",
                                    verbose = True)
        
        # load features
        stdout.write("[I::%s] load features.\n" % func)
        features = load_features(conf.feature_fn)
        save_features(features, conf.out_feature_fn)
    
        # simulate copy number variations.
        stdout.write("[I::%s] simulate copy number variations.\n" % func)
        in_sam = pysam.AlignmentFile(conf.sam_fn, "rb")
        out_sam = pysam.AlignmentFile(conf.out_sam_fn, "wb", template = in_sam)
        umi_stat = simu_cnv(
            in_sam = in_sam,
            out_sam = out_sam,
            allele_umi = allele_umi,
            cell_anno = cell_anno,
            cellreg_baf = cellreg_baf,
            cnv_profile = cnv_profile,
            features = features,
            cell_tag = conf.cell_tag,
            umi_tag = conf.umi_tag,
            debug = conf.debug
        )
        in_sam.close()
        out_sam.close()

        # index the BAM file
        stdout.write("[I::%s] index the BAM file.\n" % func)
        pysam.index(conf.out_sam_fn)

        # output UMI stat
        stdout.write("[I::%s] output UMI statistics.\n" % func)
        fp = open(conf.out_umi_stat_fn, "w")
        s_keys = ("cell", "clone", "region",
                    "A0", "A0_amp", "A0_del",
                    "A1", "A1_amp", "A1_del", 
                    "AMB", "AMB_amp", "AMB_del", 
                    "invalid", "unknown")
        s = "\t".join(s_keys) + "\n"
        fp.write(s)
        for cell, c_dat in umi_stat.items():
            clone = cell_anno[cell]
            for region, r_dat in c_dat.items():
                s = "\t".join([
                    cell, clone, region,
                    str(r_dat["A0"]), str(r_dat["A0_amp"]), str(r_dat["A0_del"]),
                    str(r_dat["A1"]), str(r_dat["A1_amp"]), str(r_dat["A1_del"]),
                    str(r_dat["AMB"]), str(r_dat["AMB_amp"]), str(r_dat["AMB_del"]),
                    str(r_dat["invalid"]), str(r_dat["unknown"])
                ]) + "\n"
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
            

def usage(fp = stderr, conf = None):
    s =  "\n" 
    s += "Version: %s\n" % (VERSION, )
    s += "Usage:   %s %s <options>\n" % (APP, COMMAND)
    s += "\n" 
    s += "Options:\n"
    s += "  --sam FILE             Indexed BAM/SAM/CRAM file.\n"
    s += "  --outdir DIR           Output dir.\n"
    s += "  --cellAnno FILE        Cell annotation file, 2 columns.\n"
    s += "  --refCellTypes STR     Reference cell types, comma separated.\n"
    s += "  --cnvProfile FILE      CNV profile file, 7 columns.\n"
    s += "  --feature FILE         Feature annotation file, typically for genes; 4 columns.\n"
    s += "  --BAFdir DIR           Dir storing cell-region BAF matrix.\n"
    s += "  --UMIdir DIR           Dir storing region-specific UMI files.\n"
    s += "  --version              Print version and exit.\n"
    s += "  --help                 Print this message and exit.\n"
    s += "\n"
    s += "Optional arguments:\n"
    s += "  --cellTAG STR          Cell barcode tag [%s]\n" % conf.CELL_TAG
    s += "  --UMItag STR           UMI tag [%s]\n" % conf.UMI_TAG
    s += "  --debug INT            Debug level, only for developer [%d]\n" % conf.DEBUG
    s += "\n"

    fp.write(s)


def simu_main(argv, conf = None):
    func = "simu_main"

    if conf is None:
        conf = Config()

    if len(argv) <= 2:
        usage(stderr, conf.defaults)
        sys.exit(1)

    opts, args = getopt.getopt(argv[2:], "", [
        "sam=",
        "outdir=",
        "cellAnno=", "refCellTypes=",
        "cnvProfile=", "feature=",
        "BAFdir=", "UMIdir=",
        "version", "help",
        
        "cellTAG=", "UMItag=",
        "debug="
    ])

    for op, val in opts:
        if len(op) > 2:
            op = op.lower()
        if op in ("--sam"): conf.sam_fn = val
        elif op in ("--outdir"): conf.out_dir = val
        elif op in ("--cellanno"): conf.cell_anno_fn = val
        elif op in ("--refcelltypes"): conf.ref_cell_types_str = val
        elif op in ("--cnvprofile"): conf.cnv_profile_fn = val
        elif op in (" --feature"): conf.feature_fn = val
        elif op in ("--bafdir"): conf.baf_dir = val
        elif op in ("--umidir"): conf.umi_dir = val
        elif op in ("--version"): stderr.write("%s\n" % VERSION); sys.exit(1)
        elif op in ("--help"): usage(); sys.exit(1)

        elif op in ("--celltag"): conf.cell_tag = val
        elif op in ("--umitag"): conf.umi_tag = val
        elif op in ("--debug"): conf.debug = int(val)
        else:
            stderr.write("[E::%s] invalid option: '%s'.\n" % (func, op))
            return(-1)

    ret = simu_core(argv, conf)
    return(ret)


COMMAND = "simu"
CONF_SEED = 123


if __name__ == "__main__":
    simu_main(sys.argv)

