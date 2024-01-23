# simulation.py - CNV simulation


import getopt
import numpy as np
import os
import pysam
import sys
from sys import stdout, stderr
import time

from .app import APP, VERSION
from .blib.region import format_chrom
from .simu.allele import load_allele_umi
from .simu.cnv import load_cnv_profile, merge_cnv_profile
from .simu.config import Config
from .simu.utils import load_cell_anno
from .utils import assert_e, assert_n


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


# Note,
# we create new UMI barcodes by simply adding distinct suffix to the original
# UMI barcode. Alternatively, we can iterate the BAM file to get all unique UMI
# barcodes first, and then create new ones. However, the latter strategy is
# inefficient and not easy to implement.

def __write_read(read, sam, umi, umi_tag, qname = None, idx = 0, umi_suffix_len = 4):
    if idx < 0 or idx >= 2 ^ umi_suffix_len:
        raise ValueError

    if qname is None:
        read.query_name += "_" + str(idx)
    else:
        read.query_name = qname + "_" + str(idx)

    if not umi:
        sam.write(read)
        return

    x = ["ACGT"[(idx >> (2 * i)) & 3] for i in range(umi_suffix_len)]
    x.reverse()
    umi_suffix = "".join(x)

    umi += umi_suffix
    read.set_tag(umi_tag, umi)
    sam.write(read)


def simu_cnv(
    in_sam, out_sam,
    cell_anno, cnv_profile, allele_umi,
    cell_tag, umi_tag
):
    func = "simu_cnv"

    clone = None
    flag = -1
    cn = None      # copy number
    for read in in_sam.fetch():
        if not read.has_tag(umi_tag):
            __write_read(read, out_sam, None, umi_tag)
            continue
        umi = read.get_tag(umi_tag)

        if not read.has_tag(cell_tag):
            __write_read(read, out_sam, umi, umi_tag)
            continue
        cell = read.get_tag(cell_tag)

        if not read.query_name:
            __write_read(read, out_sam, umi, umi_tag)
            continue
        qname = read.query_name

        if cell not in cell_anno:    # cell not in CNV clone.
            __write_read(read, out_sam, umi, umi_tag)
            continue
        clone = cell_anno[cell]

        # check whether the read has allele information.

        res = allele_umi.query(cell, umi)

        if res is None:    # ambiguous UMIs

            # check whether read is in CNV region
            # for copy gain and copy loss, ambiguous UMIs may still be forked 
            # or deleted.

            # Note:
            # Currently fork or deletion in read level is not perfect
            # as it should be done in UMI level.
            # For example, one UMI may overlap with the copy gain regions, 
            # hence all its reads should be forked, while in current scheme,
            # reads of the UMI that do not overlap CNV regions will not be
            # forked.
            # This should have little influence on the CNV simulation, as the
            # number of these reads should be very small considering the length
            # of CNV regions and UMIs.
        
            chrom = read.reference_name
            if not chrom:
                __write_read(read, out_sam, umi, umi_tag)
                continue
            positions = read.get_reference_positions()        # 0-based
            if not positions:
                __write_read(read, out_sam, umi, umi_tag)
                continue
            start, end = positions[0] + 1, positions[-1] + 1  # 1-based

            ret, cn_ale = cnv_profile.fetch(chrom, start, end + 1, clone)
            if ret < 0:
                raise ValueError
            if ret != 0:           # not in CNV region or overlap multiple CNV regions with distinct CNV profiles.
                __write_read(read, out_sam, umi, umi_tag)
                continue
            cn0, cn1 = cn_ale[:2]

            rand_f = np.random.rand()
            if rand_f < 0.5:
                cn = cn0
            else:
                cn = cn1

            if cn <= 0:
                continue
            elif cn == 1:
                __write_read(read, out_sam, umi, umi_tag)
                continue
            else:
                # update the QNAME and UMI of the forked reads.
                for i in range(cn):
                    __write_read(read, out_sam, umi, umi_tag, qname, i)

        else:
            allele, reg_id_list = res[:2]
            if len(reg_id_list) != 1:
                __write_read(read, out_sam, umi, umi_tag)
                continue
            reg_id = reg_id_list[0]

            ret, cn_ale = cnv_profile.query(reg_id, clone)
            if ret < 0:
                raise ValueError
            if ret != 0:
                __write_read(read, out_sam, umi, umi_tag)
                continue
            cn = cn_ale[allele]     # copy number

            if cn == 0:
                continue
            elif cn == 1:
                __write_read(read, out_sam, umi, umi_tag)
                continue
            else:
                for i in range(cn):
                    __write_read(read, out_sam, umi, umi_tag, qname, i)


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
        simu_cnv(
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

