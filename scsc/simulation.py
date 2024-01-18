# simulation.py - CNV simulation


import getopt
import numpy as np
import os
import pysam
import sys
from sys import stdout, stderr

from .app import APP, VERSION
from .blib.region import format_chrom
from .simu.allele import load_allele_umi
from .simu.cnv import load_cnv_profile, merge_cnv_profile
from .simu.config import Config
from .simu.utils import load_cell_anno, assert_e, assert_n


def prepare_args(conf):
    assert_e(conf.bam_fn)

    if not os.path.exists(conf.out_dir):
        os.mkdir(conf.out_dir)

    assert_e(conf.cell_anno_fn)
    conf.cell_anno = load_cell_anno(conf.cell_anno_fn)

    assert_e(conf.cnv_profile_fn)
    merged_cnv_profile_fn = os.path.join(conf.out_dir, "merged.cnv_profile.tsv")
    merge_cnv_profile(conf.cnv_profile_fn, merged_cnv_profile_fn, max_gap = 1)
    conf.cnv_profile = load_cnv_profile(merged_cnv_profile_fn, sep = "\t",
                                        verbose = True)

    assert_e(conf.umi_dir)

    assert_n(conf.cell_tag)
    assert_n(conf.umi_tag)


def __write_read(read, sam, umi, umi_tag):
    read.query_name += "_0"
    if not umi:
        sam.write(read)
    umi += "AAAA"
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
    for read in in_sam.fetch():
        if not read.has_tag(umi_tag):
            __write_read(read, out_sam, umi, umi_tag)
            continue
        umi = read.get_tag(umi_tag)

        if not read.has_tag(cell_tag):
            __write_read(read, out_sam, umi, umi_tag)
            continue
        cell = read.get_tag(cell_tag)

        if cell not in cell_anno:    # cell not in CNV clone.
            __write_read(read, out_sam, umi, umi_tag)
            continue
        clone = cell_anno[cell]

        # check whether read is in CNV region
        # even if it is ambiguous UMI, for copy gain and copy loss, it
        # may still be forked or deleted.
        
        chrom = read.reference_name
        if not chrom:
            __write_read(read, out_sam, umi, umi_tag)
            continue
        positions = read.get_reference_positions()        # 0-based
        if not positions:
            __write_read(read, out_sam, umi, umi_tag)
            continue
        start, end = positions[0] + 1, positions[-1] + 2

        ###
        res = allele_umi.query(cell, umi)
        if res is None:              # ambiguous UMIs
            __write_read(read, out_sam, umi, umi_tag)
            continue
        allele, reg_id_list = res[:2]
        if len(reg_id_set) != 1:
            flag = 2
        reg_id = reg_id_list[0]


        ret, cn_ale = cnv_profile.query(reg_id, clone)
        if ret < 0:
            raise ValueError
        if ret != 0:
            flag = 3
        cn = cn_ale[allele]     # copy number

        if cn == 0:
            continue
        elif cn == 1:
            __write_read(read, out_sam, umi, umi_tag)
            continue
        else:
            

            rand_f = np.random.rand()
            if clone_id == 0:
                n_umi_c0_amb += 1
                if rand_f < prob0:
                    umi_del[cb]["amb"].add(umi)
                else:
                    out_sam.write(read)
            else:
                n_umi_c1_amb += 1
                if rand_f < prob1:
                    umi_del[cb]["amb"].add(umi)
                else:
                    out_sam.write(read)

    print("[I::%s] #umi_c0_a0=%d; #umi_c0_a1=%d; #umi_c0_amb=%d." % (
        func, n_umi_c0_a0, n_umi_c0_a1, n_umi_c0_amb))
    print("[I::%s] #umi_c1_a0=%d; #umi_c1_a1=%d; #umi_c1_amb=%d." % (
        func, n_umi_c1_a0, n_umi_c1_a1, n_umi_c1_amb))

    n_umi_del = {}
    for cell in umi_del.keys():
        n_umi_del[cell] = {
            "clone_id": umi_del[cell]["clone_id"],
            "a0": len(umi_del[cell]["a0"]),
            "a1": len(umi_del[cell]["a1"]),
            "amb": len(umi_del[cell]["amb"])
        }

    return (n_umi_del)


def usage(fp = stderr):
    s =  "\n" 
    s += "Version: %s\n" % (VERSION, )
    s += "Usage: %s %s <options>\n" % (APP, COMMAND)
    s += "\n" 
    s += "Options:\n"
    s += "  --sid STR              Sample ID.\n"
    s += "  --sam FILE             Indexed BAM/SAM/CRAM file.\n"
    s += "  --outdir DIR           Output dir.\n"
    s += "  --cellAnno FILE        Cell annotation file, 2 columns.\n"
    s ++ "  --cnvProfile FILE      CNV profile file, 6 columns.\n"
    s += "  --UMIdir DIR           Dir storing gene-specific UMI files.\n"
    s += "  --cellTAG STR          Cell barcode tag.\n"
    s += "  --UMItag STR           UMI tag.\n"
    s += "  --version              Print version and exit.\n"
    s += "  --help                 Print this message and exit.\n"
    s += "\n"

    fp.write(s)


def simu_core(argv, conf):
    func = "simu_core"
    ret = -1

    prepare_args(conf)

    ###
    stdout.write("[I::%s] load allele-specific UMIs.\n" % func)

    fn_list = []
    for fn in os.listdir(conf.umi_dir):
        if fn.startswith("allele_umi") and fn.endswith(".tsv"):
            fn_list.append(os.path.join(conf.umi_dir, fn))

    allele_umi = load_allele_umi(fn_list, verbose = True)

    ###
    print("[I::%s] simulate copy loss in each of the two clones." % func)

    in_sam = pysam.AlignmentFile(conf.bam_fn, "rb")
    out_sam_fn = os.path.join(conf.out_dir, "out.bam")
    out_sam = pysam.AlignmentFile(out_sam_fn, "wb", template = in_sam)

    n_umi_del = simu_loss(
        in_sam = in_sam,
        out_sam = out_sam,
        chrom = conf.chrom,
        cell_tag = conf.cell_tag,
        umi_tag = conf.umi_tag,
        clone0 = clone0,
        clone1 = clone1,
        ale_umi = ale_umi,
        prob0 = conf.clone_prob[0],
        prob1 = conf.clone_prob[1],
        seed = conf.seed
    )

    in_sam.close()
    out_sam.close()

    out_fn = os.path.join(conf.out_dir, "cell_by_allele_UMI_counts.tsv")
    s = "\t".join(["cell", "clone_id", "umi_a0", "umi_a1", "umi_shared",
            "umi_del_a0", "umi_del_a1", "umi_del_amb"]) + "\n"
    for cell in sorted(n_umi_del.keys()):
        flag = True if cell in ale_umi_cnt else False
        s += "\t".join([
            cell,
            str(n_umi_del[cell]["clone_id"]),
            str(ale_umi_cnt[cell]["a0"]) if flag else "NA",
            str(ale_umi_cnt[cell]["a1"]) if flag else "NA",
            str(ale_umi_cnt[cell]["shared"]) if flag else "NA",
            str(n_umi_del[cell]["a0"]),
            str(n_umi_del[cell]["a1"]),
            str(n_umi_del[cell]["amb"]),
        ]) + "\n"
    with open(out_fn, "w") as fp:
        fp.write(s)

    print("[I::%s] All Done!" % func)


def simu_main(argv):
    func = "simu_main"

    if len(argv) <= 2:
        usage(stderr)
        sys.exit(1)

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
        if op in ("--sam"): conf.bam_fn = val
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

