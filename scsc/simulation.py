# simulation.py - CNV simulation


import getopt
import numpy as np
import os
import pysam
import sys

from .app import APP, VERSION
from .blib.region import format_chrom
from .simu.cnv import load_cnv_profile, merge_cnv_profile
from .simu.config import Config
from .simu.utils import load_cell_anno
from .simu.utils import assert_e, assert_n


def prepare_args(conf):
    assert_n(conf.sid)

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


# in_sam: input sam/bam;
# out_sam: output sam/bam;
# chrom: target chrom;
# clone0: set of cell barcodes whose allele 0 would be removed.
# clone1: set of cell barcodes whose allele 1 would be removed.
# ale_umi:
# prob0: probability to discard ambiguous UMI in clone 0.
# prob1: probability to discard ambiguous UMI in clone 1.
def simu_cnv(in_sam, out_sam,
    chrom, cell_tag, umi_tag, 
    clone0, clone1, 
    ale_umi, 
    prob0 = 0.5, prob1 = 0.5, seed = 123):

    func = "simu_cnv"

    chrom = format_chrom(chrom)

    np.random.seed(seed)
    
    n_umi_c0_a0 = n_umi_c0_a1 = n_umi_c0_amb = 0
    n_umi_c1_a0 = n_umi_c1_a1 = n_umi_c1_amb = 0
    umi_del = {}

    for read in in_sam.fetch():
        read_chrom = read.reference_name
        if not read_chrom:
            out_sam.write(read)
            continue
        read_chrom = format_chrom(read_chrom)
        if read_chrom != chrom:    # UMI is not on target chrom.
            out_sam.write(read)
            continue

        if not read.has_tag(cell_tag):
            out_sam.write(read)
            continue
        cb = read.get_tag(cell_tag)

        if not read.has_tag(umi_tag):
            out_sam.write(read)
            continue
        umi = read.get_tag(umi_tag)

        if not cb or cb not in ale_umi:
            out_sam.write(read)
            continue
        if cb not in clone0 and cb not in clone1:
            out_sam.write(read)
            continue

        clone_id = 0 if cb in clone0 else 1   # cell is in target clone.
        if cb not in umi_del:
            umi_del[cb] = {
                "clone_id": clone_id,
                "a0": set(),
                "a1": set(),
                "amb": set()
            }

        cell_umi = ale_umi[cb]
        if umi in cell_umi[0]:    # UMI is on allele 0.
            if clone_id == 0:
                n_umi_c0_a0 += 1
                umi_del[cb]["a0"].add(umi)
            else:
                n_umi_c1_a0 += 1
                out_sam.write(read)
        elif umi in cell_umi[1]:
            if clone_id == 0:
                n_umi_c0_a1 += 1
                out_sam.write(read)
            else:
                n_umi_c1_a1 += 1
                umi_del[cb]["a1"].add(umi)
        else:                     # ambiguous UMI
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


def usage(fp = sys.stderr):
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

    clone0, clone1 = get_clones(target_cells, conf.clone_frac0, conf.seed)

    for idx, clone in zip(range(2), [clone0, clone1]):
        out_fn = os.path.join(conf.out_dir, "clone%d.barcodes.tsv" % idx)
        sort_clone = sorted(list(clone))
        s = "\n".join(sort_clone) + "\n"
        with open(out_fn, "w") as fp:
            fp.write(s)

    print("[I::%s] #cells_clone0=%d; #cells_clone1=%d." % (func,
            len(clone0), len(clone1)))

    out_fn = os.path.join(conf.out_dir, "cell_anno.with_clone.2column.tsv")
    sort_cells = sorted(list(cell_anno.keys()))
    s = ""
    for cell in sort_cells:
        cell_type = cell_anno[cell]
        clone_str = ""
        if cell in clone0:
            clone_str = "_clone0"
        elif cell in clone1:
            clone_str = "_clone1"
        s += '%s\t"%s"\n' % (cell, cell_type + clone_str)
    with open(out_fn, "w") as fp:
        fp.write(s)


    ###
    print("[I::%s] load gene annotation." % func)

    gene_dict_all = load_gene_anno(conf.gene_anno_fn)
    gene_dict = {}
    for gene in gene_dict_all.keys():
        chrom = format_chrom(gene_dict_all[gene]["chrom"])
        if chrom == conf.chrom:
            gene_dict[gene] = gene_dict_all[gene]

    print("[I::%s] #all_gene=%d; #chr%s_gene=%d.\n" % (func, 
            len(gene_dict_all), conf.chrom, len(gene_dict)))

    ###
    print("[I::%s] load gene-specific UMIs generated by xcltk-plp-umi." % func)

    umi_fn_list = []
    for fn in os.listdir(conf.umi_dir):
        if fn.startswith("gene_umi") and fn.endswith(".tsv"):
            umi_fn_list.append(os.path.join(conf.umi_dir, fn))

    gene_umi = load_xcltk_umi(umi_fn_list, gene_dict)

    ###
    print("[I::%s] load global phase result." % func)

    phase = load_global_phase(conf.global_phase_fn)

    print("[I::%s] #genes having phase information: %d." % (func, len(phase)))
    print("[I::%s] #genes that need to be flipped: %d." % (func,
        np.sum([1 if flip else 0 for gene, flip in phase.items()])))

    ###
    print("[I::%s] get allele-specific UMIs." % func)

    ale_umi, ale_umi_cnt = get_allele_umi(gene_umi, phase)

    ###
    print("[I::%s] simulate copy loss in each of the two clones." % func)

    in_sam = pysam.AlignmentFile(conf.bam_fn, "rb")
    out_sam_fn = os.path.join(conf.out_dir, "%s.out.bam" % conf.sid)
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
        usage(sys.stderr)
        sys.exit(1)

    conf = Config()
    opts, args = getopt.getopt(argv[2:], "", [
        "sid=", 
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
        if op in   ("--sid"): conf.sid = val
        elif op in ("--sam"): conf.bam_fn = val
        elif op in ("--outdir"): conf.out_dir = val
        elif op in ("--cellanno"): conf.cell_anno_fn = val
        elif op in ("--cnvprofile"): conf.cnv_profile_fn = val
        elif op in ("--umidir"): conf.umi_dir = val
        elif op in ("--celltag"): conf.cell_tag = val
        elif op in ("--umitag"): conf.umi_tag = val
        elif op in ("--version"): sys.stderr.write("%s\n" % VERSION); sys.exit(1)
        elif op in ("--help"): usage(); sys.exit(1)
        else:
            sys.stderr.write("[E::%s] invalid option: '%s'.\n" % (func, op))
            return(-1)    

    ret = simu_core(argv, conf)
    return(ret)


COMMAND = "simu"
CONF_SEED = 123


if __name__ == "__main__":
    simu_main(sys.argv)

