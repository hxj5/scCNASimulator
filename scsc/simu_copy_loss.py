# simu_loss.py - data simulation that two clones harbor copy loss of different alleles.


import getopt
import numpy as np
import os
import pysam
import sys


class Config:
    def __ini__(self):
        self.sid = None

        self.bam_fn = None
        self.chrom = None
        self.cell_tag = None
        self.umi_tag = None

        self.umi_dir = None

        self.global_phase_fn = None

        self.cell_anno_fn = None
        self.target_cell_types = None

        self.clone_frac0 = None
        self.clone_prob = None
        
        self.seed = None

        self.gene_anno_fn = None

        self.out_dir = None


    def check_args(self):
        assert_n(self.sid)

        assert_e(self.bam_fn)
        assert_n(self.chrom)
        self.chrom = format_chrom(self.chrom)
        assert_n(self.cell_tag)
        assert_n(self.umi_tag)

        assert_e(self.umi_dir)

        assert_e(self.global_phase_fn)

        assert_e(self.cell_anno_fn)
        assert_n(self.target_cell_types)
        self.target_cell_types = set(self.target_cell_types.split(";"))

        assert_notnone(self.clone_frac0)
        assert_n(self.clone_prob)
        clone_prob = self.clone_prob.split(",")
        if len(clone_prob) != 2:
            raise ValueError
        self.clone_prob = [float(x) for x in clone_prob]
        
        assert_notnone(self.seed)

        assert_e(self.gene_anno_fn)

        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)

        
def assert_e(path):
    if path is None or not os.path.exists(path):
        raise OSError


def assert_n(var):
    if var is None or not var:
        raise ValueError


def assert_notnone(var):
    if var is None:
        raise ValueError


def load_cell_anno(fn):
    cell_anno = {}
    with open(fn, "r") as fp:
        for line in fp:
            items = line.strip().split("\t")
            cell, cell_type = items[:2]
            cell, cell_type = cell.strip('"'), cell_type.strip('"')
            cell_anno[cell] = cell_type
    return(cell_anno)


# @param cells A list.
def get_clones(cells, frac0, seed = 123):
    n_cells = len(cells)
    n_c0 = int(n_cells * frac0)

    np.random.seed(seed)
    clone0 = set(np.random.choice(cells, size = n_c0, replace = False))
    clone1 = set(cells).difference(clone0)

    return (clone0, clone1)


def format_chrom(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom


# e.g., xcltk/preprocess/data/annotate_genes_hg38_update_20230126.txt
# @return {<gene>:{
#     "chrom":<chrom>,
#     "start":<start>,
#     "end":<end>   
# }}
def load_gene_anno(fn):
    func = "load_gene_anno"

    gene_dict = {}
    with open(fn, "r") as fp:
        for line in fp:
            items = line.strip().split("\t")
            chrom, start, end, name = items[:4]
            if name in gene_dict: 
                print("[W::%s] '%s' is duplicate." % (func, name))
                continue
            gene_dict[name] = {
                "chrom": format_chrom(chrom),
                "start": start,
                "end": end,
            }

    return(gene_dict)


# @return {<cell>:{
#     <gene>:[
#         set(),     # umi set for allele 0
#         set(),     # umi set for allele 1
#     ]
# }}
def load_xcltk_umi(umi_fn_list, gene_dict):
    func = "load_xcltk_umi"

    gene_umi = {}
    
    for fn in umi_fn_list:
        if not os.path.exists(fn):
            print("[W::%s] '%s' does not exist in umi_dir." % (func, fn))
            raise OSError
        with open(fn, "r") as fp:
            for line in fp:
                items = line.strip().split("\t")
                cell, gene, umi, ale_idx = items[:4]
                if gene not in gene_dict:
                    continue
                ale_idx = int(ale_idx)
                if cell not in gene_umi:
                    gene_umi[cell] = {}
                if gene not in gene_umi[cell]:
                    gene_umi[cell][gene] = [set(), set()]
                gene_umi[cell][gene][ale_idx].add(umi)

    return(gene_umi)


# @return {<gene>: flip}
def load_global_phase(fn):
    func = "load_global_phase"

    phase = {}
    nl = 0
    with open(fn, "r") as fp:
        for line in fp:
            nl += 1
            if nl == 1:
                continue
            items = line.strip().split("\t")
            gene, flip = items[0], items[-1]
            flip = True if flip.strip().lower() == "true" else False
            if gene in phase:
                print("[W::%s] '%s' is duplicate." % (func, gene))
                continue
            phase[gene] = flip

    return(phase)


# @return {<cell>:[
#     set(<umi>),     # UMIs of allele 0
#     set(<umi>),     # UMIs of allele 1
# ]}
def get_allele_umi(gene_umi, phase):
    func = "get_allele_umi"

    ale_umi = {}
    ale_umi_cnt = {}
    for cell in gene_umi.keys():
        ale_umi[cell] = [set(), set()]

        for gene in gene_umi[cell].keys():
            if gene not in phase:
                print("[W::%s] '%s' has no phasing information." % (func, gene))
                continue
            flip = phase[gene]
            if flip:
                ale_umi[cell][0].update(gene_umi[cell][gene][1])
                ale_umi[cell][1].update(gene_umi[cell][gene][0])
            else:
                ale_umi[cell][0].update(gene_umi[cell][gene][0])
                ale_umi[cell][1].update(gene_umi[cell][gene][1])

        n_umi_a0 = len(ale_umi[cell][0])
        n_umi_a1 = len(ale_umi[cell][1])

        umi_shared = ale_umi[cell][0].intersection(ale_umi[cell][1])
        n_umi_shared = len(umi_shared)
        if n_umi_shared > 0:
            ale_umi[cell][0] = ale_umi[cell][0].difference(umi_shared)
            ale_umi[cell][1] = ale_umi[cell][1].difference(umi_shared)

        ale_umi_cnt[cell] = {
            "a0": n_umi_a0,
            "a1": n_umi_a1,
            "shared": n_umi_shared
        }

    return(ale_umi, ale_umi_cnt)


# in_sam: input sam/bam;
# out_sam: output sam/bam;
# chrom: target chrom;
# clone0: set of cell barcodes whose allele 0 would be removed.
# clone1: set of cell barcodes whose allele 1 would be removed.
# ale_umi:
# prob0: probability to discard ambiguous UMI in clone 0.
# prob1: probability to discard ambiguous UMI in clone 1.
def simu_loss(in_sam, out_sam,
    chrom, cell_tag, umi_tag, 
    clone0, clone1, 
    ale_umi, 
    prob0 = 0.5, prob1 = 0.5, seed = 123):

    func = "simu_loss"

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
    s += "Usage: %s <options>\n" % (APP, )  
    s += "\n" 
    s += "Options:\n"
    s += "  --sid STR              Sample ID.\n"
    s += "  --bam FILE             BAM file.\n"
    s += "  --chrom STR            Target chromosome.\n"
    s += "  --cellTAG STR          Cell barcode tag.\n"
    s += "  --UMItag STR           UMI tag.\n"
    s += "  --UMIdir DIR           Dir storing gene-specific UMI files.\n"
    s += "  --globalPhase FILE     File listing global phasing result.\n"
    s += "  --outdir DIR           Output dir.\n"
    s += "  --targetCellTypes STR  Cell types to be sampled, semicolon separated.\n"
    s += "  --cloneFrac0 FLOAT     Fraction of clone 0.\n"
    s += "  --cloneProb STR        Probability to discard ambiguous UMIs in two clones, comma separated.\n"
    s += "  --cellAnno FILE        Cell annotation file, 2 columns.\n"
    s += "  --geneAnno FILE        Gene annotation file, at least 4 columns.\n"
    s += "  --seed INT             Random seed [%d]\n" % CONF_SEED
    s += "  --version              Print version and exit.\n"
    s += "  --help                 Print this message and exit.\n"
    s += "\n"

    fp.write(s)


def main():
    func = "main"

    if len(sys.argv) <= 1:
        usage(sys.stderr)
        sys.exit(1)

    conf = Config()
    opts, args = getopt.getopt(sys.argv[1:], "", [
        "sid=", 
        "bam=", "chrom=", "cellTAG=", "UMItag=",
        "UMIdir=",
        "globalPhase=",
        "cellAnno=", "targetCellTypes=",
        "cloneFrac0=", "cloneProb=",
        "seed=", 
        "geneAnno=",
        "outdir=",
        "version", "help"
    ])

    for op, val in opts:
        if len(op) > 2:
            op = op.lower()
        if op in   ("--sid"): conf.sid = val
        elif op in ("--bam"): conf.bam_fn = val
        elif op in ("--chrom"): conf.chrom = val
        elif op in ("--celltag"): conf.cell_tag = val
        elif op in ("--umitag"): conf.umi_tag = val
        elif op in ("--umidir"): conf.umi_dir = val
        elif op in ("--globalphase"): conf.global_phase_fn = val
        elif op in ("--cellanno"): conf.cell_anno_fn = val
        elif op in ("--targetcelltypes"): conf.target_cell_types = val
        elif op in ("--clonefrac0"): conf.clone_frac0 = float(val)
        elif op in ("--cloneprob"): conf.clone_prob = val
        elif op in ("--seed"): conf.seed = int(val)
        elif op in ("--geneanno"): conf.gene_anno_fn = val
        elif op in ("--outdir"): conf.out_dir = val
        elif op in ("--version"): sys.stderr.write("%s\n" % VERSION); sys.exit(1)
        elif op in ("--help"): usage(); sys.exit(1)
        else:
            sys.stderr.write("[E::%s] invalid option: '%s'.\n" % (func, op))
            return(-1)    

    conf.check_args()

    ###
    print("[I::%s] load cell annotation." % func)
    
    cell_anno = load_cell_anno(conf.cell_anno_fn)
    target_cells = []
    for cell, cell_type in cell_anno.items():
        if cell_type in conf.target_cell_types:
            target_cells.append(cell)

    print("[I::%s] #cells_all=%d; #cells_target=%d." % (func,
            len(cell_anno), len(target_cells)))

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


APP = "simu_copy_loss.py"
VERSION = "0.0.1"

CONF_SEED = 123


if __name__ == "__main__":
    main()

