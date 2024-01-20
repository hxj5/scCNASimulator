# utils.py - utils


import numpy as np

from .blib.region import format_chrom, format_start, format_end
from .blib.zfile import zopen


# @param cells A list.
def get_clones(cells, frac0, seed = 123):
    n_cells = len(cells)
    n_c0 = int(n_cells * frac0)

    np.random.seed(seed)
    clone0 = set(np.random.choice(cells, size = n_c0, replace = False))
    clone1 = set(cells).difference(clone0)

    return (clone0, clone1)


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

