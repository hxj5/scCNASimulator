# utils.py - utils


import getopt
import numpy as np
import os
import pysam
import sys

from .region import format_chrom


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

