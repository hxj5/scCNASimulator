# utils.py - utils


import functools
import getopt
import numpy as np
import os
import pysam
import sys

from .lib.region import format_chrom
from .lib.zfile import zopen
from .cnv import CloneCNVProfile


def assert_e(path):
    if path is None or not os.path.exists(path):
        raise OSError


def assert_n(var):
    if var is None or not var:
        raise ValueError


def assert_notnone(var):
    if var is None:
        raise ValueError


### CNV profile

def load_cnv_profile(fn, sep = "\t", verbose = False):
    """Load CNV profiles from file.
    @param fn       Path to file [str]
    @param verbose  If print detailed log info [bool]
    @return         A CloneCNVProfile object if success, None otherwise.
    @note           The first 6 columns of the file should be
                    clone_id, chrom, start, end (both 1-based, inclusive), cn_ale0, cn_ale1.
    """
    func = "load_cnv_profile"
    fp = zopen(fn, "rt")
    dat = CloneCNVProfile()
    nl = 0
    if verbose:
        sys.stderr.write("[I::%s] start to load CNV profile from file '%s' ...\n" % (func, fn))
    for line in fp:
        nl += 1
        items = line.strip().split(sep)
        if len(items) < 6:
            if verbose:
                sys.stderr.write("[E::%s] too few columns of line %d.\n" % (func, nl))
            return(None)
        clone_id, chrom, start, end, cn_ale0, cn_ale1 = items[:6]
        start, end = int(start), int(end)
        cn_ale0, cn_ale1 = int(cn_ale0), int(cn_ale1)
        dat.add_cnv(chrom, start, end, cn_ale0, cn_ale1, clone_id)
    fp.close()
    return(dat)


def __cmp_two_intervals(x1, x2):
    s1, e1 = x1[:2]
    s2, e2 = x2[:2]
    if s1 == s2:
        if e1 == e2:
            return(0)
        else:
            return e1 - e2
    else:
        return s1 - s2


def merge_cnv_profile(in_fn, out_fn, max_gap = 1):
    func = "merge_cnv_profile"

    # load data
    fp = zopen(in_fn, "rt")
    dat = {}
    nl = 0
    if verbose:
        sys.stderr.write("[I::%s] start to merge CNV profile from file '%s' ...\n" % (func, in_fn))
    for line in fp:
        nl += 1
        items = line.strip().split(sep)
        if len(items) < 6:
            sys.stderr.write("[E::%s] too few columns of line %d.\n" % (func, nl))
            return(None)
        clone_id, chrom, start, end, cn_ale0, cn_ale1 = items[:6]
        start, end = int(start), int(end)
        cn_ale0, cn_ale1 = int(cn_ale0), int(cn_ale1)
        if clone_id not in dat:
            dat[clone_id] = {}
        if chrom not in dat[clone_id]:
            dat[clone_id][chrom] = {}
        ale_key = "%d_%d" % (cn_ale0, cn_ale1)
        if ale_key not in dat[clone_id][chrom]:
            dat[clone_id][chrom][ale_key] = []
        dat[clone_id][chrom][ale_key].append((start, end))
    fp.close()

    # merge (clone-specific) adjacent CNVs
    for clone_id, cl_dat in dat.items():
        for chrom, ch_dat in cl_dat.items():
            for ale_key in ch_dat.keys():
                iv_list = sorted(ch_dat[ale_key], key = functools.cmp_to_key(__cmp_two_intervals))
                s1, e1 = iv_list[0]
                new_list = []
                for s2, e2 in iv_list[1:]:
                    if s2 <= e1 + max_gap:    # overlap adjacent region
                        e1 = max(e1, e2)
                    else:                     # otherwise
                        new_list.append((s1, e1))
                        s1, e1 = s2, e2
                new_list.append((s1, e1))
                ch_dat[ale_key] = new_list

    # check whether there are (strictly) overlapping regions with distinct profiles
    for clone_id, cl_dat in dat.items():
        for chrom, ch_dat in cl_dat.items():
            iv_list = []
            for ale_key in ch_dat.keys():
                cn_ale0, cn_ale1 = [int(x) for x in ale_key.split("_")]
                iv_list.extend([(s, e, cn_ale0, cn_ale1) for s, e in ch_dat[ale_key]])
            iv_list = sorted(iv_list, key = functools.cmp_to_key(__cmp_two_intervals))
            s1, e1 = iv_list[0][:2]
            for iv in in iv_list[1:]:
                s2, e2 = iv[:2]
                if s2 <= e1:    # overlap adjacent region
                    sys.stderr.write("[E::%s] distinct CNV profiles '%s', (%d, %d) and (%d, %d).\n" % 
                        (func, chrom, s1, e1, s2, e2))
                return(None)
            cl_dat[chrom] = iv_list

    # save profile
    fp = open(out_fn, "w")
    for clone_id in sorted(dat.keys()):
        cl_dat = dat[clone_id]
        for chrom in sorted(cl_dat.keys()):
            ch_dat = cl_dat[chrom]
            for s, e, cn_ale0, cn_ale1 in ch_dat:
                fp.write("\t".join([clone_id, chrom, str(s), str(e), str(cn_ale0), str(cn_ale1)]) + "\n")
    fp.close()
                        

def save_cnv_profile(fn):
    pass


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

