# utils.py - utils


import os


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


# @return {<cell>:{
#     <region>:[
#         set(),     # umi set for allele 0
#         set(),     # umi set for allele 1
#     ]
# }}
def load_pileup_umi(fn_list, gene_dict):
    func = "load_pileup_umi"

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

