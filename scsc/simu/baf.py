# baf.py - process BAF matrix

import anndata as ad
import numpy as np
import os
import pandas as pd
import scipy as sp
from scipy import io
from scipy import sparse


class BAFCellReg:
    """calibrated / posterior BAF in each cell type group."""

    def __init__(self, adata, cell_anno, cell_type_key = "cell_type", 
        ref_cell_types = None, theo = 0.5):
        self.adata = adata
        self.cell_anno = cell_anno
        self.cell_type_key = cell_type_key
        self.ref_cell_types = ref_cell_types
        self.theo = theo

        self.baf = self.__calc_baf()

    def __calc_baf(self):
        func = "BAFCellReg::__calc_baf"

        adata = self.adata
        cell_anno = self.cell_anno
        cell_key = "cell"
        group_key = self.cell_type_key
        theo = self.theo
        pseudo_counts = [5, 10]    # for AD, DP

        if group_key not in adata.obs.columns:
            if cell_key not in adata.obs.columns:
                raise KeyError("[E::%s] cell key '%s' not in adata." % (func, cell_key))
            adata.obs[group_key] = adata.obs[cell_key].map(cell_anno)
            
        grouped = adata.obs.groupby(group_key)

        # Note:
        # for now, the simulation is assumed to be performed on normal
        # cells or normal regions of tumor cells (e.g., chr3 on GX109),
        # so here we calculate the allelic imbalance in every cell type,
        # not using values of reference cells to replace the values of
        # tumor/CNV cells.

        adata.layers["BAF"] = np.full(adata.shape, theo)
        for group, idx in grouped.indices.items():
            reg_AD = adata.X[idx, :].sum(axis = 0) + pseudo_counts[0]
            reg_DP = adata.layers["DP"][idx, :].sum(axis = 0) + pseudo_counts[1]
            adata.layers["BAF"][idx, :] = (adata.X[idx, :] + reg_AD) /      \
                (adata.layers["DP"][idx, :] + reg_DP + 1e-8)

        res = {}
        n_cell, n_reg = adata.shape
        for i in range(n_cell):
            cell = adata.obs["cell"][i]
            if cell not in res:
                res[cell] = {}
            for j in range(n_reg):
                reg_id = adata.var["reg_id"][j]
                res[cell][reg_id] = adata.layers["BAF"][i, j]

        return(res)

    def query(self, cell, reg_id):
        if cell not in self.baf or reg_id not in self.baf[cell]:
            return(None)
        else:
            return(self.baf[cell][reg_id])
        

class BAFIO:
    def __init__(self, data_dir, file_prefix, is_gzip = True, is_genotype = False):
        self.data_dir = data_dir
        self.file_prefix = file_prefix
        self.is_gzip = is_gzip
        self.is_genotype = is_genotype

    def load_data(self):
        prefix = self.file_prefix

        regions = load_regions(os.path.join(self.data_dir, prefix + ".region.tsv"))
        samples = load_samples(os.path.join(self.data_dir, prefix + ".samples.tsv"))
        AD_mtx = load_matrix(os.path.join(self.data_dir, prefix + ".AD.mtx"))
        DP_mtx = load_matrix(os.path.join(self.data_dir, prefix + ".DP.mtx"))
        OTH_mtx = load_matrix(os.path.join(self.data_dir, prefix + ".OTH.mtx"))
        
        adata = ad.AnnData(
            X = AD_mtx, 
            obs = regions,
            var = samples)           
        adata.layers["DP"] = DP_mtx
        adata.layers["OTH"] = OTH_mtx

        self.adata = adata
        return(adata)

    def save_data(self, out_dir):
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        
        adata = self.adata
        prefix = self.file_prefix

        save_regions(adata.obs, 
            fn = os.path.join(out_dir, prefix + ".region.tsv"))
        save_samples(adata.var, 
            fn = os.path.join(out_dir, prefix + ".samples.tsv"))
        save_matrix(adata.X, os.path.join(out_dir, prefix + ".AD.mtx"))
        save_matrix(adata.layers["DP"], os.path.join(out_dir, prefix + ".DP.mtx"))
        save_matrix(adata.layers["OTH"], os.path.join(out_dir, prefix + ".OTH.mtx"))


def load_matrix(fn):
    mtx = None
    try:
        mtx = sp.io.mmread(fn)
    except:
        mtx = io.mmread(fn)
    mtx = mtx.toarray()    # convert from sparse matrix to ndarray to support slicing.
    return(mtx)


def load_regions(fn):
    df = pd.read_csv(fn, header = None, sep = "\t")
    df.columns = ["chrom", "start", "end", "reg_id"]
    return(df)


def load_samples(fn):
    df = pd.read_csv(fn, header = None)
    df.columns = ["cell"]
    return(df)


def save_matrix(mtx, fn):
    mtx = sparse.csr_matrix(mtx)   # convert from ndarray to sparse matrix to be fully compatible with .mtx format.
    io.mmwrite(fn, mtx)


def save_regions(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)


def save_samples(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)
