# Utils

import os
from sys import stdout, stderr

from ..blib.region import Region, RegionSet
from ..blib.zfile import zopen
from .region import load_snp_from_tsv, load_snp_from_vcf


def _fmt_line(ln, k):
        items = ln.split("\t")
        items[0] = str(int(items[0]) + k)
        return("\t".join(items))


# internal use only!
def merge_mtx(in_fn_list, in_format,
              out_fn, out_fmode, out_format,
              nrow_list, ncol, nrecord, remove = False):
    if len(in_fn_list) != len(nrow_list):
        return(-1)
    bufsize = 1048576    # 1M

    is_bytes = "b" in out_fmode
    out_fp = zopen(out_fn, out_fmode, out_format, is_bytes = is_bytes)

    nrow_total = sum(nrow_list)
    header  = "%%MatrixMarket matrix coordinate integer general\n"
    header += "%%\n"
    header += "%d\t%d\t%d\n" % (nrow_total, ncol, nrecord)
    if is_bytes:
        header = bytes(header, "utf8")
    out_fp.write(header)

    nline = 0
    k = 0
    for in_fn, nrow in zip(in_fn_list, nrow_list):
        with zopen(in_fn, "rt", in_format) as in_fp:
            while True:
                lines = in_fp.readlines(bufsize)
                if not lines:
                    break
                nline += len(lines)
                lines = [_fmt_line(ln, k) for ln in lines]
                s = "".join(lines)
                if is_bytes:
                    s = bytes(s, "utf8")
                out_fp.write(s)
        k += nrow
    out_fp.close()
    if nline != nrecord:
        return(-1)
    if remove:
        for in_fn in in_fn_list:
            os.remove(in_fn)
    return(0) 


# internal use only!
def merge_tsv(in_fn_list, in_format, 
              out_fn, out_fmode, out_format, 
              remove = False):
    bufsize = 1048576   # 1M
    is_bytes = "b" in out_fmode
    out_fp = zopen(out_fn, out_fmode, out_format, is_bytes)
    in_fmode = "rb" if is_bytes else "rt"
    for in_fn in in_fn_list:
        with zopen(in_fn, in_fmode, in_format) as in_fp:
            while True:
                dat = in_fp.read(bufsize)
                if not dat:
                    break
                out_fp.write(dat)
    out_fp.close()
    if remove:
        for in_fn in in_fn_list:
            os.remove(in_fn)
    return(0)


# internal use only!
def rewrite_mtx(in_fn, in_format, 
                out_fn, out_fmode, out_format, 
                nrow, ncol, nrecord, remove = False):
    if in_fn == out_fn:
        return(-1)
    bufsize = 1048576   # 1M
  
    is_bytes = "b" in out_fmode
    out_fp = zopen(out_fn, out_fmode, out_format, is_bytes = is_bytes)
    header  = "%%MatrixMarket matrix coordinate integer general\n"
    header += "%%\n"
    header += "%d\t%d\t%d\n" % (nrow, ncol, nrecord)
    if is_bytes:
        header = bytes(header, "utf8")
    out_fp.write(header)

    nline = 0
    in_fmode = "rb" if is_bytes else "rt"
    sep = b"" if is_bytes else ""
    with zopen(in_fn, in_fmode, in_format) as in_fp:
        while True:
            lines = in_fp.readlines(bufsize)
            if not lines:
                break
            nline += len(lines)
            s = sep.join(lines)
            out_fp.write(s)
    out_fp.close()
    if nline != nrecord:
        return(-1)
    if remove:
        os.remove(in_fn)
    return(0)


class XClonePhasedGene(Region):
    def __init__(self, chrom, start, end, name, allele_flip):
        super().__init__(chrom, start, end, name)
        self.allele_flip = allele_flip


def load_xclone_phased_gene(fn, sep = "\t", verbose = False):
    func = "load_xclone_phased_gene"
    fp = zopen(fn, "r")
    nl = 0
    rs = RegionSet(is_uniq = True)
    for line in fp:
        nl += 1
        if nl == 1:
            continue
        items = line.strip().split(sep)
        if len(items) < 5:
            raise IOError
        gene_name, chrom, start, end, allele_flip = items[0], items[2], items[3], items[4], items[-1]
        start, end = int(start), int(end)
        if allele_flip == "True":
            allele_flip = True
        elif allele_flip == "False":
            allele_flip = False
        else:
            raise ValueError
        reg = XClonePhasedGene(chrom, start, end + 1, gene_name, allele_flip)
        rs.add(reg)
    fp.close()
    return(rs)


def __format_xclone_snp(snp):
    s = "\t".join([snp.chrom, str(snp.pos), snp.ref, snp.alt, str(snp.ref_idx), str(snp.alt_idx)])
    return(s)


def get_fine_xclone_phased_snp(snp_fn, gene_fn, out_fn, verbose = False):
    func = "get_fine_xclone_phased_snp"
    snp_rs = load_snp_from_vcf(snp_fn, verbose = verbose)
    gene_rs = load_xclone_phased_gene(gene_fn, sep = "\t", verbose = verbose)
    out_fp = open(out_fn, "w")

    snp_list = snp_rs.get_regions(sort = True)
    for snp in snp_list:
        hits = gene_rs.fetch(snp.chrom, snp.start, snp.end)
        if not hits:
            s = __format_xclone_snp(snp)
            out_fp.write(s + "\n")
        else:     # CHECK ME! what if more than 1 hits?
            reg = hits[0]
            if reg.allele_flip:
                snp.ref_idx = 1 - snp.ref_idx
                snp.alt_idx = 1 - snp.alt_idx
            s = __format_xclone_snp(snp)
            out_fp.write(s + "\n")
    out_fp.close()

