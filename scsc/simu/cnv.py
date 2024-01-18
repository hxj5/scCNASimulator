# cnv.py - CNV routines


import functools
import sys
from sys import stdout, stderr

from ..blib.region import format_chrom, format_start, format_end  \
                    Region, RegionSet 
from ..blib.zfile import zopen


class CNVRegCN(Region):
    """Allele-specific copy numbers of CNV region
    @param chrom    Chromosome name [str]
    @param start    1-based start pos, inclusive [int]
    @param end      1-based end pos, exclusive [int]
    @param cn_ale0  Copy Number of the first allele [int]
    @param cn_ale1  Copy Number of the second allele [int]
    """
    def __init__(self, chrom, start, end, cn_ale0, cn_ale1):
        super().__init__(self, chrom, start, end)
        self.cn_ale0 = cn_ale0
        self.cn_ale1 = cn_ale1


class CNVProfile:
    def __init__(self):
        self.rs = RegionSet()

    def add_cnv(self, chrom, start, end, cn_ale0, cn_ale1):
        """Add a new CNV profile.
        @param chrom    Chromosome name [str]
        @param start    1-based start pos, inclusive [int]
        @param end      1-based end pos, exclusive [int]
        @param cn_ale0  Copy Number of the first allele [int]
        @param cn_ale1  Copy Number of the second allele [int]
        @return         0 success, 1 discarded as duplicate, -1 error [int]
        """
        reg = CNVRegCN(chrom, start, end, cn_ale0, cn_ale1)
        ret = self.rs.add(reg)
        return(ret)

    def fetch(self, chrom, start, end):
        """Get the CNV profile for the query region.
        @param chrom    Chromosome name [str]
        @param start    1-based start pos, inclusive [int]
        @param end      1-based end pos, exclusive [int]
        @return         A tuple of two elements:
                        ret: return code [int]
                          -1, error;
                          0, overlap with only 1 region;
                          1, overlap with zero or more than 1 regions.
                        profile: a tuple of copy numbers for the first and 
                          second alleles [tuple<int, int>]
                          None if zero or more than 1 hits.
        """
        hits = self.rs.fetch(chrom, start, end)
        if not hits:
            return((1, None))
        elif len(hits) == 1:
            reg = hits[0]
            return((0, (reg.cn_ale0, reg.cn_ale1)))
        else:
            return((1, None))

    def get_all(self):
        reg_list = self.rs.get_regions(sort = True)
        dat = {
                "chrom":[],
                "start":[],
                "end":[],
                "cn_ale0":[],
                "cn_ale1":[]
            }
        for reg in reg_list:
            dat["chrom"].append(reg.chrom)
            dat["start"].append(reg.start)
            dat["end"].append(reg.end - 1)
            dat["cn_ale0"].append(reg.cn_ale0)
            dat["cn_ale1"].append(reg.cn_ale1)
        return(dat)     
        

class CloneCNVProfile:
    def __init__(self):
        self.dat = {}

    def add_cnv(self, chrom, start, end, cn_ale0, cn_ale1, clone_id):
        if clone_id not in self.dat:
            self.dat[clone_id] = CNVProfile()
        cp = self.dat[clone_id]
        ret = cp.add_cnv(chrom, start, end, cn_ale0, cn_ale1)
        return(ret)

    def fetch(self, chrom, start, end, clone_id):
        """Get the CNV profile for the query region and cell.
        @param chrom    Chromosome name [str]
        @param start    1-based start pos, inclusive [int]
        @param end      1-based end pos, exclusive [int]
        @return         see @return of @func CNVProfile::fetch()
                        ret: 11 if clone_id is invalid;
        """
        if clone_id in self.dat:
            cp = self.dat[clone_id]    # cnv profile
            ret, hits = cp.fetch(chrom, start, end)
            return((ret, hits))
        else:
            return((11, None))

    def get_all(self):
        dat_list = {
                "clone_id":[],
                "chrom":[],
                "start":[],
                "end":[],
                "cn_ale0":[],
                "cn_ale1":[]
            }
        for clone_id in sorted(self.dat.keys()):
            cp = self.dat[clone_id]
            cp_dat = cp.get_all()
            n = len(cp_dat["chrom"])
            dat_list["clone_id"].extend([clone_id] * n)
            dat_list["chrom"].extend(cp_dat["chrom"])
            dat_list["start"].extend(cp_dat["start"])
            dat_list["end"].extend(cp_dat["end"])
            dat_list["cn_ale0"].extend(cp_dat["cn_ale0"])
            dat_list["cn_ale1"].extend(cp_dat["cn_ale1"])
        return(dat_list)


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
        stderr.write("[I::%s] start to load CNV profile from file '%s' ...\n" % (func, fn))
    for line in fp:
        nl += 1
        items = line.strip().split(sep)
        if len(items) < 6:
            if verbose:
                stderr.write("[E::%s] too few columns of line %d.\n" % (func, nl))
            return(None)
        clone_id, chrom, start, end, cn_ale0, cn_ale1 = items[:6]
        start, end = format_start(start), format_end(end)
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
    """Merge adjacent regions with the same CNV profiles

    Merge adjacent regions with the same allele-specific copy number
    profile in each CNV clone.

    @param in_fn    Path to input file [str]
    @param out_fn   Path to output file [str]
    @param max_gap  The maximum gap length that is allowed between two
                    adjacent regions. `1` for strict adjacence.
    @return  Void
    """
    func = "merge_cnv_profile"

    # load data
    fp = zopen(in_fn, "rt")
    dat = {}
    nl = 0
    if verbose:
        stderr.write("[I::%s] start to merge CNV profile from file '%s' ...\n" % (func, in_fn))
    for line in fp:
        nl += 1
        items = line.strip().split(sep)
        if len(items) < 6:
            stderr.write("[E::%s] too few columns of line %d.\n" % (func, nl))
            return(None)
        clone_id, chrom, start, end, cn_ale0, cn_ale1 = items[:6]
        start, end = format_start(start), format_end(end)
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
                    stderr.write("[E::%s] distinct CNV profiles '%s', (%d, %d) and (%d, %d).\n" % 
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
                        

def save_cnv_profile(dat, fn, verbose = False):
    """Save CNV profile to file
    @param dat      A CloneCNVProfile object.
    @param fn       The output file [str]
    @param verbose  Whether to print detailed logging information [bool]
    @return Void.
    """
    func = "save_cnv_profile"
    fp = zopen(fn, "wt")
    if verbose:
        stderr.write("[I::%s] start to save CNV profile to file '%s' ...\n" % (func, fn))
    cp = dat.get_all()
    for i in range(len(cp["chrom"])):
        s = "\t".join([
                cp["clone_id"][i],
                cp["chrom"][i],
                str(cp["start"][i]),
                str(cp["end"][i]),
                str(cp["cn_ale0"][i]),
                str(cp["cn_ale1"][i])
            ]) + "\n"
        fp.write(s)
    fp.close()

