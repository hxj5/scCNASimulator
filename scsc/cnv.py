# cnv.py - CNV routines


from .blib.region import Region, RegionSet


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

