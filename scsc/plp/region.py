# Region routine
# Author: Xianjie Huang


from ..lib.region import Region, RegionSet


class SNP(Region):
    """Phased SNP
    @param chrom    Chromosome name [str]
    @param pos      1-based position [int]
    @param ref      The ref base [str]
    @param alt      The alt base [str]
    @param ref_idx  The GT index for ref base, 0 or 1 [int]
    @param alt_idx  The GT index for alt base, 1 or 0 [int]   
    """
    def __init__(self, chrom, pos, ref, alt, ref_idx, alt_idx):
        super().__init__(chrom, pos, pos + 1)
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.ref_idx = ref_idx
        self.alt_idx = alt_idx
        self.gt = {ref:ref_idx, alt:alt_idx}

    def get_id(self):
        return "%s_%d" % (self.chrom, self.pos)

    def get_region_allele_index(self, base):
        return self.gt[base] if base in self.gt else -1
    

class SNPSet(RegionSet):
    """A set of phased SNPs"""
    def __init__(self, is_uniq = False):
        super().__init__(is_uniq)

    def add(self, snp):
        return super().add(snp)


class BlockRegion(Region):
    """Block Region
    @param chrom    Chromosome name [str]
    @param start    1-based start pos, inclusive [int]
    @param end      1-based end pos, exclusive [int]
    @param name     Name of the block [str]
    @param snp_list A list of SNPs located within the block [list of SNP objects]
    """
    def __init__(self, chrom, start, end, name = None, snp_list = None):
        super().__init__(chrom, start, end)
        self.name = name
        self.snp_list = snp_list


if __name__ == "__main__":
    pass

