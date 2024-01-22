# Region routine
# Author: Xianjie Huang


from sys import stdout, stderr

from ..blib.region import Region, RegionSet
from ..blib.zfile import zopen


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


def load_snp_from_tsv(fn, verbose = False):
    """Load phased SNPs from TSV file.
    @param fn       Path to TSV file containing 6 columns without header [str]:
                    <chrom> <pos> <ref> <alt> <ref_hap> <alt_hap>
    @param verbose  If print detailed log info [bool]
    @return         A SNPSet object if success, None otherwise.
    """
    func = "load_snp_from_tsv"
    fp = zopen(fn, "rt")
    snp_set = SNPSet()
    nl = 0
    if verbose:
        stderr.write("[I::%s] start to load SNPs from tsv '%s' ...\n" % (func, fn))
    for line in fp:
        nl += 1
        #if nl == 1:
        #    continue
        parts = line.rstrip().split("\t")
        if len(parts) < 6:
            if verbose:
                stderr.write("[W::%s] too few columns of line %d.\n" % (func, nl))
            continue
        ref, alt = parts[2].upper(), parts[3].upper()
        if len(ref) != 1 or ref not in "ACGTN":
            if verbose:
                stderr.write("[W::%s] invalid REF base of line %d.\n" % (func, nl))
            continue
        if len(alt) != 1 or alt not in "ACGTN":
            if verbose:
                stderr.write("[W::%s] invalid ALT base of line %d.\n" % (func, nl))
            continue
        a1, a2 = parts[4], parts[5]
        if (a1 == "0" and a2 == "1") or (a1 == "1" and a2 == "0"):
            snp = SNP(
                chrom = parts[0], 
                pos = int(parts[1]), 
                ref = ref, 
                alt = alt, 
                ref_idx = int(a1), 
                alt_idx = int(a2)
            )
            if snp_set.add(snp) < 0:
                if verbose:
                    stderr.write("[E::%s] failed to add SNP of line %d.\n" % (func, nl))
                return None
        else:
            if verbose:
               stderr.write("[W::%s] invalid GT of line %d.\n" % (func, nl))
            continue          
    fp.close()
    return snp_set


def load_snp_from_vcf(fn, verbose = False):
    """Load phased SNPs from VCF file.
    @param fn       Path to VCF file [str]
    @param verbose  If print detailed log info [bool]
    @return         A SNPSet object if success, None otherwise.
    """
    func = "load_snp_from_vcf"
    fp = zopen(fn, "rt")
    snp_set = SNPSet()
    nl = 0
    if verbose:
        stderr.write("[I::%s] start to load SNPs from vcf '%s' ...\n" % (func, fn))
    for line in fp:
        nl += 1
        if line[0] in ("#", "\n"):
            continue
        parts = line.rstrip().split("\t")
        if len(parts) < 10:
            if verbose:
                stderr.write("[W::%s] too few columns of line %d.\n" % (func, nl))
            continue
        ref, alt = parts[3].upper(), parts[4].upper()
        if len(ref) != 1 or ref not in "ACGTN":
            if verbose:
                stderr.write("[W::%s] invalid REF base of line %d.\n" % (func, nl))
            continue
        if len(alt) != 1 or alt not in "ACGTN":
            if verbose:
                stderr.write("[W::%s] invalid ALT base of line %d.\n" % (func, nl))
            continue          
        fields = parts[8].split(":")
        if "GT" not in fields:
            if verbose:
                stderr.write("[W::%s] GT not in line %d.\n" % (func, nl))
            continue
        idx = fields.index("GT")
        values = parts[9].split(":")
        if len(values) != len(fields):
            if verbose:
               stderr.write("[W::%s] len(fields) != len(values) in line %d.\n" % (func, nl))
            continue
        gt = values[idx]
        sep = ""
        if "|" in gt: 
            sep = "|"
        elif "/" in gt: 
            sep = "/"
        else:
            if verbose:
               stderr.write("[W::%s] invalid delimiter of line %d.\n" % (func, nl))
            continue
        a1, a2 = gt.split(sep)[:2]
        if (a1 == "0" and a2 == "1") or (a1 == "1" and a2 == "0"):
            snp = SNP(
                chrom = parts[0], 
                pos = int(parts[1]), 
                ref = ref, 
                alt = alt, 
                ref_idx = int(a1), 
                alt_idx = int(a2)
            )
            if snp_set.add(snp) < 0:
                if verbose:
                    stderr.write("[E::%s] failed to add SNP of line %d.\n" % (func, nl))
                return None
        else:
            if verbose:
               stderr.write("[W::%s] invalid GT of line %d.\n" % (func, nl))
            continue          
    fp.close()
    return snp_set


class BlockRegion(Region):
    """Block Region
    @param chrom    Chromosome name [str]
    @param start    1-based start pos, inclusive [int]
    @param end      1-based end pos, exclusive [int]
    @param name     Name of the block [str]
    @param snp_list A list of SNPs located within the block [list of SNP objects]
    """
    def __init__(self, chrom, start, end, name = None, snp_list = None):
        super().__init__(chrom, start, end, name)
        self.name = name
        self.snp_list = snp_list


def load_region_from_txt(fn, sep = "\t", verbose = False):
    """Load regions from file.
    @param fn       Path to file [str]
    @param verbose  If print detailed log info [bool]
    @return         A list of BlockRegion objects if success, None otherwise.
    @note           The first 4 columns of the file should be
                    chrom, start, end (both 1-based, inclusive), name
    """
    func = "load_region_from_txt"
    fp = zopen(fn, "rt")
    reg_list = []
    rs = RegionSet(is_uniq = True)
    nl = 0
    if verbose:
        stdout.write("[I::%s] start to load regions from file '%s' ...\n" % (func, fn))
    for line in fp:
        nl += 1
        parts = line.rstrip().split(sep)
        if len(parts) < 4:
            if verbose:
                stderr.write("[E::%s] too few columns of line %d.\n" % (func, nl))
            return None           
        chrom, start, end, name = parts[:4]
        start, end = int(start), int(end)
        reg = BlockRegion(chrom, start, end + 1, name)
        ret = rs.add(reg)
        if ret != 0:
            if verbose:
                stderr.write("[W::%s] retcode %d for adding region '%s:%s-%s'.\n" % \
                    (func, ret, chrom, start, end))
    fp.close()
    reg_list = rs.get_regions()
    if verbose:
        stdout.write("[I::%s] %d regions loaded.\n" % (func, len(reg_list)))

    return reg_list


if __name__ == "__main__":
    pass

