# utils.py - utils

from ..utils.region import Region, RegionSet
from ..utils.zfile import zopen


def load_cell_anno(fn):
    cell_anno = {}
    with open(fn, "r") as fp:
        for line in fp:
            items = line.strip().split("\t")
            cell, cell_type = items[:2]
            cell, cell_type = cell.strip('"'), cell_type.strip('"')
            cell_anno[cell] = cell_type
    return(cell_anno)


def save_cell_anno(anno, fn, sort_cells = True):
    cell_list = list(anno.keys())
    if sort_cells:
        cell_list = sorted(cell_list)
    with open(fn, "w") as fp:
        for cell in cell_list:
            cell_type = anno[cell]
            s = "%s\t%s\n" % (cell, cell_type)
            fp.write(s)


def load_features(fn, sep = "\t", verbose = False):
    func = "load_features"
    fp = zopen(fn, "r")
    nl = 0
    rs = RegionSet(is_uniq = True)
    for line in fp:
        nl += 1
        items = line.strip().strip('"').split(sep)
        if len(items) < 4:
            raise IOError("[E::%s] missing fields in line %d." % (func, nl))
        chrom, start, end, reg_id = items[:4]
        start, end = int(start), int(end)
        reg_id = reg_id.strip().strip('"')
        reg = Region(chrom, start, end, reg_id)
        rs.add(reg)
    fp.close()
    return(rs)


def save_features(data, fn):
    func = "save_features"
    regions = data.get_regions(sort = True)
    with open(fn, "w") as fp:
        for reg in regions:
            s = "%s\t%s\t%s\t%s\n" % (reg.chrom, reg.start + 1, reg.end, reg.get_id())
            fp.write(s)
