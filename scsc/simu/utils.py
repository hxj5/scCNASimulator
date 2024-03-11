# utils.py - utils


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
    with open(fn, "w") as fp:
        cell_list = list(anno.keys())
        if sort_cells:
            cell_list = sorted(cell_list)
        for cell in cell_list:
            cell_type = anno[cell]
            s = "%s\t%s\n" % (cell, cell_type)
            fp.write(s)
