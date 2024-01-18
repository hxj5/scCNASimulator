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

