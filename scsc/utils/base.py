# base.py


import os


def assert_e(path):
    if path is None or not os.path.exists(path):
        raise OSError


def assert_n(var):
    if var is None or (isinstance(var, str) and len(var) == 0):
        raise ValueError


def assert_notnone(var):
    if var is None:
        raise ValueError


