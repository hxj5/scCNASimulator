# allele.py - allele-specific information


import os
import sys
from sys import stdout, stderr


class AlleleUMI:
    def __init__(self):
        self.dat = {}

    def add(self, cell, reg_id, umi, allele):
        """Add one allele-specific UMI
        @param cell     The cell barcode [str]
        @param reg_id   The region ID [str]
        @param umi      The UMI [str]
        @param allele   The index of allele, one of {0, 1} [int]
        @return         The retcode [int]
                        0 if UMI is aligned to one allele;
                        1 if UMI is aligned to both alleles and will be discarded;
                        -1 if error;
        """
        if cell not in self.dat:
            self.dat[cell] = {}

        if umi not in self.dat[cell]:
            self.dat[cell][umi] = [0, {}]    # state, ale_dat

        ale_dat = self.dat[cell][umi][1]
        n_ale = len(ale_dat)

        if n_ale == 0:
            assert self.dat[cell][umi][0] == 0
            ale_dat[allele] = set([reg_id])
        elif n_ale == 1:
            assert self.dat[cell][umi][0] == 0
            if allele in ale_dat:
                ale_dat[allele].add(reg_id)
            else:     # UMI aligned to both alleles
                ale_dat[allele] = set([reg_id])
                self.dat[cell][umi][0] = 1
        else:
            assert self.dat[cell][umi][0] == 1
            ale_dat[allele].add(reg_id)

        return(self.dat[cell][umi][0])

    def query(self, cell, umi):
        """Query the allele index
        @param cell     The cell barcode [str]
        @param umi      The UMI [str]
        @return         A tuple of two elements if there is matched cell and UMI,
                          allele: allele index, one of {0, 1};
                          A list of region IDs.
                        Otherwise, None.
        """
        if cell not in self.dat:
            return(None)
        if umi not in self.dat[cell]:
            return(None)
        state, ale_dat = self.dat[cell][umi]
        if state == 0:
            assert len(ale_dat) == 1
            allele = list(ale_dat.keys())[0]
            return((allele, list(ale_dat[allele])))
        else:
            return(None)
        

def load_allele_umi(fn_list, verbose = False):
    func = "load_allele_umi"

    au = AlleleUMI()

    for fn in fn_list:
        if not os.path.exists(fn):
            stderr.write("[E::%s] '%s' does not exist in umi_dir.\n" % (func, fn))
            raise OSError
        fp = open(fn, "r")
        for line in fp:
            items = line.strip().split("\t")
            cell, reg_id, umi, ale_idx = items[:4]
            cell = cell.strip('"')
            reg_id = reg_id.strip('"')
            umi = umi.strip('"')
            ale_idx = int(ale_idx)
            ret = au.add(cell, reg_id, umi, ale_idx)
            if ret != 0:
                if verbose:
                    stderr.write("[W::%s] UMI '%s' is not uniquely aligned to one allele.\n" % \
                            (func, umi))
        fp.close()

    return(au)

