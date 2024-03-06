# core.py - core part of simulation part


import numpy as np
from sys import stdout, stderr


class UMICount:
    """
    UMI counting in CNV regions.
    """
    def __init__(self):
        self.data = {}

        self.ALE_A0 = 0
        self.ALE_A1 = 1
        self.ALE_AMB = 2
        self.ALE_INV = -1    # invalid

        self.ST_DEL = -1
        self.ST_NEU = 0
        self.ST_AMP = 1

    def __set_allele(self, cell, umi, allele = None):
        self.__set_default(cell, umi)
        self.data[cell][umi]["allele"] = allele

    def __set_default(self, cell, umi):
        if cell not in self.data:
            self.data[cell] = {}
        if umi not in self.data[cell]:
            self.data[cell][umi] = {
                "allele": None,
                "cn": None,
                "region": None,
                "state": None
            }

    def get_allele(self, cell, umi):
        if cell not in self.data or umi not in self.data[cell]:
            return(None)
        else:
            return(self.data[cell][umi]["allele"])

    def get_copy_number(self, cell, umi):
        if cell not in self.data or umi not in self.data[cell]:
            return(None)
        else:
            return(self.data[cell][umi]["cn"])

    def get_region(self, cell, umi):
        if cell not in self.data or umi not in self.data[cell]:
            return(None)
        else:
            return(self.data[cell][umi]["region"])

    def is_allele_a0(self, allele):
        return allele == self.ALE_A0

    def is_allele_a1(self, allele):
        return allele == self.ALE_A1

    def is_allele_ambiguous(self, allele):
        return allele == self.ALE_AMB

    def is_allele_invalid(self, allele):
        return allele == self.ALE_INV

    def is_state_amplify(self, state):
        return state == self.ST_AMP

    def is_state_delete(self, state):
        return state == self.ST_DEL

    def is_state_neutral(self, state):
        return state == self.ST_NEU

    def isin(self, cell, umi):
        return(cell in self.data and umi in self.data[cell])

    def set_allele_a0(self, cell, umi):
        self.__set_allele(cell, umi, self.ALE_A0)

    def set_allele_a1(self, cell, umi):
        self.__set_allele(cell, umi, self.ALE_A1)

    def set_allele_ambiguous(self, cell, umi):
        self.__set_allele(cell, umi, self.ALE_AMB)

    def set_allele_invalid(self, cell, umi):
        self.__set_allele(cell, umi, self.ALE_INV)

    def set_copy_number(self, cell, umi, cn):
        self.__set_default(cell, umi)
        self.data[cell][umi]["cn"] = cn

    def set_region(self, cell, umi, reg_id):
        self.__set_default(cell, umi)
        self.data[cell][umi]["region"] = reg_id

    def set_state_amplify(self, cell, umi):
        self.__set_default(cell, umi)
        self.data[cell][umi]["state"] = self.ST_AMP

    def set_state_delete(self, cell, umi):
        self.__set_default(cell, umi)
        self.data[cell][umi]["state"] = self.ST_DEL

    def set_state_neutral(self, cell, umi):
        self.__set_default(cell, umi)
        self.data[cell][umi]["state"] = self.ST_NEU

    def stat(self):
        func = "UMICount::stat"
        res = {}
        for cell, c_dat in self.data.items():
            if cell not in res:
                res[cell] = {}

            for umi, u_dat in c_dat.items():
                allele, cn, reg_id, state = [u_dat[k] for k in \
                    ("allele", "cn", "region", "state")]
                if reg_id is None:
                    reg_id = "None"
                if reg_id not in res[cell]:
                    res[cell][reg_id] = {
                        "A0": 0,
                        "A0_amp": 0,
                        "A0_del": 0,
                        "A1": 0,
                        "A1_amp": 0,
                        "A1_del": 0,
                        "AMB": 0,
                        "AMB_amp": 0,
                        "AMB_del": 0,
                        "invalid": 0,
                        "unknown": 0
                    }
                if self.is_allele_a0(allele):
                    res[cell][reg_id]["A0"] += 1
                    if self.is_state_amplify(state):
                        res[cell][reg_id]["A0_amp"] += 1
                    elif self.is_state_delete(state):
                        res[cell][reg_id]["A0_del"] += 1
                elif self.is_allele_a1(allele):
                    res[cell][reg_id]["A1"] += 1
                    if self.is_state_amplify(state):
                        res[cell][reg_id]["A1_amp"] += 1
                    elif self.is_state_delete(state):
                        res[cell][reg_id]["A1_del"] += 1
                elif self.is_allele_ambiguous(allele):
                    res[cell][reg_id]["AMB"] += 1
                    if self.is_state_amplify(state):
                        res[cell][reg_id]["AMB_amp"] += 1
                    elif self.is_state_delete(state):
                        res[cell][reg_id]["AMB_del"] += 1
                elif self.is_allele_invalid(allele):
                    res[cell][reg_id]["invalid"] += 1
                else:
                    res[cell][reg_id]["unknown"] += 1
        return(res)


# Note,
# we create new UMI barcodes by simply adding distinct suffix to the original
# UMI barcode. Alternatively, we can iterate the BAM file to get all unique UMI
# barcodes first, and then create new ones. However, the latter strategy is
# inefficient and not easy to implement.

def __write_read(read, sam, umi, umi_tag, qname = None, idx = 0, umi_suffix_len = 4):
    func = "__write_read"
    if idx < 0 or idx >= 2 ^ umi_suffix_len:
        raise ValueError("[E::%s] invalid idx '%s'." % (func, idx))

    if qname is None:
        read.query_name += "_" + str(idx)
    else:
        read.query_name = qname + "_" + str(idx)

    if not umi:
        sam.write(read)
        return

    x = ["ACGT"[(idx >> (2 * i)) & 3] for i in range(umi_suffix_len)]
    x.reverse()
    umi_suffix = "".join(x)

    umi += umi_suffix
    read.set_tag(umi_tag, umi)
    sam.write(read)


def simu_cnv(
    in_sam, out_sam,
    cell_anno, cnv_profile, allele_umi,
    cell_tag, umi_tag,
    debug = 0
):
    func = "simu_cnv"

    cnv_clones = cnv_profile.get_clones()
    cnv_clones = set(cnv_clones)

    if debug:
        stderr.write("[D::%s] name of CNV clones:\n" % func)
        stderr.write(str(cnv_clones) + "\n")

    clone = None
    allele = None
    cn = None           # copy number
    reg_id_list = None

    uc = UMICount()

    if debug:
        stderr.write("[D::%s] begin to iterate reads.\n" % func)

    for read in in_sam.fetch():
        if not read.has_tag(umi_tag):
            __write_read(read, out_sam, None, umi_tag)
            continue
        umi = read.get_tag(umi_tag)

        if not read.has_tag(cell_tag):
            __write_read(read, out_sam, umi, umi_tag)
            continue
        cell = read.get_tag(cell_tag)

        if not read.query_name:
            __write_read(read, out_sam, umi, umi_tag)
            continue
        qname = read.query_name

        if cell not in cell_anno:    
            __write_read(read, out_sam, umi, umi_tag)
            continue
        clone = cell_anno[cell]
        if clone not in cnv_clones:    # cell not in CNV clone.
            __write_read(read, out_sam, umi, umi_tag)
            continue

        # Note:
        # Currently fork or deletion in read level is not perfect
        # as it should be done in UMI level.
        # For example, one UMI may overlap with the copy loss regions, 
        # and all its reads are selected to be discarded, 
        # while in current scheme, reads of the UMI that do not overlap 
        # CNV regions will not be discarded, hence this UMI is unexpectedly 
        # kept in the output BAM file.
        # This should have little influence on the CNV simulation, 
        # as the number of these reads should be very small considering 
        # the length of CNV regions and UMIs.

        # check whether the read has allele information.
        rec_allele = uc.get_allele(cell, umi)

        if debug:
            if cell == "GATGAAAAGTGTGAAT-1" and umi in ("CCAGGACTTT", "TAGTTTGATG"):
                stderr.write("[D::%s] rec_allele:'%s' ('%s'-'%s').\n" % \
                    (func, rec_allele, cell, umi))

        if rec_allele is None:
            res = allele_umi.query(cell, umi)

            if debug:
                if cell == "GATGAAAAGTGTGAAT-1" and umi in ("CCAGGACTTT", "TAGTTTGATG"):
                    stderr.write("[D::%s] res:'%s' ('%s'-'%s').\n" % \
                        (func, str(res), cell, umi))

            if res is None:     # no allele info; region info to be checked.
                chrom = read.reference_name
                if not chrom:
                    __write_read(read, out_sam, umi, umi_tag)
                    continue
                positions = read.get_reference_positions()        # 0-based
                if not positions:
                    __write_read(read, out_sam, umi, umi_tag)
                    continue
                start, end = positions[0] + 1, positions[-1] + 1  # 1-based

                n, profile = cnv_profile.fetch(chrom, start, end + 1, clone)
                if n < 0:
                    raise ValueError("[E::%s] CNVProfile fetch failed for ambiguous UMIs ('%s'-'%s')." % \
                        (func, cell, umi))
                elif n == 0:           # not in CNV regions
                    __write_read(read, out_sam, umi, umi_tag)
                    continue
                elif n == 1:
                    cn0, cn1, reg_id = profile[0][:3]
                else:         # overlap multiple CNV regions (with distinct CNV profiles).
                    uc.set_allele_invalid(cell, umi)
                    __write_read(read, out_sam, umi, umi_tag)
                    continue

                rand_f = np.random.rand()
                if rand_f < 0.5:
                    cn = cn0
                else:
                    cn = cn1

                uc.set_allele_ambiguous(cell, umi)
                uc.set_copy_number(cell, umi, cn)
                uc.set_region(cell, umi, reg_id)

            else:               # UMI has allele info and overlaps CNV regions.
                allele, reg_id_list = res[:2]

                if debug:
                    if cell == "GATGAAAAGTGTGAAT-1" and umi in ("CCAGGACTTT", "TAGTTTGATG"):
                        stderr.write("[D::%s] allele:'%s' ('%s'-'%s').\n" % \
                            (func, allele, cell, umi))

                if len(reg_id_list) == 1:
                    if allele not in (0, 1):
                        raise ValueError("[E::%s] allele '%s' should be 0 or 1 ('%s'-'%s')." % \
                            (func, allele, cell, umi))
                else:      # the UMI has allele info; but overlap with multiple regions.
                    uc.set_allele_invalid(cell, umi)
                    __write_read(read, out_sam, umi, umi_tag)
                    continue
                reg_id = reg_id_list[0]

                if debug:
                    if cell == "GATGAAAAGTGTGAAT-1" and umi in ("CCAGGACTTT", "TAGTTTGATG"):
                        stderr.write("[D::%s] region:'%s'; clone:'%s' ('%s'-'%s').\n" % \
                            (func, reg_id, clone, cell, umi))
        
                n, profile = cnv_profile.query(reg_id, clone)

                if debug:
                    if cell == "GATGAAAAGTGTGAAT-1" and umi in ("CCAGGACTTT", "TAGTTTGATG"):
                        stderr.write("[D::%s] n:'%s'; profile:'%s' ('%s'-'%s').\n" % \
                            (func, n, str(profile), cell, umi))

                if debug:
                    if cell == "GATGAAAAGTGTGAAT-1" and umi in ("CCAGGACTTT", "TAGTTTGATG"):
                        stderr.write("[D::%s] cnv_profile.clones='%s'.\n" % \
                            (func, str(cnv_profile.dat.keys())))
                        cp = cnv_profile.dat["immune_cnv_clone_1"]
                        stderr.write("[D::%s] rs.cid='%s'.\n" % \
                            (func, str(cp.rs.cid.keys())))

                if n < 0:
                    raise ValueError("[E::%s] CNVProfile query failed for ambiguous UMIs ('%s'-'%s')." % \
                        (func, cell, umi))
                elif n == 0:      # not in CNV regions
                    __write_read(read, out_sam, umi, umi_tag)
                    continue
                elif n == 1:
                    cn0, cn1, reg_id = profile[0][:3]
                else:         # overlap multiple CNV regions (with distinct CNV profiles).
                    uc.set_allele_invalid(cell, umi)
                    __write_read(read, out_sam, umi, umi_tag)
                    continue

                if debug:
                    if cell == "GATGAAAAGTGTGAAT-1" and umi in ("CCAGGACTTT", "TAGTTTGATG"):
                        stderr.write("[D::%s] allele2:'%s' ('%s'-'%s').\n" % \
                            (func, allele, cell, umi))
                
                if allele == 0:
                    uc.set_allele_a0(cell, umi)
                    uc.set_copy_number(cell, umi, cn0)
                else:      # must be 1.
                    uc.set_allele_a1(cell, umi)
                    uc.set_copy_number(cell, umi, cn1)
                uc.set_region(cell, umi, reg_id)

        allele = uc.get_allele(cell, umi)
        assert allele is not None
        if uc.is_allele_invalid(allele):
            __write_read(read, out_sam, umi, umi_tag)
            continue
        else:
            cn = uc.get_copy_number(cell, umi)
            assert cn is not None
            if cn == 0:
                uc.set_state_delete(cell, umi)
                continue
            elif cn == 1:
                uc.set_state_neutral(cell, umi)
                __write_read(read, out_sam, umi, umi_tag)
                continue
            else:
                uc.set_state_amplify(cell, umi)
                for i in range(cn):
                    __write_read(read, out_sam, umi, umi_tag, qname, i)

    if debug:
        stderr.write("[D::%s] begin to do UMI stat.\n" % func)

    uc_stat = uc.stat()
    return(uc_stat)

