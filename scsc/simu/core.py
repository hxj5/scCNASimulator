# core.py - core part of simulation part


import numpy as np

class UMICount:
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
        res = {}
        for cell, c_dat in self.data.items():
            if cell not in res:
                res[cell] = {
                    "A0": 0,
                    "A0_amp": 0,
                    "A0_del": 0,
                    "A1": 0,
                    "A1_amp": 0,
                    "A1_del": 0,
                    "AMB": 0,
                    "AMB_amp": 0,
                    "AMB_del": 0,
                    "invalid": 0
                }
            for umi, u_dat in c_dat.items():
                ale, st = u_dat["allele"], u_dat["state"]
                if ale == self.ALE_A0:
                    res[cell]["A0"] += 1
                    if st == self.ST_AMP:
                        res[cell]["A0_amp"] += 1
                    elif st == self.ST_DEL:
                        res[cell]["A0_del"] += 1
                elif ale == self.ALE_A1:
                    res[cell]["A1"] += 1
                    if st == self.ST_AMP:
                        res[cell]["A1_amp"] += 1
                    elif st == self.ST_DEL:
                        res[cell]["A1_del"] += 1
                elif ale == self.ALE_AMB:
                    res[cell]["AMB"] += 1
                    if st == self.ST_AMP:
                        res[cell]["AMB_amp"] += 1
                    elif st == self.ST_DEL:
                        res[cell]["AMB_del"] += 1
                else:
                    res[cell]["invalid"] += 1
        return(res)


# Note,
# we create new UMI barcodes by simply adding distinct suffix to the original
# UMI barcode. Alternatively, we can iterate the BAM file to get all unique UMI
# barcodes first, and then create new ones. However, the latter strategy is
# inefficient and not easy to implement.

def __write_read(read, sam, umi, umi_tag, qname = None, idx = 0, umi_suffix_len = 4):
    if idx < 0 or idx >= 2 ^ umi_suffix_len:
        raise ValueError

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
    cell_tag, umi_tag
):
    func = "simu_cnv"

    cnv_clones = cnv_profile.get_clones()
    cnv_clones = set(cnv_clones)

    clone = None
    cn = None           # copy number
    uc = UMICount()

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
        reg_id_list = None
        allele = uc.get_allele(cell, umi)
        if allele is None:
            res = allele_umi.query(cell, umi)
            if res is None:    # ambiguous UMIs
                uc.set_allele_ambiguous(cell, umi)
            else:
                allele, reg_id_list = res[:2]
                if len(reg_id_list) == 1:
                    if allele == 0:
                        uc.set_allele_a0(cell, umi)
                    else:
                        uc.set_allele_a1(cell, umi)
                else:
                    uc.set_allele_invalid(cell, umi)
                    __write_read(read, out_sam, umi, umi_tag)
                    continue
        elif allele < 0:         # invalid allele
            __write_read(read, out_sam, umi, umi_tag)
            continue

        if allele == 0 or allele == 1:
            cn = uc.get_copy_number(cell, umi)
            if cn is None:
                if reg_id_list is None:
                    res = allele_umi.query(cell, umi)
                    ale, reg_id_list = res[:2]
                assert len(reg_id_list) == 1
                reg_id = reg_id_list[0]
    
                ret, cn_ale = cnv_profile.query(reg_id, clone)
                if ret < 0:
                    uc.set_copy_number(cell, umi, -1)
                    raise ValueError
                if ret != 0:
                    uc.set_copy_number(cell, umi, -1)
                    __write_read(read, out_sam, umi, umi_tag)
                    continue
                cn = cn_ale[allele]     # copy number
                uc.set_copy_number(cell, umi, cn)
    
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


        else:           # ambiguous UMIs
            cn = uc.get_copy_number(cell, umi)
            if cn is None:
                chrom = read.reference_name
                if not chrom:
                    __write_read(read, out_sam, umi, umi_tag)
                    continue
                positions = read.get_reference_positions()        # 0-based
                if not positions:
                    __write_read(read, out_sam, umi, umi_tag)
                    continue
                start, end = positions[0] + 1, positions[-1] + 1  # 1-based

                ret, cn_ale = cnv_profile.fetch(chrom, start, end + 1, clone)
                if ret < 0:
                    uc.set_copy_number(cell, umi, -1)
                    raise ValueError
                if ret != 0:           # not in CNV region or overlap multiple CNV regions with distinct CNV profiles.
                    uc.set_copy_number(cell, umi, -1)
                    __write_read(read, out_sam, umi, umi_tag)
                    continue
                cn0, cn1 = cn_ale[:2]

                rand_f = np.random.rand()
                if rand_f < 0.5:
                    cn = cn0
                else:
                    cn = cn1

                uc.set_copy_number(cell, umi, cn)

            if cn <= 0:
                uc.set_state_delete(cell, umi)
                continue
            elif cn == 1:
                uc.set_state_neutral(cell, umi)
                __write_read(read, out_sam, umi, umi_tag)
                continue
            else:
                # update the QNAME and UMI of the forked reads.
                uc.set_state_amplify(cell, umi)
                for i in range(cn):
                    __write_read(read, out_sam, umi, umi_tag, qname, i)

    uc_stat = uc.stat()
    return(uc_stat)

