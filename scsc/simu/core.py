# core.py - core part of simulation part


class UMICount:
    def __init__(self):
        self.data = {}

        self.A0 = 0
        self.A1 = 1
        self.AMB = 2

        self.DEL = -1
        self.NEU = 0
        self.AMP = 1

    def get_allele(self, cell, umi):
        if cell not in self.data or umi not in self.data[cell]:
            return(None)
        else:
            return(self.data[cell][umi]["allele"])

    def get_cn(self, cell, umi):
        if cell not in self.data or umi not in self.data[cell]:
            return(None)
        else:
            return(self.data[cell][umi]["cn"])

    def isin(self, cell, umi):
        return(cell in self.data and umi in self.data[cell])

    def set_allele(self, cell, umi, allele = None):
        if cell not in self.data:
            self.data[cell] = {}
        if umi not in self.data[cell]:

    def set_cn(self, cell, umi, cn):
        pass

    def to_amplify(self, cell, umi):
        pass

    def to_delete(self, cell, umi):
        pass


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

    clone = None
    cn = None           # copy number
    amb_umi = {}        # ambiguous UMIs that have been processed

    umi_stat = {}

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

        if cell not in cell_anno:    # cell not in CNV clone.
            __write_read(read, out_sam, umi, umi_tag)
            continue
        clone = cell_anno[cell]

        # check whether the read has allele information.

        res = allele_umi.query(cell, umi)

        if res is None:    # ambiguous UMIs

            # check whether read is in CNV region
            # for copy gain and copy loss, ambiguous UMIs may still be forked 
            # or deleted.

            # Note:
            # Currently fork or deletion in read level is not perfect
            # as it should be done in UMI level.
            # For example, one UMI may overlap with the copy loss regions, 
            # and all its reads are selected to be discarded, 
            # while in current scheme, reads of the UMI that do not overlap 
            # CNV regions will not be discarded, 
            # hence this UMI is unexpectedly kept in the output BAM file.
            # This should have little influence on the CNV simulation, 
            # as the number of these reads should be very small considering 
            # the length of CNV regions and UMIs.

            if cell in amb_umi and umi in amb_umi[cell]:    # already processed the UMI
                cn = amb_umi[cell][umi]
            else:
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
                    raise ValueError
                if ret != 0:           # not in CNV region or overlap multiple CNV regions with distinct CNV profiles.
                    __write_read(read, out_sam, umi, umi_tag)
                    continue
                cn0, cn1 = cn_ale[:2]

                rand_f = np.random.rand()
                if rand_f < 0.5:
                    cn = cn0
                else:
                    cn = cn1

                if cell not in amb_umi:
                    amb_umi[cell] = {}
                amb_umi[cell][umi] = cn     # as umi not in amb_umi[cell]

            if cn <= 0:
                continue
            elif cn == 1:
                __write_read(read, out_sam, umi, umi_tag)
                continue
            else:
                # update the QNAME and UMI of the forked reads.
                for i in range(cn):
                    __write_read(read, out_sam, umi, umi_tag, qname, i)

        else:
            allele, reg_id_list = res[:2]
            if len(reg_id_list) != 1:
                __write_read(read, out_sam, umi, umi_tag)
                continue
            reg_id = reg_id_list[0]

            ret, cn_ale = cnv_profile.query(reg_id, clone)
            if ret < 0:
                raise ValueError
            if ret != 0:
                __write_read(read, out_sam, umi, umi_tag)
                continue
            cn = cn_ale[allele]     # copy number

            if cn == 0:
                continue
            elif cn == 1:
                __write_read(read, out_sam, umi, umi_tag)
                continue
            else:
                for i in range(cn):
                    __write_read(read, out_sam, umi, umi_tag, qname, i)
