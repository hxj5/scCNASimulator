# pileup.py - pileup reads or UMIs from scRNA-seq BAM file.
#
# This file and files in the folder "plp" were initially designed for the 
# "xcltk pileup" module in python package "xcltk" for XClone preprocessing.
#
# Note: 
# - now it only supports 10x data.


import getopt
import multiprocessing
import os
import pickle
import sys
import time

from sys import stdout, stderr

from .app import APP, VERSION
from .plp.config import Config
from .plp.core import sp_count
from .plp.region import load_region_from_txt, load_snp_from_vcf, \
                        load_snp_from_tsv
from .plp.thread import ThreadData
from .plp.utils import merge_mtx, merge_tsv, rewrite_mtx
from .utils.xlog import log_debug, log_err, log_info
from .utils.zfile import zopen, ZF_F_GZIP, ZF_F_PLAIN


def prepare_config(conf):
    """Prepare configures for downstream analysis
    @param conf  A Config object.
    @return      0 if success, -1 otherwise.
    @note        This function should be called after cmdline is parsed.
    """
    if conf.sam_fn:
        if not os.path.isfile(conf.sam_fn):
            log_err("sam file '%s' does not exist." % conf.sam_fn)
            return(-1)
    else:
        log_err("sam file(s) needed!")
        return(-1)

    if conf.barcode_fn:
        if os.path.isfile(conf.barcode_fn):
            with zopen(conf.barcode_fn, "rt") as fp:
                conf.barcodes = sorted([x.strip() for x in fp])   # UPDATE!! use numpy or pandas to load
            if len(set(conf.barcodes)) != len(conf.barcodes):
                log_err("duplicate barcodes!")
                return(-1)
        else:
            log_err("barcode file '%s' does not exist." % conf.barcode_fn)
            return(-1)
    else:       
        conf.barcodes = None
        log_err("barcode file needed!")
        return(-1)

    if not conf.out_dir:
        log_err("out dir needed!")
        return(-1)
    if not os.path.isdir(conf.out_dir):
        os.mkdir(conf.out_dir)
    conf.out_region_fn = os.path.join(conf.out_dir, conf.out_prefix + "region.tsv")
    conf.out_sample_fn = os.path.join(conf.out_dir, conf.out_prefix + "samples.tsv")
    conf.out_ad_fn = os.path.join(conf.out_dir, conf.out_prefix + "AD.mtx")
    conf.out_dp_fn = os.path.join(conf.out_dir, conf.out_prefix + "DP.mtx")
    conf.out_oth_fn = os.path.join(conf.out_dir, conf.out_prefix + "OTH.mtx")

    if conf.region_fn:
        if os.path.isfile(conf.region_fn): 
            conf.reg_list = load_region_from_txt(conf.region_fn, verbose = True)
            if not conf.reg_list:
                log_err("failed to load region file.")
                return(-1)
            log_info("count %d regions in %d single cells." % (
                len(conf.reg_list), len(conf.barcodes)))
        else:
            log_err("region file '%s' does not exist." % conf.region_fn)
            return(-1)
    else:
        log_err("region file needed!")
        return(-1)

    if conf.snp_fn:
        if os.path.isfile(conf.snp_fn):
            if conf.snp_fn.endswith(".vcf") or conf.snp_fn.endswith(".vcf.gz") \
                    or conf.snp_fn.endswith(".vcf.bgz"):
                conf.snp_set = load_snp_from_vcf(conf.snp_fn, verbose = True)
            else:
                conf.snp_set = load_snp_from_tsv(conf.snp_fn, verbose = True)
            if not conf.snp_set or conf.snp_set.get_n() <= 0:
                log_err("failed to load snp file.")
                return(-1)
            else:
                log_info("%d SNPs loaded." % conf.snp_set.get_n())       
        else:
            log_err("snp file '%s' does not exist." % conf.snp_fn)
            return(-1)      
    else:
        log_err("SNP file needed!")
        return(-1)

    if conf.cell_tag and conf.cell_tag.upper() == "NONE":
        conf.cell_tag = None
    if conf.cell_tag and conf.barcodes:
        pass       
    elif (not conf.cell_tag) ^ (not conf.barcodes):
        log_err("should not specify cell_tag or barcodes alone.")
        return(-1)
    else:
        log_err("should specify cell_tag and barcodes.")
        return(-1)        

    if conf.umi_tag:
        if conf.umi_tag.upper() == "AUTO":
            if conf.barcodes is None:
                conf.umi_tag = None
            else:
                conf.umi_tag = conf.defaults.UMI_TAG_BC
        elif conf.umi_tag.upper() == "NONE":
            conf.umi_tag = None
    else:
        log_err("umi tag needed!")
        return(-1)

    if conf.barcodes:
        with open(conf.out_sample_fn, "w") as fp:
            fp.write("".join([b + "\n" for b in conf.barcodes]))

    if conf.excl_flag < 0:
        if conf.use_umi():
            conf.excl_flag = conf.defaults.EXCL_FLAG_UMI
        else:
            conf.excl_flag = conf.defaults.EXCL_FLAG_XUMI

    return(0)


def show_progress(rv = None):
    return(rv)


def pileup_core(argv, conf):
    ret = -1
    cmdline = "unexpected command line"

    start_time = time.time()

    try:
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time))
        log_info("start time: %s." % time_str)

        cmdline = " ".join(argv)
        log_info("CMD: %s" % cmdline)

        if prepare_config(conf) < 0:
            raise ValueError("errcode -2")
        log_info("program configuration:")
        conf.show(fp = stderr, prefix = "\t")

        conf.umi_dir = os.path.join(conf.out_dir, "umi")
        if not os.path.exists(conf.umi_dir):
            os.mkdir(conf.umi_dir)

        # extract SNPs for each region
        if conf.debug > 0:
            log_debug("extract SNPs for each region.")
        reg_list = []
        for reg in conf.reg_list:
            snp_list = conf.snp_set.fetch(reg.chrom, reg.start, reg.end)
            if snp_list and len(snp_list) > 0:
                reg.snp_list = snp_list
                reg_list.append(reg)
            else:
                if conf.debug > 0:
                    log_debug("no SNP fetched for region '%s'." % reg.name)
        log_info("%d regions extracted with SNPs." % len(reg_list))

        if not conf.output_all_reg:
            conf.reg_list = reg_list

        # split region list and save to file       
        m_reg = len(conf.reg_list)
        m_thread = conf.nproc if m_reg >= conf.nproc else m_reg

        reg_fn_list = []
        n_reg = m_reg // m_thread
        r_reg = m_reg - n_reg * m_thread
        k_reg = 0
        i_thread = 0
        while k_reg <= m_reg - 1:
            t_reg = n_reg + 1 if i_thread < r_reg else n_reg
            reg_fn = conf.out_prefix + "region.pickle." + str(i_thread)
            reg_fn = os.path.join(conf.out_dir, reg_fn)
            reg_fn_list.append(reg_fn)
            with open(reg_fn, "wb") as fp:
                pickle.dump(conf.reg_list[k_reg:(k_reg + t_reg)], fp)
            k_reg += t_reg
            i_thread += 1
        for reg in conf.reg_list:  # save memory
            del reg
        conf.reg_list.clear()
        conf.reg_list = None
        conf.snp_set.destroy()
        conf.snp_set = None

        thdata_list = []
        ret_sp = -1
        pool = multiprocessing.Pool(processes = m_thread)
        mp_result = []
        for i in range(m_thread):
            thdata = ThreadData(
                idx = i, conf = conf,
                reg_obj = reg_fn_list[i], is_reg_pickle = True,
                out_region_fn = conf.out_region_fn + "." + str(i),
                out_ad_fn = conf.out_ad_fn + "." + str(i),
                out_dp_fn = conf.out_dp_fn + "." + str(i),
                out_oth_fn = conf.out_oth_fn + "." + str(i),
                out_fn = None
            )
            thdata_list.append(thdata)
            if conf.debug > 0:
                log_debug("data of thread-%d before sp_count:" % i)
                thdata.show(fp = stderr, prefix = "\t")
            mp_result.append(pool.apply_async(
                    func = sp_count, 
                    args = (thdata, ), 
                    callback = show_progress))   # TODO: error_callback?
        pool.close()
        pool.join()
        mp_result = [res.get() for res in mp_result]
        retcode_list = [item[0] for item in mp_result]
        thdata_list = [item[1] for item in mp_result]
        if conf.debug > 0:
            log_debug("returned values of multi-processing:")
            log_debug("\t%s" % str(retcode_list))

        # check running status of each sub-process
        for thdata in thdata_list:         
            if conf.debug > 0:
                log_debug("data of thread-%d after sp_count:" %  thdata.idx)
                thdata.show(fp = stderr, prefix = "\t")
            if thdata.ret < 0:
                raise ValueError("errcode -3")

        # merge results
        if merge_tsv(
            [td.out_region_fn for td in thdata_list], ZF_F_GZIP, 
            conf.out_region_fn, "wb", ZF_F_PLAIN, 
            remove = True
        ) < 0:
            raise ValueError("errcode -15")

        nr_reg_list = [td.nr_reg for td in thdata_list]

        if merge_mtx(
            [td.out_ad_fn for td in thdata_list], ZF_F_GZIP, 
            conf.out_ad_fn, "w", ZF_F_PLAIN,
            nr_reg_list, len(conf.barcodes), sum([td.nr_ad for td in thdata_list]),
            remove = True
        ) < 0:
            raise ValueError("errcode -17")

        if merge_mtx(
            [td.out_dp_fn for td in thdata_list], ZF_F_GZIP, 
            conf.out_dp_fn, "w", ZF_F_PLAIN,
            nr_reg_list, len(conf.barcodes), sum([td.nr_dp for td in thdata_list]),
            remove = True
        ) < 0:
            raise ValueError("errcode -19")

        if merge_mtx(
            [td.out_oth_fn for td in thdata_list], ZF_F_GZIP, 
            conf.out_oth_fn, "w", ZF_F_PLAIN,
            nr_reg_list, len(conf.barcodes), sum([td.nr_oth for td in thdata_list]),
            remove = True
        ) < 0:
            raise ValueError("errcode -21")

    except ValueError as e:
        log_err(str(e))
        log_err("Running program failed.")
        log_err("Quiting ...")
        ret = -1

    else:
        log_info("All Done!")
        ret = 0

    finally:
        log_info("CMD: %s" % cmdline)

        end_time = time.time()
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
        log_info("end time: %s" % time_str)
        log_info("time spent: %.2fs" % end_time - start_time)

    return(ret)


def usage(fp = stderr, conf = None):
    s =  "\n"
    s += "Version: %s\n" % VERSION
    s += "Usage:   %s %s <options>\n" % (APP, COMMAND)
    s += "\n" 
    s += "Options:\n"
    s += "  -s, --sam FILE         Indexed sam/bam/cram file.\n"
    s += "  -O, --outdir DIR       Output directory for sparse matrices.\n"
    s += "  -b, --barcode FILE     A plain file listing all effective cell barcode.\n"
    s += "  -R, --region FILE      A TSV file listing target regions. The first 4 columns shoud be:\n"
    s += "                         chrom, start, end (both 1-based and inclusive), name.\n"
    s += "  -P, --phasedSNP FILE   A TSV or VCF file listing phased SNPs (i.e., containing phased GT).\n"
    s += "  -h, --help             Print this message and exit.\n"
    s += "\n"
    s += "Optional arguments:\n"
    s += "  -p, --nproc INT        Number of processes [%d]\n" % conf.NPROC
    s += "      --cellTAG STR      Tag for cell barcodes [%s]\n" % conf.CELL_TAG
    s += "      --UMItag STR       Tag for UMI, set to None when reads only [%s]\n" % conf.UMI_TAG
    s += "      --minCOUNT INT     Mininum aggragated count for SNP [%d]\n" % conf.MIN_COUNT
    s += "      --minMAF FLOAT     Mininum minor allele fraction for SNP [%f]\n" % conf.MIN_MAF
    s += "      --outputAllReg     If set, output all inputted regions.\n"
    s += "      --countDupHap      If set, UMIs aligned to both haplotypes will be counted.\n"
    s += "  -D, --debug INT        Used by developer for debugging [%d]\n" % conf.DEBUG
    s += "\n"
    s += "Read filtering:\n"
    s += "  --inclFLAG INT    Required flags: skip reads with all mask bits unset [%d]\n" % conf.INCL_FLAG
    s += "  --exclFLAG INT    Filter flags: skip reads with any mask bits set [%d\n" % conf.EXCL_FLAG_UMI
    s += "                    (when use UMI) or %d (otherwise)]\n" % conf.EXCL_FLAG_XUMI
    s += "  --minLEN INT      Minimum mapped length for read filtering [%d]\n" % conf.MIN_LEN
    s += "  --minMAPQ INT     Minimum MAPQ for read filtering [%d]\n" % conf.MIN_MAPQ
    s += "  --countORPHAN     If use, do not skip anomalous read pairs.\n"
    s += "\n"

    fp.write(s)


def pileup_main(argv, conf = None):
    """Command-Line interface.
    @param argv   A list of cmdline parameters [list]
    @param conf   The plp.Config object.
    @return       0 if success, -1 otherwise [int]
    """
    if conf is None:
        conf = Config()

    if len(argv) <= 2:
        usage(stderr, conf.defaults)
        sys.exit(1)

    opts, args = getopt.getopt(
        args = argv[2:], 
        shortopts = "-s:-O:-b:-R:-P:-h-p:-D:", 
        longopts = [
            "sam=", "outdir=", "barcode=",
            "region=", "phasedSNP=", 
            "help",

            "nproc=", 
            "cellTAG=", "UMItag=", 
            "minCOUNT=", "minMAF=", "outputAllReg", "countDupHap",
            "debug=",

            "inclFLAG=", "exclFLAG=", "minLEN=", "minMAPQ=", "countORPHAN"
        ])

    for op, val in opts:
        if len(op) > 2:
            op = op.lower()
        if op in   ("-s", "--sam"): conf.sam_fn = val
        elif op in ("-O", "--outdir"): conf.out_dir = val
        elif op in ("-b", "--barcode"): conf.barcode_fn = val
        elif op in ("-R", "--region"): conf.region_fn = val
        elif op in ("-P", "--phasedsnp"): conf.snp_fn = val
        elif op in ("-h", "--help"): usage(stderr, conf.defaults); sys.exit(1)

        elif op in ("-p", "--nproc"): conf.nproc = int(val)
        elif op in (      "--celltag"): conf.cell_tag = val
        elif op in (      "--umitag"): conf.umi_tag = val
        elif op in (      "--mincount"): conf.min_count = int(val)
        elif op in (      "--minmaf"): conf.min_maf = float(val)
        elif op in (      "--outputallreg"): conf.output_all_reg = True
        elif op in (      "--countduphap"): conf.no_dup_hap = False
        elif op in ("-D", "--debug"): conf.debug = int(val)

        elif op in ("--inclflag"): conf.incl_flag = int(val)
        elif op in ("--exclflag"): conf.excl_flag = int(val)
        elif op in ("--minlen"): conf.min_len = int(val)
        elif op in ("--minmapq"): conf.min_mapq = float(val)
        elif op in ("--countorphan"): conf.no_orphan = False

        else:
            log_err("invalid option: '%s'." % op)
            return(-1)

    ret = pileup_core(argv, conf)
    return(ret)


COMMAND = "pileup"


if __name__ == "__main__":
    pileup_main(sys.argv)

