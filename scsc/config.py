# config.py - global configuration


from .region import format_chrom


class Config:
    def __ini__(self):
        self.sid = None

        self.bam_fn = None
        self.chrom = None
        self.cell_tag = None
        self.umi_tag = None

        self.umi_dir = None

        self.global_phase_fn = None

        self.cell_anno_fn = None
        self.target_cell_types = None

        self.clone_frac0 = None
        self.clone_prob = None
        
        self.seed = None

        self.gene_anno_fn = None

        self.out_dir = None


    def check_args(self):
        assert_n(self.sid)

        assert_e(self.bam_fn)
        assert_n(self.chrom)
        self.chrom = format_chrom(self.chrom)
        assert_n(self.cell_tag)
        assert_n(self.umi_tag)

        assert_e(self.umi_dir)

        assert_e(self.global_phase_fn)

        assert_e(self.cell_anno_fn)
        assert_n(self.target_cell_types)
        self.target_cell_types = set(self.target_cell_types.split(";"))

        assert_notnone(self.clone_frac0)
        assert_n(self.clone_prob)
        clone_prob = self.clone_prob.split(",")
        if len(clone_prob) != 2:
            raise ValueError
        self.clone_prob = [float(x) for x in clone_prob]
        
        assert_notnone(self.seed)

        assert_e(self.gene_anno_fn)

        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)

        
APP = "scsc"
VERSION = "0.0.1"
AUTHOR = "Xianjie Huang"

