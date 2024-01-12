
Manual
------

Input:

#. Indexed BAM file for 10x scRNA-seq data;
#. Phased SNPs;
#. A list of genomic regions;
#. Clone-specific CNV profile;
#. Cells in each clone;

Output:

#. In stage 1, the cell x allele UMI list;
#. In stage 2, new indexed BAM file;


Phased SNPs
~~~~~~~~~~~
The input phased SNPs should be stored in either a tab delimited text file
(i.e., TSV file) or a VCF file.

Input is in TSV format
++++++++++++++++++++++

If it is in a TSV file, then the first six columns of the file should be
``chrom``, ``pos``, ``ref``, ``alt``, ``ref_hap``, ``alt_hap``, where

chrom : str
    the chromosome name, e.g., chr1.

pos : int
    the genomic position of the SNP, 1-based.

ref : str
    the reference (REF) allele of the SNP, one of ``{'A', 'C', 'G', 'T'}``.

alt : str
    the alternative (ALT) allele of the SNP, one of ``{'A', 'C', 'G', 'T'}``.

ref_hap : int
    the haplotype index of the reference allele, one of ``{0, 1}``.

alt_hap : int
    the haplotype index of the alternative allele, one of ``{0, 1}``.

Note that **this input file should not contain a header line**. 
An example is as follows:

.. code-block::

  chr1     120    A     C       0       1
  chr1     260    C     T       1       0
  chr2     580    A     G       1       0

Input is in VCF format
++++++++++++++++++++++

If it is in VCF format, the file should contain the ``GT`` field in 
``FORMAT`` (i.e., the 9th column).
The corresponding phased genotype could be delimited by either ``'/'`` or
``'|'``, e.g., "0/1", or "0|1".

.. note::
   As reference phasing, e.g., with Eagle2, is not perfect, one UMI may 
   cover two SNPs with conflicting haplotype states.


A list of genomic regions
~~~~~~~~~~~~~~~~~~~~~~~~~
The input regions shold be stored in a tab delimited text file with the first
four columns being ``chrom``, ``start``, ``end``, ``name``, where

chrom : str
    the chromosome name, e.g., chr1.

start : int
    the start genomic position of the region, 1-based and inclusive.

end : int
    the end genomic position of the region, 1-based and inclusive.

name : str
    the name of the region.

Note that **this input file should not contain a header line**. 
An example is as follows:

.. code-block::

  chr1     1000    3600    gene_1
  chr1     6100    9500    gene_2
  chr2     3500    7200    gene_3


.. note::
   Currently, scsc extract the allele-specific UMIs in each genomic region
   (typically gene), and uses multi-thread to process the genes in parallel,
   for higher computational efficiency, in terms of both running speed and 
   memory usage.
   This strategy should generally work well and is expected to give highly
   concordant results to sequentially pileup reads.
   However, it may still lead to some issues, e.g., discarding UMIs due to poor
   gene annotation. 


Clone-specific CNV profile
~~~~~~~~~~~~~~~~~~~~~~~~~~



Cells in each clone
~~~~~~~~~~~~~~~~~~~


The cell x allele UMI list
~~~~~~~~~~~~~~~~~~~~~~~~~~


