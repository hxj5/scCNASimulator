
A list of genomic regions
~~~~~~~~~~~~~~~~~~~~~~~~~
The input regions shold be stored in a tab delimited text file with the first
4 columns being ``chrom``, ``start``, ``end``, ``name``, where

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

