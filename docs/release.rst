
..
   History
   =======


Release v0.0.1 (27/01/2024)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Implement the three modules ``pileup``, ``simu``, and ``pipeline`` for 
CNV simulation in 10x scRNA-seq data.
* the ``pileup`` module pileups allele-specific UMIs in each single cell.
* the ``simu`` module simulates CNVs by forking or discarding UMIs (both the
  allele-specific and ambiguous UMIs) based on the given a clone-specific 
  CNV profile, and output a new indexed BAM file.
* the ``pipeline`` module is a wrapper that runs sequentially both ``pileup``
  and ``simu`` modules.

Note the simulation of this version has following features:
* it requires clone-specific and allele-specific CNV profiles as input.
* it assumes the copy numbers of target CNV regions are both 1 for two 
  alleles, hence the simulation is recommended to be performed on normal or
  reference cells.
* it only supports 10x scRNA-seq (UMI-based) data.
* when forking or discarding UMIs, especially the ambiguous UMIs, it is
  performed in (input-)region level instead of gene level.

Here, the ambiguous UMIs are those without allele information.

