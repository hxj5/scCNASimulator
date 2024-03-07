
..
   History
   =======

Release v0.0.2 (07/03/2024)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Bug fix:

* simu: no allele-specific UMIs extracted. It is because (1) the region IDs
  are not properly loaded into ``set()``; (2) the region IDs are not used by
  ``CNVRegCN`` as region names, hence the allele information cannot be 
  extracted by ``query()`` functions.

Feature enhancement:

* simu: add ``UMICount`` to count the processed allele-specific UMIs in each
  clonal CNV region.
* simu: add a ``--debug`` option.
* simu: remove the leading and tailing ``"`` when loading the region IDs and
  clone IDs.


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

