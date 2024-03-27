scCNASimulator
==============

scCNASimulator: CNA Simulation in droplet-based scRNA-seq data
--------------------------------------------------------------

The tool is a python package designed for CNA simulation in droplet-based 
scRNA-seq data.
It mainly takes a indexed BAM file and clone-specific CNA profile as input,
and output a new indexed BAM file containing the desired CNA alignments.

The tool has three modules: ``pileup``, ``simu``, and ``pipeline``.

* The ``pileup`` module pileups allele-specific UMIs in each single cell.
* The ``simu`` module simulates CNAs by forking or discarding UMIs (both the
  allele-specific and ambiguous UMIs) based on the given a clone-specific 
  CNA profile, and output a new indexed BAM file.
* The ``pipeline`` module is a wrapper that runs sequentially both ``pileup``
  and ``simu`` modules.

Currently, the simulation has following features:

* It requires clone-specific and allele-specific CNV profiles as input.
* It assumes the copy numbers of target CNV regions are both 1 for two 
  alleles, hence the simulation is recommended to be performed on normal or
  reference cells.
* It only supports droplet-based scRNA-seq data.
* When forking or discarding UMIs, especially the ambiguous UMIs, it is
  performed in (input-)region level instead of gene level.

Here, the ambiguous UMIs are those without allele information.


News
~~~~

Release notes are at `docs/release.rst <./docs/release.rst>`_.


Installation
------------

.. code-block:: bash

   pip install -U git+https://github.com/hxj5/scCNASimulator


Manual
------

The full manual is at `docs/manual.rst <./docs/manual.rst>`_.


FAQ and feedback
----------------

For troubleshooting, please have a look of `docs/FAQ.rst <./docs/FAQ.rst>`_, 
and we welcome reporting any issue_ for bugs, questions and new feature 
requests.


Acknowledgement
---------------


.. _issue: https://github.com/hxj5/scCNASimulator/issues

