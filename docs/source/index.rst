Welcome to sequencing_tools's documentation!
============================================

This is `Douglas's repo <https://github.com/wckdouglas/sequencing_tools>`_ for useful bioinformatics stuff. 

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Installation
============

Installation can be done via conda::

    $ conda config --add channels bioconda
    $ conda create -q -n test-environment python=3.6 \
        cython numpy networkx seaborn pyBigwig six pysam ujson pytest scipy   \
        matplotlib samtools future pytest-cov codecov
    $ git clone https://github.com/wckdouglas/sequencing_tools.git
    $ cd sequencing_tools
    $ pip install .

fastq_tools
===========
This module contains functions to interact with fastq files

.. automodule:: sequencing_tools.fastq_tools._fastq_tools
    :members:

.. autoclass:: sequencing_tools.fastq_tools.pe_align.ConsensusBuilder
    :members:


fasta_tools
===========
This module contains functions to interact with fasta files

.. automodule:: sequencing_tools.fasta_tools
    :members:

.. autoclass:: sequencing_tools.fasta_tools.IUPAC
    :members:
    
.. autoclass:: sequencing_tools.fasta_tools.MultiAlignments
    :members:


bam_tools
=========
This module contains functions to interact with bam files and bed files

.. automodule:: sequencing_tools.bam_tools._bam_tools
    :members:

consensus_tools
===============

.. autoclass:: sequencing_tools.consensus_tools._consensus_tools.ErrorCorrection
    :members:


gene_tools
==========
This module contains functions to interact with bed files, refflat, and gtf files.

.. autoclass:: sequencing_tools.gene_tools.Bed12Record
    :members:

.. autoclass:: sequencing_tools.gene_tools.GTFRecord
    :members:

.. autoclass:: sequencing_tools.gene_tools.transcriptome.Transcriptome
    :members:

.. autoclass:: sequencing_tools.gene_tools.transcriptome.Transcript
    :members:

.. autoclass:: sequencing_tools.gene_tools.transcriptome.Exon
    :members:

viz_tools
=========
This module contains some functions that work with figures

.. automodule:: sequencing_tools.viz_tools
    :members:

stats_tools
===========
This module contains some math functions

.. automodule:: sequencing_tools.stats_tools._stats_tools
    :members:

.. autoclass:: sequencing_tools.stats_tools._stats_tools.Bootstrap
    :members:

.. autoclass:: sequencing_tools.stats_tools.regression.GradientDescent
    :members:


utils
=====

.. automodule:: sequencing_tools.utils
    :members:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
