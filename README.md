This is the workflow for analyzing direct Oxford Nanopore RNA sequencing data

To run the workflow, the following directories should exist:
- fast5: directory that contains raw FAST5 files
- annotations: directory that contains reference transcriptome

Workflow can be run using snakemake -c{number of cores} command, e.g., snakemake -c1
