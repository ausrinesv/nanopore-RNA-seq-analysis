Workflow for analyzing direct Oxford Nanopore RNA sequencing data

To run the workflow, the following directories should exist:
- tools: contains necessary tools for the analysis
- fast5: contains raw FAST5 files
- annotations: contains reference transcriptome

Workflow can be run using snakemake -c{number of cores} command, e.g., snakemake -c1
