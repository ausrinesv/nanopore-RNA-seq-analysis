
# Pipeline for processing data from direct nanopore RNA sequencing data

rule all:
    input:
        "guppy_out",
        "merged/merged.fastq",
        "merged/merged.bam",
        "merged/merged.blow5",
        "merged/merged.bam.bai",
        "eventalign.txt",
        "m6anet_out/dataprep",
        "m6anet_out/predictions",
        "counts.tsv",
        "m6anet_out/predictions"

rule guppy_basecalling:
    input:
        fast5_dir = "fast5"
    output:
        out_dir = directory("guppy_out")
    params:
        flowcell = "FLO-FLG001",
        kit = "SQK-RNA002",
        guppy = "tools/ont-guppy-cpu/bin/guppy_basecaller",
        ref = "annotations/Saccharomyces_cerevisiae_YHR174W_mRNA_sequence.fa"
    shell:
        """
        {params.guppy} -i {input.fast5_dir} -s {output.out_dir} --flowcell {params.flowcell} --kit {params.kit} --bam_out --align_ref {params.ref} --minimap_opt_string -p0 --minimap_opt_string -N10
        """

rule merge_fastq_files:
    input:
        fastq_dir = "guppy_out"
    output:
        merged_fastq = "merged/merged.fastq"
    shell:
        """
        mkdir -p merged
        cat {input.fastq_dir}/pass/*.fastq > {output.merged_fastq}
        """

rule merge_bam_files:
    input:
        bam_dir = "guppy_out"
    output:
        merged_bam = "merged/merged.bam"
    shell:
        """
        samtools merge {output.merged_bam} {input.bam_dir}/pass/*.bam
        """

rule fast5_to_blow5:
    input:
        fast5_dir = "fast5"
    output:
        blow5_dir = directory("blow5")
    params:
        slow5tools = "tools/slow5tools-v1.1.0/slow5tools"
    shell:
        """
        {params.slow5tools} f2s {input.fast5_dir} -d {output.blow5_dir}
        """

rule merge_slow5_files:
    input:
        blow5_dir = "blow5"
    output:
        merged_file = "merged/merged.blow5"
    params:
        slow5tools = "tools/slow5tools-v1.1.0/slow5tools"
    shell:
        """
        {params.slow5tools} merge {input.blow5_dir} -o {output.merged_file}
        """

rule create_index:
    input:
        merged_blow5 = "merged/merged.blow5",
        merged_bam = "merged/merged.bam"
    output:
        "merged/merged.blow5.idx",
        "merged/merged.fastq.index",
        "merged/merged.fastq.index.fai",
        "merged/merged.fastq.index.gzi",
        "merged/merged.fastq.index.readdb",
        "merged/merged.bam.bai"
    params:
        nanopolish = "tools/nanopolish/nanopolish"
    shell:
        """
        {params.nanopolish} index {input.merged_fastq} --slow5 {input.merged_blow5}
        samtools index {input.merged_bam}
        """

rule eventalign:
    input:
        fastq = "merged/merged.fastq",
        bam = "merged/merged.bam",
        bam_idx = "merged/merged.bam.bai"
    params:
        nanopolish = "tools/nanopolish/nanopolish",
        ref = "annotations/Saccharomyces_cerevisiae_YHR174W_mRNA_sequence.fa"
    output:
        "eventalign.txt"
    shell:
        """
        {params.nanopolish} eventalign --reads {input.fastq} --bam {input.bam} --genome {params.ref} --scale-events --signal-index --summary eventalign_summary.txt > {output}
        """

rule m6anet_dataprep:
    input:
        "eventalign.txt"
    output:
        directory("m6anet_out/dataprep")
#        "m6anet_out/data.json",
#        "m6anet_out/data.log",
#        "m6anet_out/data.info",
#        "m6anet_out/eventalign.index"
    shell:
        """
        m6anet dataprep --eventalign {input} --out_dir {output} --n_processes 2
        """
        

rule m6anet_inference:
    input:
        "m6anet_out/dataprep"
    output:
        directory("m6anet_out/predictions")
#        "m6anet_out/data.indiv_proba.csv",
#        "m6anet_out/data.site_proba.csv"
    shell:
        """
        m6anet inference --input_dir {input} --out_dir {output} --n_processes 2 --num_iterations 1000
        """

rule estimate_counts:
    input:
        "merged/merged.bam"
    output:
        "counts.tsv"
    shell:
        """
        NanoCount -i {input} -o {output}
        """
