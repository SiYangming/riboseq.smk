rule umi_tools_extract_se:
    input:
        fastq1="results/fastq_files/{sample}_cutadapt.fastq.gz",
    output:
        fastq="results/fastq_files/{sample}_UMI_clipped.fastq.gz",
    log:
        "logs/umi_tools/extract/{sample}.log",
    conda:
        "../envs/umi_tools.yaml"
    threads: THREADS
    params:
        extract_method=umi_extract_method,
        bc_pattern=umi_bc_pattern,
        extra=config["umi_tools"].get("extra_params", ""),
    script:
        "../scripts/umi_tools_extract.py"


rule umi_tools_extract_pe:
    input:
        fastq1="results/fastq_files/{sample}_R1_cutadapt.fastq.gz",
        fastq2="results/fastq_files/{sample}_R2_cutadapt.fastq.gz",
    output:
        fastq1="results/fastq_files/{sample}_R1_UMI_clipped.fastq.gz",
        fastq2="results/fastq_files/{sample}_R2_UMI_clipped.fastq.gz",
    log:
        "logs/umi_tools/extract/{sample}_pe.log",
    conda:
        "../envs/umi_tools.yaml"
    threads: THREADS
    params:
        extract_method=umi_extract_method,
        bc_pattern=umi_bc_pattern,
        extra=config["umi_tools"].get("extra_params", ""),
    script:
        "../scripts/umi_tools_extract.py"


rule umi_tools_dedup_pc:
    input:
        bam="results/BAM_files/{sample}_pc_sorted.bam",
        bai="results/BAM_files/{sample}_pc_sorted.bam.bai",
    output:
        bam="results/BAM_files/{sample}_pc_deduplicated.bam",
    log:
        "logs/umi_tools/dedup/{sample}_pc.log",
    conda:
        "../envs/umi_tools.yaml"
    threads: THREADS
    params:
        stats_prefix="logs/umi_tools/dedup/{sample}_pc",
        extra=lambda wc: (
            config["umi_tools"].get("extra_params", "")
            + (" --paired" if is_pe(wc.sample) else "")
        ).strip(),
    script:
        "../scripts/umi_tools_dedup.py"


rule umi_tools_dedup_genome:
    input:
        bam="results/BAM_files/{sample}_genome_sorted.bam",
        bai="results/BAM_files/{sample}_genome_sorted.bam.bai",
    output:
        bam="results/BAM_files/{sample}_genome_deduplicated.bam",
    log:
        "logs/umi_tools/dedup/{sample}_genome.log",
    conda:
        "../envs/umi_tools.yaml"
    threads: THREADS
    params:
        stats_prefix="logs/umi_tools/dedup/{sample}_genome",
        extra=lambda wc: (
            config["umi_tools"].get("extra_params", "")
            + (" --paired" if is_pe(wc.sample) else "")
        ).strip(),
    script:
        "../scripts/umi_tools_dedup.py"
