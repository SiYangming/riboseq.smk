def umi_extract_method(wildcards):
    if sample_type(wildcards.sample) in RPF_TYPES:
        return config["rpf"].get("umi_extract_method", "regex")
    return config["totals"].get("umi_extract_method", "string")


def umi_bc_pattern(wildcards):
    if sample_type(wildcards.sample) in RPF_TYPES:
        return config["rpf"]["umi_bc_pattern"]
    return config["totals"]["umi_bc_pattern"]


rule umi_tools_extract_se:
    input:
        fastq1="results/fastq_files/{sample}_cutadapt.fastq.gz"
    output:
        fastq="results/fastq_files/{sample}_UMI_clipped.fastq.gz"
    params:
        extract_method=umi_extract_method,
        bc_pattern=umi_bc_pattern,
        extra=config["umi_tools"].get("extra_params", "")
    conda:
        "../envs/umi_tools.yaml"
    log:
        "logs/umi_tools/extract/{sample}.log"
    threads:
        THREADS
    script:
        "../scripts/umi_tools_extract.py"


rule umi_tools_extract_pe:
    input:
        fastq1="results/fastq_files/{sample}_R1_cutadapt.fastq.gz",
        fastq2="results/fastq_files/{sample}_R2_cutadapt.fastq.gz"
    output:
        fastq1="results/fastq_files/{sample}_R1_UMI_clipped.fastq.gz",
        fastq2="results/fastq_files/{sample}_R2_UMI_clipped.fastq.gz"
    params:
        extract_method=umi_extract_method,
        bc_pattern=umi_bc_pattern,
        extra=config["umi_tools"].get("extra_params", "")
    conda:
        "../envs/umi_tools.yaml"
    log:
        "logs/umi_tools/extract/{sample}_pe.log"
    threads:
        THREADS
    script:
        "../scripts/umi_tools_extract.py"


rule umi_tools_dedup_pc:
    input:
        bam="results/BAM_files/{sample}_pc_sorted.bam",
        bai="results/BAM_files/{sample}_pc_sorted.bam.bai"
    output:
        bam="results/BAM_files/{sample}_pc_deduplicated.bam"
    params:
        stats_prefix="logs/umi_tools/dedup/{sample}_pc",
        extra=lambda wc: (
            config["umi_tools"].get("extra_params", "")
            + (" --paired" if is_pe(wc.sample) else "")
        ).strip()
    conda:
        "../envs/umi_tools.yaml"
    log:
        "logs/umi_tools/dedup/{sample}_pc.log"
    threads:
        THREADS
    script:
        "../scripts/umi_tools_dedup.py"


rule umi_tools_dedup_genome:
    input:
        bam="results/BAM_files/{sample}_genome_sorted.bam",
        bai="results/BAM_files/{sample}_genome_sorted.bam.bai"
    output:
        bam="results/BAM_files/{sample}_genome_deduplicated.bam"
    params:
        stats_prefix="logs/umi_tools/dedup/{sample}_genome",
        extra=lambda wc: (
            config["umi_tools"].get("extra_params", "")
            + (" --paired" if is_pe(wc.sample) else "")
        ).strip()
    conda:
        "../envs/umi_tools.yaml"
    log:
        "logs/umi_tools/dedup/{sample}_genome.log"
    threads:
        THREADS
    script:
        "../scripts/umi_tools_dedup.py"
