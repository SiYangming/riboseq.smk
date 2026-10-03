rule samtools_sam_to_bam:
    input:
        sam="results/SAM_files/{sample}_pc.sam"
    output:
        bam=temp("results/BAM_files/{sample}_pc_unsorted.bam")
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/view/{sample}_pc.log"
    threads:
        THREADS
    script:
        "../scripts/samtools_sam_to_bam.py"


rule samtools_sort_pc:
    input:
        lambda wc: (
            f"results/BAM_files/{wc.sample}_pc.bam"
            if sample_type(wc.sample) in RPF_TYPES
            else f"results/BAM_files/{wc.sample}_pc_unsorted.bam"
        )
    output:
        "results/BAM_files/{sample}_pc_sorted.bam"
    params:
        extra=config["samtools"].get("sort_extra_params", ""),
        mem_overhead_factor=0.1
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/sort/{sample}_pc.log"
    threads:
        THREADS
    resources:
        mem_mb=8192
    script:
        "../scripts/samtools_sort.py"


rule samtools_index_pc:
    input:
        bam="results/BAM_files/{sample}_pc_sorted.bam"
    output:
        bai="results/BAM_files/{sample}_pc_sorted.bam.bai"
    params:
        extra=config["samtools"].get("index_extra_params", "")
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/index/{sample}_pc.log"
    threads:
        THREADS
    script:
        "../scripts/samtools_index.py"


rule samtools_sort_pc_dedup:
    input:
        "results/BAM_files/{sample}_pc_deduplicated.bam"
    output:
        "results/BAM_files/{sample}_pc_deduplicated_sorted.bam"
    params:
        extra=config["samtools"].get("sort_extra_params", ""),
        mem_overhead_factor=0.1
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/sort/{sample}_pc_dedup.log"
    threads:
        THREADS
    resources:
        mem_mb=8192
    script:
        "../scripts/samtools_sort.py"


rule samtools_index_pc_dedup:
    input:
        bam="results/BAM_files/{sample}_pc_deduplicated_sorted.bam"
    output:
        bai="results/BAM_files/{sample}_pc_deduplicated_sorted.bam.bai"
    params:
        extra=config["samtools"].get("index_extra_params", "")
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/index/{sample}_pc_dedup.log"
    threads:
        THREADS
    script:
        "../scripts/samtools_index.py"


rule samtools_sort_genome:
    input:
        "results/STAR/{sample}_Aligned.out.bam"
    output:
        "results/BAM_files/{sample}_genome_sorted.bam"
    params:
        extra=config["samtools"].get("sort_extra_params", ""),
        mem_overhead_factor=0.1
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/sort/{sample}_genome.log"
    threads:
        THREADS
    resources:
        mem_mb=8192
    script:
        "../scripts/samtools_sort.py"


rule samtools_index_genome:
    input:
        bam="results/BAM_files/{sample}_genome_sorted.bam"
    output:
        bai="results/BAM_files/{sample}_genome_sorted.bam.bai"
    params:
        extra=config["samtools"].get("index_extra_params", "")
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/index/{sample}_genome.log"
    threads:
        THREADS
    script:
        "../scripts/samtools_index.py"
