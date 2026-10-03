rule stage_fastq:
    input:
        lambda wc: fq_raw(wc.sample, wc.mate)
    output:
        "results/fastq_files/raw/{sample}_{mate}.fastq.gz"
    wildcard_constraints:
        mate="R1|R2"
    script:
        "../scripts/stage_fastq.py"


rule fastqc_raw:
    input:
        fastq="results/fastq_files/raw/{sample}_{mate}.fastq.gz"
    output:
        html="results/fastQC_files/raw/{sample}_{mate}_fastqc.html",
        zip="results/fastQC_files/raw/{sample}_{mate}_fastqc.zip"
    params:
        extra=config["fastqc"].get("extra", ""),
        mem_overhead_factor=config["fastqc"].get("mem_overhead_factor", 0.1)
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/raw/{sample}_{mate}.log"
    threads:
        THREADS
    resources:
        mem_mb=config["fastqc"].get("mem_mb", 8192)
    script:
        "../scripts/fastqc.py"


rule fastqc_file:
    input:
        fastq="results/fastq_files/{prefix}.fastq.gz"
    output:
        html="results/fastQC_files/{prefix}_fastqc.html",
        zip="results/fastQC_files/{prefix}_fastqc.zip"
    params:
        extra=config["fastqc"].get("extra", ""),
        mem_overhead_factor=config["fastqc"].get("mem_overhead_factor", 0.1)
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/{prefix}.log"
    threads:
        THREADS
    resources:
        mem_mb=config["fastqc"].get("mem_mb", 8192)
    script:
        "../scripts/fastqc.py"
