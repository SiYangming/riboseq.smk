rule cutadapt_se:
    input:
        fastq1="results/fastq_files/raw/{sample}_R1.fastq.gz",
    output:
        fastq="results/fastq_files/{sample}_cutadapt.fastq.gz",
        qc="results/fastq_files/{sample}_cutadapt.qc.txt",
    log:
        "logs/cutadapt/{sample}.log",
    conda:
        "../envs/cutadapt.yaml"
    threads: THREADS
    params:
        adapters=cutadapt_adapters,
        extra=cutadapt_extra,
    script:
        "../scripts/cutadapt.py"


rule cutadapt_pe:
    input:
        fastq1="results/fastq_files/raw/{sample}_R1.fastq.gz",
        fastq2="results/fastq_files/raw/{sample}_R2.fastq.gz",
    output:
        fastq1="results/fastq_files/{sample}_R1_cutadapt.fastq.gz",
        fastq2="results/fastq_files/{sample}_R2_cutadapt.fastq.gz",
        qc="results/fastq_files/{sample}_pe_cutadapt.qc.txt",
    log:
        "logs/cutadapt/{sample}_pe.log",
    conda:
        "../envs/cutadapt.yaml"
    threads: THREADS
    params:
        adapters=cutadapt_adapters,
        extra=cutadapt_extra,
    script:
        "../scripts/cutadapt.py"
