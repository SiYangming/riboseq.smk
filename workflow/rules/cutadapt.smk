def cutadapt_adapters(wildcards):
    if sample_type(wildcards.sample) in RPF_TYPES:
        return f"-a {config['rpf']['adaptor']}"
    adaptor = config["totals"]["adaptor"]
    if is_pe(wildcards.sample):
        return f"-a {adaptor} -A {adaptor}"
    return f"-a {adaptor}"


def cutadapt_extra(wildcards):
    extra = config["cutadapt"].get("extra", "") or "--nextseq-trim=20"
    if "--nextseq-trim" not in extra:
        extra = f"--nextseq-trim=20 {extra}".strip()
    if sample_type(wildcards.sample) in RPF_TYPES:
        return f"{extra} -m {config['rpf']['min_len']} -M {config['rpf']['max_len']}"
    return f"{extra} -m {config['totals']['min_len']}"


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
