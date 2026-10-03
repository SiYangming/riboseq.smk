def bbmap_filter_extra(wildcards):
    mem = config["bbmap"].get("memory", "Xmx=4g")
    extra = config["bbmap"].get("filter_extra", "ambiguous=best nodisk")
    return f"{mem} {extra}".strip()


def bbmap_pc_extra(wildcards):
    mem = config["bbmap"].get("memory", "Xmx=4g")
    extra = config["bbmap"].get(
        "pc_extra", "ambiguous=best nodisk trimreaddescription=t"
    )
    return f"{mem} {extra}".strip()


rule bbmap_rrna:
    input:
        fastq=lambda wc: umi_fq(wc.sample),
    output:
        mapped="results/fastq_files/{sample}_rRNA.fastq.gz",
        unmapped="results/fastq_files/{sample}_non_rRNA.fastq.gz",
    log:
        "logs/bbmap/{sample}_rRNA.log",
    conda:
        "../envs/bbmap.yaml"
    threads: THREADS
    params:
        command="bbmap.sh",
        extra=bbmap_filter_extra,
        ref=config["rrna_fasta"],
    script:
        "../scripts/bbtools.py"


rule bbmap_trna:
    input:
        fastq="results/fastq_files/{sample}_non_rRNA.fastq.gz",
    output:
        mapped="results/fastq_files/{sample}_tRNA.fastq.gz",
        unmapped="results/fastq_files/{sample}_non_rRNA_tRNA.fastq.gz",
    log:
        "logs/bbmap/{sample}_tRNA.log",
    conda:
        "../envs/bbmap.yaml"
    threads: THREADS
    params:
        command="bbmap.sh",
        extra=bbmap_filter_extra,
        ref=config["trna_fasta"],
    script:
        "../scripts/bbtools.py"


rule bbmap_pc:
    input:
        fastq="results/fastq_files/{sample}_non_rRNA_tRNA.fastq.gz",
        ref=lambda wc: most_abundant_fasta(),
    output:
        bam="results/BAM_files/{sample}_pc.bam",
        mapped="results/fastq_files/{sample}_pc.fastq.gz",
        unmapped="results/fastq_files/{sample}_unaligned.fastq.gz",
    log:
        "logs/bbmap/{sample}_pc.log",
    conda:
        "../envs/bbmap.yaml"
    threads: THREADS
    params:
        command="bbmap.sh",
        extra=bbmap_pc_extra,
        ref=lambda wc: most_abundant_fasta(),
    script:
        "../scripts/bbtools.py"
