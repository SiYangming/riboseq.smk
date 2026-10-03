rule rsem_prepare_reference:
    input:
        reference_genome=config["pc_fasta"]
    output:
        seq=RSEM_INDEX + ".seq",
        grp=RSEM_INDEX + ".grp",
        ti=RSEM_INDEX + ".ti",
        bt2_1=RSEM_INDEX + ".1.bt2",
        bt2_2=RSEM_INDEX + ".2.bt2",
        bt2_3=RSEM_INDEX + ".3.bt2",
        bt2_4=RSEM_INDEX + ".4.bt2",
        bt2_rev1=RSEM_INDEX + ".rev.1.bt2",
        bt2_rev2=RSEM_INDEX + ".rev.2.bt2"
    params:
        extra=config["rsem"].get("prepare_reference_params", "--bowtie2"),
        gtf=""
    conda:
        "../envs/rsem.yaml"
    log:
        "logs/rsem/prepare_reference.log"
    threads:
        THREADS
    script:
        "../scripts/rsem_prepare_reference.py"


rule rsem_filter_pairs:
    input:
        bam=lambda wc: pc_final_bam(wc.sample)
    output:
        out=temp("results/BAM_files/{sample}_pc_proper_pairs.bam")
    params:
        extra="-f 2 -b",
        region=""
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/view_pairs/{sample}.log"
    threads:
        THREADS
    script:
        "../scripts/samtools_view.py"


rule rsem_name_sort:
    input:
        "results/BAM_files/{sample}_pc_proper_pairs.bam"
    output:
        temp("results/BAM_files/{sample}_pc_deduplicated_name_sorted.bam")
    params:
        extra="-n",
        mem_overhead_factor=0.1
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/namesort/{sample}.log"
    threads:
        THREADS
    resources:
        mem_mb=8192
    script:
        "../scripts/samtools_sort.py"


def rsem_alignments(wildcards):
    if is_pe(wildcards.sample):
        return "results/BAM_files/{sample}_pc_deduplicated_name_sorted.bam".format(
            sample=wildcards.sample
        )
    return pc_final_bam(wildcards.sample)


rule rsem_calculate_expression:
    input:
        bam=rsem_alignments,
        idx=multiext(RSEM_INDEX, ".seq", ".grp", ".ti")
    output:
        genes="results/rsem/{sample}.genes.results",
        isoforms="results/rsem/{sample}.isoforms.results"
    params:
        index=RSEM_INDEX,
        extra=config["rsem"].get("calculate_expression_params", ""),
        mean=config["rsem"].get("fragment_length_mean", 300),
        sd=config["rsem"].get("fragment_length_sd", 100),
        strandedness=config["rsem"].get("strandedness", "forward"),
        paired_end=lambda wc: is_pe(wc.sample)
    conda:
        "../envs/rsem.yaml"
    log:
        "logs/rsem/{sample}.log"
    threads:
        THREADS
    script:
        "../scripts/rsem_calculate_expression.py"
