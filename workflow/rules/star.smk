def star_index_input(_wildcards=None):
    files = {"fasta": config.get("genome_fasta") or config["pc_fasta"]}
    gtf = config.get("gtf") or ""
    if gtf:
        files["gtf"] = gtf
    return files


rule star_index:
    input:
        unpack(star_index_input)
    output:
        directory(STAR_INDEX)
    params:
        extra=config["star"].get("index_extra", ""),
        sjdbOverhang=config["star"].get("sjdb_overhang", 100)
    conda:
        "../envs/star.yaml"
    log:
        "logs/star/index.log"
    threads:
        THREADS
    script:
        "../scripts/star_index.py"


def star_align_input(wildcards):
    reads = {"idx": STAR_INDEX}
    if is_pe(wildcards.sample):
        reads["fq1"] = umi_fq(wildcards.sample, "R1")
        reads["fq2"] = umi_fq(wildcards.sample, "R2")
    else:
        reads["fq1"] = umi_fq(wildcards.sample)
    return reads


rule star_align:
    input:
        unpack(star_align_input)
    output:
        aln="results/STAR/{sample}_Aligned.out.bam",
        sj="results/STAR/{sample}_SJ.out.tab",
        log_out="results/STAR/{sample}.Log.out",
        log_final="results/STAR/{sample}.Log.final.out"
    params:
        extra=config["star"].get(
            "align_extra",
            "--outSAMtype BAM Unsorted --outSAMprimaryFlag AllBestScore "
            "--outFilterMultimapNmax 5 --outFilterMismatchNmax 5 --alignEndsType EndToEnd",
        )
    conda:
        "../envs/star.yaml"
    log:
        "logs/star/align/{sample}.log"
    threads:
        THREADS
    script:
        "../scripts/star_align.py"
