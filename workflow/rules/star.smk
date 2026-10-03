rule star_index:
    input:
        unpack(star_index_input),
    output:
        directory(STAR_INDEX),
    log:
        "logs/star/index.log",
    conda:
        "../envs/star.yaml"
    threads: THREADS
    params:
        extra=config["star"].get("index_extra", ""),
        sjdbOverhang=config["star"].get("sjdb_overhang", 100),
    script:
        "../scripts/star_index.py"


rule star_align:
    input:
        unpack(star_align_input),
    output:
        aln="results/STAR/{sample}_Aligned.out.bam",
        sj="results/STAR/{sample}_SJ.out.tab",
        log_out="results/STAR/{sample}.Log.out",
        log_final="results/STAR/{sample}.Log.final.out",
    log:
        "logs/star/align/{sample}.log",
    conda:
        "../envs/star.yaml"
    threads: THREADS
    params:
        extra=config["star"].get(
            "align_extra",
            "--outSAMtype BAM Unsorted --outSAMprimaryFlag AllBestScore "
            "--outFilterMultimapNmax 5 --outFilterMismatchNmax 5 --alignEndsType EndToEnd",
        ),
    script:
        "../scripts/star_align.py"
