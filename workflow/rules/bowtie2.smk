rule bowtie2_align_totals:
    input:
        sample=lambda wc: align_reads(wc.sample),
        idx=multiext(
            RSEM_INDEX,
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        )
    output:
        sam="results/SAM_files/{sample}_pc.sam"
    params:
        extra=config["bowtie2"].get(
            "bowtie2_extra_params",
            "--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1",
        )
    conda:
        "../envs/bowtie2.yaml"
    log:
        "logs/bowtie2/{sample}.log"
    threads:
        THREADS
    script:
        "../scripts/bowtie2_align.py"
