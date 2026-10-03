rule faidx_pc:
    input:
        fasta=most_abundant_fasta(),
    output:
        most_abundant_fasta() + ".fai",
    log:
        "logs/samtools/faidx_pc.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input.fasta} > {log} 2>&1"


rule count_by_length:
    input:
        bam=lambda wc: pc_final_bam(wc.sample),
        bai=lambda wc: pc_final_bam(wc.sample) + ".bai",
        fasta=lambda wc: most_abundant_fasta(),
        fai=lambda wc: most_abundant_fasta() + ".fai",
    output:
        "results/Counts_files/{sample}_pc_L{length}_Off0.counts",
    log:
        "logs/counts/{sample}_L{length}.log",
    wildcard_constraints:
        length="[0-9]+",
    conda:
        "../envs/python.yaml"
    params:
        sdir=SCRIPT_DIR,
        outdir=lambda wildcards, output: os.path.dirname(str(output[0])),
    shell:
        "PYTHONPATH={params.sdir} python {params.sdir}/counting_script.py "
        "-bam {input.bam} -fasta {input.fasta} -len {wildcards.length} -offset 0 "
        "-out_file {wildcards.sample}_pc_L{wildcards.length}_Off0.counts "
        "-out_dir {params.outdir} > {log} 2>&1"


rule sum_region_counts:
    input:
        counts="results/Counts_files/{sample}_pc_L{length}_Off0.counts",
        regions=config["region_lengths"],
    output:
        "results/Analysis/region_counts/{sample}_pc_L{length}_Off0_region_counts.csv",
    log:
        "logs/counts/region/{sample}_L{length}.log",
    conda:
        "../envs/python.yaml"
    params:
        sdir=SCRIPT_DIR,
        offset=config.get("region_offset", 15),
    shell:
        "PYTHONPATH={params.sdir} python {params.sdir}/summing_region_counts.py "
        "{wildcards.sample}_pc_L{wildcards.length}_Off0.counts {params.offset} {input.regions} "
        "-in_dir results/Counts_files -out_dir results/Analysis/region_counts > {log} 2>&1"


rule sum_spliced_counts:
    input:
        counts="results/Counts_files/{sample}_pc_L{length}_Off0.counts",
        regions=config["region_lengths"],
    output:
        "results/Analysis/spliced_counts/{sample}_pc_L{length}_Off0_start_site.csv",
    log:
        "logs/counts/spliced/{sample}_L{length}.log",
    conda:
        "../envs/python.yaml"
    params:
        sdir=SCRIPT_DIR,
        n=config.get("splice_nt", 50),
    shell:
        "PYTHONPATH={params.sdir} python {params.sdir}/summing_spliced_counts.py "
        "{wildcards.sample}_pc_L{wildcards.length}_Off0.counts {params.n} {input.regions} "
        "-in_dir results/Counts_files -out_dir results/Analysis/spliced_counts > {log} 2>&1"


rule periodicity:
    input:
        counts="results/Counts_files/{sample}_pc_L{length}_Off0.counts",
        regions=config["region_lengths"],
    output:
        "results/Analysis/periodicity/{sample}_pc_L{length}_Off0_periodicity.csv",
    log:
        "logs/counts/periodicity/{sample}_L{length}.log",
    conda:
        "../envs/python.yaml"
    params:
        sdir=SCRIPT_DIR,
        offset=config.get("region_offset", 15),
    shell:
        "PYTHONPATH={params.sdir} python {params.sdir}/periodicity.py "
        "{wildcards.sample}_pc_L{wildcards.length}_Off0.counts {input.regions} "
        "-offset {params.offset} -in_dir results/Counts_files "
        "-out_dir results/Analysis/periodicity > {log} 2>&1"


rule count_final:
    input:
        bam=lambda wc: pc_final_bam(wc.sample),
        bai=lambda wc: pc_final_bam(wc.sample) + ".bai",
        fasta=lambda wc: most_abundant_fasta(),
        fai=lambda wc: most_abundant_fasta() + ".fai",
    output:
        "results/Counts_files/{sample}_pc_final.counts",
    log:
        "logs/counts/{sample}_final.log",
    conda:
        "../envs/python.yaml"
    params:
        sdir=SCRIPT_DIR,
        lengths=config.get("final_lengths", "27,28,29,30,31,32"),
        offsets=config.get("final_offsets", "12,12,12,12,13,13"),
        outdir=lambda wildcards, output: os.path.dirname(str(output[0])),
    shell:
        "PYTHONPATH={params.sdir} python {params.sdir}/counting_script.py "
        "-bam {input.bam} -fasta {input.fasta} -len {params.lengths} -offset {params.offsets} "
        "-out_file {wildcards.sample}_pc_final.counts -out_dir {params.outdir} > {log} 2>&1"


rule sum_cds_counts:
    input:
        counts="results/Counts_files/{sample}_pc_final.counts",
        regions=config["region_lengths"],
    output:
        "results/Analysis/CDS_counts/{sample}_pc_final_counts_all_frames.csv",
    log:
        "logs/counts/cds/{sample}.log",
    conda:
        "../envs/python.yaml"
    params:
        sdir=SCRIPT_DIR,
    shell:
        "PYTHONPATH={params.sdir} python {params.sdir}/summing_CDS_counts.py "
        "{wildcards.sample}_pc_final.counts {input.regions} -remove_end_codons "
        "-in_dir results/Counts_files -out_dir results/Analysis/CDS_counts > {log} 2>&1"


rule sum_utr5_counts:
    input:
        counts="results/Counts_files/{sample}_pc_final.counts",
        regions=config["region_lengths"],
    output:
        "results/Analysis/UTR5_counts/{sample}_pc_final_UTR5_counts.csv",
    log:
        "logs/counts/utr5/{sample}.log",
    conda:
        "../envs/python.yaml"
    params:
        sdir=SCRIPT_DIR,
    shell:
        "PYTHONPATH={params.sdir} python {params.sdir}/summing_UTR5_counts.py "
        "{wildcards.sample}_pc_final.counts {input.regions} "
        "-in_dir results/Counts_files -out_dir results/Analysis/UTR5_counts > {log} 2>&1"


rule counts_to_csv:
    input:
        counts="results/Counts_files/{sample}_pc_final.counts",
        fasta=lambda wc: most_abundant_fasta(),
    output:
        "results/Counts_files/csv_files/{sample}_pc_final_counts.csv",
    log:
        "logs/counts/csv/{sample}.log",
    conda:
        "../envs/python.yaml"
    params:
        sdir=SCRIPT_DIR,
    shell:
        "PYTHONPATH={params.sdir} python {params.sdir}/counts_to_csv.py "
        "{wildcards.sample}_pc_final.counts {input.fasta} -one_csv "
        "-in_dir results/Counts_files -out_dir results/Counts_files/csv_files > {log} 2>&1"


rule count_codon_occupancy:
    input:
        counts="results/Counts_files/{sample}_pc_final.counts",
        fasta=lambda wc: most_abundant_fasta(),
        regions=config["region_lengths"],
    output:
        "results/Analysis/codon_counts/{sample}_pc_final_codon_counts_20_-10.csv",
    log:
        "logs/counts/codon/{sample}.log",
    conda:
        "../envs/python.yaml"
    params:
        sdir=SCRIPT_DIR,
    shell:
        "PYTHONPATH={params.sdir} python {params.sdir}/count_codon_occupancy.py "
        "{wildcards.sample}_pc_final.counts {input.fasta} {input.regions} 20 -10 "
        "-in_dir results/Counts_files -out_dir results/Analysis/codon_counts > {log} 2>&1"


rule calculate_most_abundant:
    input:
        expand("results/rsem/{sample}.isoforms.results", sample=totals_samples()),
    output:
        txt="results/Analysis/most_abundant_transcripts/most_abundant_transcripts.txt",
        csv="results/Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv",
    log:
        "logs/R/calculate_most_abundant_transcript.log",
    conda:
        "../envs/r_analysis.yaml"
    params:
        **r_params(),
        parent_dir=results_parent,
        r_script="calculate_most_abundant_transcript.R",
    script:
        "../scripts/run_r.py"


rule filter_most_abundant_fasta:
    input:
        fasta=config["pc_fasta"],
        ids="results/Analysis/most_abundant_transcripts/most_abundant_transcripts.txt",
    output:
        "results/Analysis/most_abundant_transcripts/most_abundant_transcripts.fa",
    log:
        "logs/counts/filter_fasta.log",
    conda:
        "../envs/python.yaml"
    params:
        sdir=SCRIPT_DIR,
    shell:
        "PYTHONPATH={params.sdir} python {params.sdir}/filter_FASTA.py "
        "{input.fasta} {input.ids} {output} > {log} 2>&1"
