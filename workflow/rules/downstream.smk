DOWNSTREAM_R = [
    "QC/RPFs_read_counts.R",
    "QC/periodicity.R",
    "QC/region_counts.R",
    "QC/heatmaps.R",
    "QC/offset_plots.R",
    "QC/offset_aligned_single_nt_plots.R",
    "DESeq2/DESeq2_Totals.R",
    "DESeq2/DESeq2_RPFs.R",
    "DESeq2/DESeq2_TE.R",
    "DESeq2/DE_analysis_plots.R",
    "meta_plots/bin_data.R",
    "meta_plots/plot_binned_data.R",
    "meta_plots/normalise_individual_transcript_counts.R",
    "meta_plots/plot_individual_mRNAs.R",
    "codon_occupancy.R",
    "gsea/fgsea.R",
]


rule downstream:
    input:
        expand(
            "results/Counts_files/csv_files/{sample}_pc_final_counts.csv",
            sample=rpf_samples(),
        ),
        expand("results/rsem/{sample}.isoforms.results", sample=totals_samples()),
    output:
        touch("results/Analysis/downstream.done"),
    log:
        "logs/R/downstream.log",
    conda:
        "../envs/r_analysis.yaml"
    params:
        **r_params(),
        parent_dir=results_parent,
        scripts=DOWNSTREAM_R,
    script:
        "../scripts/run_r.py"
