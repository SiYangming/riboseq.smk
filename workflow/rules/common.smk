import os

import pandas as pd
from snakemake.utils import validate

samples = (
    pd.read_csv(config["sample_sheet"], sep="\t", dtype=str)
    .set_index("sample", drop=False)
    .sort_index()
)
validate(samples, schema="../schemas/samples.schema.yaml")
validate(config, schema="../schemas/config.schema.yaml")

ENTRY = config["entrypoint"]
THREADS = int(config.get("threads", 4))
SKIP_UMI = bool(config.get("skip_umi", False))
USE_STAR = bool(config.get("use_star", True))
RPF_TYPES = {"riboseq", "rpf", "rpfs"}
TOTAL_TYPES = {"rnaseq", "totals", "total", "rna"}
RPF_LENGTHS = list(
    range(int(config.get("rpf_min_len", 25)), int(config.get("rpf_max_len", 35)) + 1)
)
RSEM_INDEX = "results/index/rsem/rsem"
STAR_INDEX = "results/index/star"
SCRIPT_DIR = os.path.join(str(workflow.basedir), "scripts")


def _cell(sample, col):
    if col not in samples.columns:
        return ""
    val = samples.loc[sample, col]
    if pd.isna(val):
        return ""
    text = str(val).strip()
    return "" if text in {".", "nan", "None"} else text


def is_pe(sample):
    return bool(_cell(sample, "fastq_2"))


def sample_type(sample):
    return str(_cell(sample, "type") or "riboseq").lower()


def rpf_samples():
    return [s for s in samples["sample"] if sample_type(s) in RPF_TYPES]


def totals_samples():
    return [s for s in samples["sample"] if sample_type(s) in TOTAL_TYPES]


def mates(sample):
    return ["R1", "R2"] if is_pe(sample) else ["R1"]


def fq_raw(sample, mate):
    if mate == "R2":
        return _cell(sample, "fastq_2")
    return _cell(sample, "fastq_1")


def staged_fq(sample, mate):
    return f"results/fastq_files/raw/{sample}_{mate}.fastq.gz"


def cutadapt_fq(sample, mate=None):
    if is_pe(sample):
        mate = mate or "R1"
        return f"results/fastq_files/{sample}_{mate}_cutadapt.fastq.gz"
    return f"results/fastq_files/{sample}_cutadapt.fastq.gz"


def umi_fq(sample, mate=None):
    if SKIP_UMI:
        return cutadapt_fq(sample, mate)
    if is_pe(sample):
        mate = mate or "R1"
        return f"results/fastq_files/{sample}_{mate}_UMI_clipped.fastq.gz"
    return f"results/fastq_files/{sample}_UMI_clipped.fastq.gz"


def align_reads(sample):
    if is_pe(sample):
        return [umi_fq(sample, "R1"), umi_fq(sample, "R2")]
    return [umi_fq(sample)]


def pc_sorted_bam(sample):
    return f"results/BAM_files/{sample}_pc_sorted.bam"


def pc_final_bam(sample):
    if SKIP_UMI:
        return pc_sorted_bam(sample)
    return f"results/BAM_files/{sample}_pc_deduplicated_sorted.bam"


def most_abundant_fasta():
    configured = config.get("most_abundant_fasta") or ""
    if configured:
        return configured
    if totals_samples():
        return "results/Analysis/most_abundant_transcripts/most_abundant_transcripts.fa"
    return config["pc_fasta"]


def r_params():
    return {
        "rdir": os.path.join(str(workflow.basedir), "scripts", "R"),
        "project_root": os.path.abspath("."),
        "info_csv": config.get("info_csv", ""),
        "fasta_dir": config.get("fasta_dir", "reference"),
        "rpf_names": " ".join(rpf_samples()),
        "totals_names": " ".join(totals_samples()),
    }


def results_parent(wildcards, output):
    return "results"


def pipeline_targets():
    targets = []
    for sample in samples["sample"]:
        for mate in mates(sample):
            targets.append(f"results/fastQC_files/raw/{sample}_{mate}_fastqc.html")
    if ENTRY == "qc":
        return targets
    if ENTRY in ("totals", "all"):
        for sample in totals_samples():
            targets.append(f"results/rsem/{sample}.isoforms.results")
            if USE_STAR:
                targets.append(f"results/BAM_files/{sample}_genome_sorted.bam")
        if totals_samples():
            targets.append(
                "results/Analysis/most_abundant_transcripts/most_abundant_transcripts.fa"
            )
    if ENTRY in ("rpfs", "all"):
        for sample in rpf_samples():
            targets.append(
                f"results/Counts_files/csv_files/{sample}_pc_final_counts.csv"
            )
            targets.append(
                f"results/Analysis/codon_counts/{sample}_pc_final_codon_counts_20_-10.csv"
            )
    if ENTRY in ("downstream", "all"):
        targets.append("results/Analysis/downstream.done")
    return targets


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


def umi_extract_method(wildcards):
    if sample_type(wildcards.sample) in RPF_TYPES:
        return config["rpf"].get("umi_extract_method", "regex")
    return config["totals"].get("umi_extract_method", "string")


def umi_bc_pattern(wildcards):
    if sample_type(wildcards.sample) in RPF_TYPES:
        return config["rpf"]["umi_bc_pattern"]
    return config["totals"]["umi_bc_pattern"]


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


def rsem_alignments(wildcards):
    if is_pe(wildcards.sample):
        return f"results/BAM_files/{wildcards.sample}_pc_deduplicated_name_sorted.bam"
    return pc_final_bam(wildcards.sample)


def star_index_input(_wildcards=None):
    files = {"fasta": config.get("genome_fasta") or config["pc_fasta"]}
    gtf = config.get("gtf") or ""
    if gtf:
        files["gtf"] = gtf
    return files


def star_align_input(wildcards):
    reads = {"idx": STAR_INDEX}
    if is_pe(wildcards.sample):
        reads["fq1"] = umi_fq(wildcards.sample, "R1")
        reads["fq2"] = umi_fq(wildcards.sample, "R2")
    else:
        reads["fq1"] = umi_fq(wildcards.sample)
    return reads
