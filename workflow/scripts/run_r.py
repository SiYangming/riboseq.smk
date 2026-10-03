"""Run one or more R scripts from workflow/scripts/R with Ribo-seq env vars."""

__author__ = "Yangming Si"
__copyright__ = "Copyright 2026, Yangming Si"
__email__ = "siyangming1991@163.com"
__license__ = "MIT"

import os

from snakemake.shell import shell

rdir = snakemake.params.rdir
parent = snakemake.params.parent_dir
project = snakemake.params.project_root
info_csv = snakemake.params.get("info_csv", "")
fasta_dir = snakemake.params.get("fasta_dir", "")
rpf_names = snakemake.params.get("rpf_names", "")
totals_names = snakemake.params.get("totals_names", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

scripts = snakemake.params.get("scripts")
if scripts:
    script_list = list(scripts)
else:
    script_list = [snakemake.params.r_script]

os.makedirs(parent, exist_ok=True)
os.makedirs(os.path.join(parent, "Analysis", "DESeq2_output"), exist_ok=True)
os.makedirs(os.path.join(parent, "Analysis", "most_abundant_transcripts"), exist_ok=True)

joined = " && ".join(f"Rscript {s}" for s in script_list)
env_prefix = (
    "cd {rdir:q} && "
    "export RIBO_SEQ_PARENT_DIR={parent:q} && "
    "export RIBO_SEQ_PROJECT_ROOT={project:q} && "
    "export RIBO_SEQ_INFO_CSV={info_csv:q} && "
    "export RIBO_SEQ_FASTA_DIR={fasta_dir:q} && "
    "export RIBO_SEQ_RPF_FILENAMES={rpf_names:q} && "
    "export RIBO_SEQ_TOTALS_FILENAMES={totals_names:q} && "
)
shell(env_prefix + joined + " {log}")
