"""Snakemake wrapper for bowtie2 index."""

__author__ = "Yangming Si"
__copyright__ = "Copyright 2026, Yangming Si"
__email__ = "siyangming1991@163.com"
__license__ = "MIT"


import os
from snakemake.shell import shell
import os, sys
import exec_wrapper

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Docker wrapper integration
docker_prefix, tool_bin = exec_wrapper.exec_wrapper_binary(
    snakemake.config,
    "bowtie2",
    "index_bin",
    "bowtie2-build",
)

index = os.path.commonprefix(snakemake.output).rstrip(".")


shell(
    "{docker_prefix}{tool_bin}"
    " --threads {snakemake.threads}"
    " {extra}"
    " {snakemake.input.ref}"
    " {index}"
    " {log}"
)
