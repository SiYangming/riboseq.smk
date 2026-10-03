"""Copy or gzip-stage a FASTQ into the workflow tree."""

__author__ = "Yangming Si"
__copyright__ = "Copyright 2026, Yangming Si"
__email__ = "siyangming1991@163.com"
__license__ = "MIT"

import gzip
import os
import shutil

src = str(snakemake.input[0])
dst = str(snakemake.output[0])
log_path = str(snakemake.log[0]) if snakemake.log else None
os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
if log_path:
    os.makedirs(os.path.dirname(log_path) or ".", exist_ok=True)

if src.endswith(".gz") and dst.endswith(".gz"):
    try:
        os.link(os.path.abspath(src), dst)
        action = "hardlink"
    except OSError:
        shutil.copy2(src, dst)
        action = "copy"
elif dst.endswith(".gz"):
    with open(src, "rb") as inf, gzip.open(dst, "wb") as outf:
        shutil.copyfileobj(inf, outf)
    action = "gzip"
else:
    shutil.copy2(src, dst)
    action = "copy"

if log_path:
    with open(log_path, "w", encoding="utf-8") as logf:
        logf.write(f"{action}: {src} -> {dst}\n")
