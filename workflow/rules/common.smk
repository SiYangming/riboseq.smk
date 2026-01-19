rule all:
	input: []

# import basic packages
import os
import pandas as pd
from snakemake.utils import validate

# Config parsing
sample_sheet_path = config.get("samples", {}).get("sheet")
if sample_sheet_path is None:
    sample_sheet_path = config.get("sample_sheet")

if not sample_sheet_path:
    raise ValueError("Sample sheet path not found in configuration (samples.sheet or sample_sheet)")

ext = os.path.splitext(sample_sheet_path)[1].lower()
sep = "," if ext == ".csv" else "\t"

# read sample sheet
samples = (
	pd.read_csv(sample_sheet_path, sep=sep, dtype={"sample": str})
	.set_index("sample", drop=False)
	.sort_index()
)

# validate sample sheet and config file
validate(samples, schema="../schemas/samples.schema.yaml")
validate(config, schema="../schemas/config.schema.yaml")
