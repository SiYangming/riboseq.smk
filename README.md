# Snakemake workflow: `riboseq.smk`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/SiYangming/riboseq.smk/workflows/Tests/badge.svg?branch=main)](https://github.com/SiYangming/riboseq.smk/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

Ribosome profiling (RPFs) plus matched total RNA-seq (Totals) and Downstream R analysis, matching [SiYangming/Ribo-seq](https://github.com/SiYangming/Ribo-seq).

## Usage

```bash
bash run_smk.sh --directory .test --cores 2
```

`exec_mode` in `config/config.yaml`: `native` | `conda` | `container`.

`entrypoint`: `all` | `totals` | `rpfs` | `downstream` | `qc`.

`container` uses Apptainer (`apptainer exec` / `--sdm apptainer`) and each tool's `container_image`. There is no Docker mode.

Full chr20 smoke data: `bash .test/fetch_testdata.sh` (from [SiYangming/Ribo-seq](https://github.com/SiYangming/Ribo-seq) `test/`).

## Authors

- Yangming Si

## References

> Köster, J. et al. _Sustainable data analysis with Snakemake_. F1000Research, 2021.
> Bushell-lab/Ribo-seq pipeline (fork: SiYangming/Ribo-seq).
