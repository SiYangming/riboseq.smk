## Workflow overview

Ribo-seq (RPFs) plus matched total RNA-seq (Totals), then Downstream R analysis. Logic matches [SiYangming/Ribo-seq](https://github.com/SiYangming/Ribo-seq) (`run.sh` pipelines RPFs / Totals / Downstream).

## Input data

Sample sheet (`config/samples.tsv`):

| sample | fastq_1 | fastq_2 | type | treatment |
| ------ | ------- | ------- | ---- | --------- |
| rpf1 | path/to/rpf.fastq.gz | | riboseq | control |
| tot1 | path/to/R1.fastq.gz | path/to/R2.fastq.gz | rnaseq | treated |

`type`: `riboseq` (RPF, typically SE) or `rnaseq` (Totals, SE or PE). Empty `fastq_2` means single-end.

Set FASTA / GTF / region-length paths in `config/config.yaml`. `entrypoint`: `all` | `totals` | `rpfs` | `downstream` | `qc`.

`exec_mode`: `native` | `conda` | `container` (Apptainer; set each tool `container_image`).

Minimal FastQC smoke: `bash run_smk.sh --directory .test --cores 2` with `.test/config` (`entrypoint: qc`). Full chr20 data: `bash .test/fetch_testdata.sh`.
