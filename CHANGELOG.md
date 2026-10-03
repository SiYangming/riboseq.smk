# Changelog

## 1.0.0 (2026-10-03)


### Bug Fixes

* add conda env to stage_fastq for snakemake lint ([3af69cc](https://github.com/SiYangming/riboseq.smk/commit/3af69cc97c4b6cd98f0526666e44808c0c6fd985))
* disambiguate FastQC rules and vendor tiny CI FASTQs ([22bc2d1](https://github.com/SiYangming/riboseq.smk/commit/22bc2d18c49261950d954f864eb5b97e9f693277))
* move helper functions into common.smk for snakemake lint ([2669afc](https://github.com/SiYangming/riboseq.smk/commit/2669afc8cc2f8c2f38b07d0444d255f87ee3bc6c))
* pass snakefmt, prettier, and snakemake lint on CI ([df134a6](https://github.com/SiYangming/riboseq.smk/commit/df134a639a8bd9a753a06eccd3df82f853c75ad0))

## [0.1.0](https://github.com/SiYangming/riboseq.smk/releases/tag/v0.1.0) (2026-10-03)

### Features

* Initial Snakemake port of SiYangming/Ribo-seq (RPFs, Totals, Downstream) with native, conda, and container exec modes.
