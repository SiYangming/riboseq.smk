# riboseq.smk：将 riboseq.sh_r_py 与 nf-core/riboseq 迁移到 Snakemake 的规划

本目录用于承载基于 **Snakemake** 重构的 Ribo-seq / RNA-seq 工作流程。

当前你有两套主要流程：

- `/Users/siyangming/nextflow_nf_core/riboseq.smk/riboseq.sh_r_py`：
  - 经典 Shell + Python + R Ribo-seq 流程（RPFs + Totals），以 `common_variables.sh` 为中心，手动分步脚本。
  - 主要工具：fastqc、cutadapt、UMI-tools、bbmap、bowtie2、samtools、RSEM，以及一系列 Python 计数脚本与 R 下游分析脚本。
- `/Users/siyangming/nextflow_nf_core/riboseq.nf`：
  - nf-core/riboseq 的 Nextflow 实现，模块化、容器化，支持 Ribo-seq + RNA-seq + TE 分析。
  - 主要使用 nf-core modules（fastp/TrimGalore、fastqc、STAR、Salmon、SortMeRNA/bbsplit、ribowaltz、ribotish、ribotricer、anota2seq、multiqc 等），通过 `assets/samplesheet.csv` + `conf/*.config` 驱动。

本规划文件的目标：在不破坏现有两个流程的前提下，设计一套统一的 Snakemake workflow 迁移方案，并按优先级拆分为多个阶段，便于逐步实现和回归测试。

---

## 一、现有两套流程的执行逻辑对比（摘要）

### 1. riboseq.sh_r_py：Shell + Python + R 流程

- **输入组织方式**：
  - 在 `Shell_scripts/common_variables.sh` 中写死 `RPF_filenames` 和 `Totals_filenames`，以及 `parent_dir`、FASTQ/FASTQ_C、BAM、Analysis、plots 等所有路径。
  - R 端 `R_scripts/common_variables.R` 需要与之保持一致。

- **Totals（RNA-seq）主线**：
  - `Totals_0_QC.sh`：raw FASTQ → fastqc 质量评估。
  - `Totals_1_adaptor_removal.sh`：cutadapt 去接头 + 质控剪切 + fastqc。
  - `Totals_2_extract_UMIs.sh`：UMI-tools extract 提取 UMI + fastqc。
  - `Totals_3a_align_reads_transcriptome.sh`：bowtie2 比对到转录组索引（兼容 RSEM）+ samtools sort/index。
  - `Totals_4a_deduplication_transcriptome.sh`：UMI-tools dedup 去重复。
  - `Totals_5_isoform_quantification.sh`：RSEM gene/isoform 定量。
  - `Totals_6a_write_most_abundant_transcript_fasta.sh` + `calculate_most_abundant_transcript.R`：生成每个基因的最丰度转录本列表和 FASTA。
  - `Totals_6b_extract_read_counts.sh`：调用 `extract_read_counts.py` 从日志提取 reads 流程统计。

- **RPFs（Ribo-seq）主线**：
  - `RPFs_0_QC.sh`：raw FASTQ → fastqc。
  - `RPFs_1_adaptor_removal.sh`：cutadapt 去接头 + fastqc。
  - `RPFs_2_extract_UMIs.sh`：UMI-tools extract（5'/3' 各 4 nt UMI）+ fastqc。
  - `RPFs_3_align_reads.sh`：
    - bbmap → rRNA（拆分 rRNA / non-rRNA），
    - bbmap → tRNA（拆分 tRNA / non-rRNA_tRNA），
    - bbmap → 最丰度转录本 FASTA（pc），输出 BAM + pc/unaligned FASTQ + fastqc。
  - `RPFs_4_deduplication.sh`：UMI-tools dedup + samtools sort/index。
  - `RPFs_5_Extract_counts_all_lengths.sh`：
    - `samtools faidx` most_abundant_fasta；
    - 对每个样本、每个 read length（25–35）调用 `counting_script.py` 生成 `*_pc_L{len}_Off0.counts`。
  - `RPFs_6a_summing_region_counts.sh`：`summing_region_counts.py` + `region_lengths`，生成 5'UTR/CDS/3'UTR 区域计数。
  - `RPFs_6b_summing_spliced_counts.sh`：`summing_spliced_counts.py`，生成 start/stop/UTR 边界 profile。
  - `RPFs_6c_periodicity.sh`：`periodicity.py`，统计 CDS 内 f0/f1/f2 frame 计数。
  - `RPFs_6d_extract_read_counts.sh`：`extract_read_counts.py {sample} RPFs`，统计各步骤 reads 数量。
  - `RPFs_7_Extract_final_counts.sh`：针对最终选定的长度和 offset 调用 `counting_script.py`，生成 `*_pc_final.counts`。
  - `RPFs_8a/8b/8c/8d`：`summing_CDS_counts.py`、`summing_UTR5_counts.py`、`counts_to_csv.py`、`count_codon_occupancy.py`，得到 CDS/UTR5/codon 级别的计数和矩阵。

- **R 下游分析主线**（依赖 Shell + Python 输出）：
  - DESeq2：`DESeq2_Totals.R` / `DESeq2_RPFs.R` / `DESeq2_TE.R` + `DE_analysis_plots.R`。
  - QC：`QC/` 下多个脚本，通过 read_counts / region_counts / periodicity 等进行可视化。
  - meta_plots / codon_occupancy / GSEA / feature_properties 等模块进一步挖掘。

### 2. nf-core/riboseq：Nextflow 流程

Nextflow 顶层入口是 `/Users/siyangming/nextflow_nf_core/riboseq.nf/main.nf`：

- 通过 `PIPELINE_INITIALISATION` 做参数校验、日志设置等；
- `PREPARE_GENOME` 子工作流根据 `--fasta/--gtf/--gff/--additional_fasta/--gencode` 等生成或加载：
  - genome/transcript FASTA，
  - STAR / Salmon / SortMeRNA / bbsplit 等索引，
  - rRNA 数据库、chrom sizes 等；
- 主工作流 `RIBOSEQ`（在 `workflows/riboseq`，此处不展开全部细节）根据 samplesheet 运行：
  - fastq lint（fq lint）验证 FASTQ；
  - fastqc、Trimgalore/fastp adapter trimming；
  - rRNA 移除（SortMeRNA 或 bbsplit）；
  - STAR 基因组比对 + Salmon pseudo alignment；
  - 支持 UMI（`--with_umi`，使用 umi-tools）；
  - ribowaltz P-site 识别与诊断；
  - ribotish QC + ORF prediction；
  - ribotricer ORF 检测；
  - anota2seq TE 分析（结合 Ribo-seq + RNA-seq）；
  - multiqc 整合所有 QC 结果；
- `PIPELINE_COMPLETION` 负责邮件 / slack 报告和结束信息。

输入 samplesheet 结构在 `docs/usage.md` 中约定，核心字段：

- `sample, fastq_1, fastq_2, strandedness, type, ...`，其中 `type` 区分 `riboseq` / `rnaseq` / `tiseq`。

总体来说，nf-core/riboseq 是一个**高度模块化、物种与平台无关、支持 UMI 与 TE 的通用化流程**；而 riboseq.sh_r_py 更偏向特定实验/文库类型的手工脚本组合。

---

## 二、Snakemake 迁移的总体设计思路

在 Snakemake 中希望达到的目标是：

1. 统一样本与参数配置：
   - 使用 `config/config.yaml` + `config/samplesheet.csv`（或 `samples_rpf.tsv` / `samples_totals.tsv`）统一描述样本、分组和参考信息；
   - 尽量兼容 nf-core 的 samplesheet 列（`sample, fastq_1, fastq_2, strandedness, type, treatment, pair`），方便在 Snakemake 和 Nextflow 之间对比和复用。

2. 统一目录结构与输出：
   - 在 `riboseq.smk/` 下建立 Snakemake 标准结构：
     - `workflow/Snakefile`
     - `workflow/rules/`：拆分成 `preprocess.smk`、`align.smk`、`counts.smk`、`riboseq_legacy.smk`、`riboseq_nfcore_like.smk` 等；
     - `workflow/envs/`：不同工具的 conda 环境（或 container 定义）；
     - `workflow/config/`：`config.yaml`, `samples_*.tsv`；
     - `workflow/scripts/`：必要的辅助 Python 脚本（大部分可直接引用已有 `Python_scripts/`）。
   - 输出目录建议保留与现有 Shell 流程兼容的结构（`Counts_files/`、`Analysis/`、`plots/`），便于 R 脚本在过渡期继续使用；同时对 nf-core 风格结果（multiqc 等）在 `results/` 下增加子目录。

3. 逻辑上分为两条主干，公用尽可能多的前处理：
   - **“legacy” 主干**：严格复现 `riboseq.sh_r_py` 的生物学逻辑，包括 UMI 处理策略、bbmap + bowtie2 + RSEM、counting_script + 一系列 Python / R 下游脚本；
   - **“nf-core-like” 主干**：尽量复刻 nf-core/riboseq 的模块组合（fastp/TrimGalore + STAR + Salmon + ribowaltz + ribotish + ribotricer + anota2seq + multiqc），方便对照验证和未来过渡；
   - 共用 reference 准备、fastqc、samplesheet 解析等 rule，以减少重复实现。

4. 参考处理策略：
   - Snakemake 中添加一个 `prepare_reference` 工作流：
     - 若用户提供 `fasta` + `gtf`：则调用本地脚本或移植 nf-core 的逻辑生成 STAR/Salmon/SortMeRNA/bbsplit 索引；
     - 若用户提供现成索引：则直接在 `config.yaml` 中指定，Snakemake 仅做存在性检查；
   - 与 nf-core 类似，支持 `--gencode` / `--skip_gtf_filter` 等配置在 `config.yaml` 中开关；
   - 保留 `filter_gtf.py`、`gtf2bed`、`preprocess_transcripts_fasta_gencode` 等 local modules 的功能，但以 Snakemake rule 形式包装。

5. 充分利用现有 Python / R 代码：
   - Snakemake 的 rule 尽量直接调用现有脚本（`Python_scripts/*.py`、`R_scripts/*.R`），通过命令行参数传入 input/output，而不是大改脚本本身；
   - 统一在 rule 的 `conda:` 字段中指定对应环境，确保依赖可重现；
   - 未来如需将逻辑收敛到 Snakemake 原生 `script:`，可以再做二次重构，但第一阶段以“包装”为主。

---

## 三、与 riboseq.sh_r_py 的 Snakemake 映射

按功能，将 Shell + Python + R 流程拆分为 Snakemake rule 组：

1. **配置与目录创建**

- 使用 `config/config.yaml` 替代 `common_variables.sh`：
  - `parent_dir`、所有子目录路径；
  - adapters（RPF/Totals）、UMI 模式、read length/offset 范围；
  - FASTA / index / region_lengths 路径。
- 通过 `directory()` / `touch` 目标自动创建 `fastq_files/`、`BAM_files/`、`Analysis/`、`plots/` 等目录，无需再运行 `makeDirs.sh`。

2. **Totals（RNA-seq）Snakemake rules**

- `rule totals_fastqc_raw`：raw FASTQ → fastqc。
- `rule totals_cutadapt`：raw FASTQ → cutadapt 修剪 → `*_cutadapt.fastq.gz` + log。
- `rule totals_umi_extract`（可选）：UMI-tools extract → `*_UMI_clipped.fastq.gz` + log + fastqc。
- `rule totals_align_transcriptome`：bowtie2 → `${sample}_pc_sorted.bam` + `.bai` + log。
- `rule totals_dedup_transcriptome`：UMI-tools dedup → `${sample}_pc_deduplicated.bam` + log。
- `rule totals_rsem_quant`：RSEM → `${sample}.genes.results` / `.isoforms.results`。
- `rule totals_extract_read_counts`：`extract_read_counts.py {sample} Totals -log_dir ...` → `${sample}_read_counts.csv`。
- `rule totals_calculate_most_abundant_transcripts`：调用 `calculate_most_abundant_transcript.R` → `most_abundant_transcripts_IDs.csv` / `.txt`。
- `rule totals_write_most_abundant_fasta`：`filter_FASTA.py` → `most_abundant_transcripts.fa`，供 RPF 流程和 legacy counting 使用。

3. **RPFs（Ribo-seq）Snakemake rules**

- `rule rpf_fastqc_raw`：raw FASTQ → fastqc。
- `rule rpf_cutadapt`：raw FASTQ → cutadapt 修剪 → `*_cutadapt.fastq.gz` + log。
- `rule rpf_umi_extract`：UMI-tools extract（regex pattern）→ `*_UMI_clipped.fastq.gz` + log + fastqc。
- `rule rpf_align_rrna` / `rule rpf_align_trna`：bbmap 比对到 rRNA/tRNA → 对齐和未对齐 FASTQ + logs。
- `rule rpf_align_pc`：
  - input：`*_non_rRNA_tRNA.fastq.gz` + `most_abundant_transcripts.fa`；
  - output：`${sample}_pc_sorted.bam` + `.bai` + `pc/unaligned.fastq.gz`。
- `rule rpf_dedup`：UMI-tools dedup + samtools sort/index → `${sample}_pc_deduplicated_sorted.bam` + `.bai`。
- `rule rpf_counts_all_lengths`：
  - 先确保 `most_abundant_transcripts.fa.fai` 存在；
  - 对每个 `{length}` 调用 `counting_script.py`，输出 `*_pc_L{len}_Off0.counts`；
- `rule rpf_region_counts`：`summing_region_counts.py` → `_region_counts.csv`。
- `rule rpf_spliced_counts`：`summing_spliced_counts.py` → `_start_site.csv` 等；
- `rule rpf_periodicity`：`periodicity.py` → `_periodicity.csv`。
- `rule rpf_extract_read_counts`：`extract_read_counts.py {sample} RPFs` → `${sample}_read_counts.csv`。
- `rule rpf_final_counts`：使用 config 中最终 `lengths`/`offsets`，调用 `counting_script.py` → `*_pc_final.counts`。
- `rule rpf_cds_counts` / `rule rpf_utr5_counts` / `rule rpf_counts_to_csv` / `rule rpf_codon_occupancy`：分别包装 `summing_CDS_counts.py`、`summing_UTR5_counts.py`、`counts_to_csv.py`、`count_codon_occupancy.py`。

4. **R 分析模块的 Snakemake 封装（后续阶段）**

- 以简单包装为主：
  - `rule deseq2_totals` / `rule deseq2_rpfs` / `rule deseq2_te` / `rule deseq2_plots`；
  - `rule meta_normalise_counts` / `rule meta_bin_data` / `rule meta_plot_binned` / `rule meta_plot_individual`；
  - `rule fgsea` / `rule rrvgo` / `rule codon_occupancy_plot` 等；
- 每个 rule 的 `input:` 明确指向前面 rules 产生的 counts/annotation 文件，不再依赖隐式路径。

---

## 四、与 nf-core/riboseq 的 Snakemake 映射

为保持与 nf-core/riboseq 的思路一致，在 Snakemake 中单独设计一组“nf-core-like” rules：

1. **样本与对比信息**

- 直接复用 nf-core 的 samplesheet 格式：`assets/samplesheet.csv`；
- 在 `config/config.yaml` 中增加 `samplesheet: "assets/samplesheet.csv"` 和 `contrasts: path/to/contrasts.csv` 字段；
- 在 Snakemake 的 `config` 解析中读取 `type`（riboseq/rnaseq）、`treatment`、`pair` 信息，为 TE 分析和 pairing 提供输入。

2. **参考基因组准备（对应 PREPARE_GENOME 子工作流）**

- `rule prepare_reference_fasta_gtf`：接受 `fasta` + `gtf`（或 `gff`），输出清洗后的 GTF、FASTA、fai、chrom_sizes 等；
- `rule build_star_index`：STAR genomeGenerate；
- `rule build_salmon_index`：salmon index；
- `rule build_sortmerna_index`：SortMeRNA index；
- `rule build_bbsplit_index`：bbsplit index；
- `rule preprocess_transcript_fasta_gencode`：迁移模块 `modules/local/preprocess_transcripts_fasta_gencode` 的逻辑，生成转录本 FASTA + tx2gene 等；
- 所有索引在 Snakemake 中作为显式输出，由后续 alignment / pseudo alignment / rRNA 去除 rule 依赖。

3. **FASTQ 预处理和 QC**

- `rule fqlint_raw` / `rule fqlint_after_trimming`：包装 nf-core `fq/lint` 模块逻辑（调用 `fq` 二进制），做输入数据 lint；
- `rule fastqc_raw` / `rule fastqc_trimmed`：复用 nf-core fastqc 模块的参数组合；
- `rule trimmer`：根据 `config["trimmer"]` 选择 TrimGalore! 或 fastp：
  - TrimGalore 路径直接参考 `docs/usage.md` 中推荐参数；
  - fastp 保留常用参数并暴露 `extra_fastp_args` 给 config；
- 这些步骤与 legacy 主干的 cutadapt/fastqc 可以部分共享（通过 wrapper 或 rule 复用）。

4. **rRNA/qc 对应的 rules**

- `rule rrna_depletion_sortmerna`：SortMeRNA 将 reads 分为 rRNA / non-rRNA；
- `rule rrna_depletion_bbsplit`：bbsplit alternative，实现与 nf-core `bbmap/bbsplit` module 类似的行为；
- 输出与 MultiQC 兼容的 log 格式，方便后续汇总。

5. **比对与定量**

- `rule star_align_genome`：STAR 比对到 genome，并将 alignment 投射到 transcriptome（与 nf-core 相同）；
- `rule salmon_quant`：Salmon quant 对 transcript-level quantification；
- `rule kallisto_index` / `rule kallisto_quant`（可选）：迁移 Kallisto modules；
- UMI 支持：
  - `rule umi_extract`：根据 `with_umi` / `umitools_extract_method` / `umitools_bc_pattern` 等参数调用 umi-tools；
  - `rule umi_dedup`：对应 umi-tools dedup，支持 grouping_method 配置。

6. **Ribo 特异模块**

- `rule ribowaltz`：包装 `modules/nf-core/ribowaltz` 模块逻辑：
  - 输入 STAR BAM + annotation（GTF/tx2gene）；
  - 输出 P-site 识别与周期性 QC 结果；
- `rule ribotish_quality` / `rule ribotish_predict`：调用 ribotish 对 ORF 做 QC 和预测；
- `rule ribotricer_prepareorfs` / `rule ribotricer_detectorfs`：准备 ORF 和 Ribo-seq reads，运行 ribotricer 检测 ORF；
- 这些 rule 可以作为 nf-core-like 主干的一部分，与 legacy counting_script 路径并行存在。

7. **TE 分析与报告（anota2seq + multiqc）**

- `rule anota2seq_run`：
  - 输入：`type == rnaseq` 和 `type == riboseq` 的 quant/counts + contrasts 表；
  - 输出：TE 分析结果（表格 + plots），行为上对齐 nf-core/anota2seq module；
- `rule multiqc`：整合 fastqc、trimmer、alignment、rRNA depletion、ribowaltz、ribotish、ribotricer 等所有 logs，生成 MultiQC 报告。

---

## 五、按优先级划分的分阶段实施计划

下面给出一个面向实际开发的分阶段计划，优先级从高到低（P0 > P1 > P2 > P3 > P4）。

### 阶段 P0（最高优先级）：Snakemake 框架与配置打底

**目标**：在不实现任何具体分析规则的前提下，搭建好 Snakemake 项目的骨架与配置系统。

- 在 `riboseq.smk/` 下创建：
  - `workflow/Snakefile`（包含一个空的 `rule all` 和 config 载入逻辑）；
  - `config/config.yaml`（只写入基础结构即可）；
  - `testdata/samplesheet_local.csv`（作为默认示例样本表，兼容 nf-core 字段）；
  - `workflow/envs/`（先留空，用于后续添加 conda 环境）。
- 将 `common_variables.sh` 和 `R_scripts/common_variables.R` 中的关键信息抽象进 `config.yaml`：
  - `paths`、`adapters`、`umi`、`reference`、`rpf_lengths`、`rpf_offsets` 等 sections；
  - 记录 nf-core 风格参数（`with_umi`、`skip_gtf_filter`、`remove_ribo_rna`、`aligner` 等）。
- 目标产物：
  - `snakemake -n` 能成功读取 config，并打印空 DAG，无报错。

### 阶段 P1（高优先级）：先迁移 riboseq.sh_r_py 的 Totals + RPF 主干

**目标**：在 Snakemake 中重现当前 Shell + Python 流程的核心功能，使现有 R 脚本可以直接复用新输出。

- Totals 流程：
  - 实现 `totals_fastqc_raw` / `totals_cutadapt` / `totals_umi_extract` / `totals_align_transcriptome` / `totals_dedup_transcriptome` / `totals_rsem_quant` / `totals_extract_read_counts` rules；
  - 输出路径保持与原脚本兼容（`BAM_files/`、`rsem/`、`logs/` 等）；
  - 使用 `envs/` 下的 conda 环境配置 fastqc/cutadapt/umi-tools/bowtie2/samtools/RSEM。

- RPFs 流程：
  - 实现 `rpf_fastqc_raw` / `rpf_cutadapt` / `rpf_umi_extract` / `rpf_align_rrna` / `rpf_align_trna` / `rpf_align_pc` / `rpf_dedup` rules；
  - 实现 `rpf_counts_all_lengths` / `rpf_region_counts` / `rpf_spliced_counts` / `rpf_periodicity` / `rpf_extract_read_counts` / `rpf_final_counts` / `rpf_cds_counts` / `rpf_utr5_counts` / `rpf_counts_to_csv` / `rpf_codon_occupancy` 包装 rule；
  - 所有计数脚本直接调用现有 Python 文件，先不改动脚本逻辑。

- 验证：
  - 使用 `testdata/GSE182201` 作为输入，比较 Snakemake 输出与运行 Shell 脚本（Totals_all.sh / RPFs_all.sh）时生成的关键文件是否一致（或在可接受差异范围内）。

### 阶段 P2（高优先级）：迁移 nf-core/riboseq 的参考准备与预处理逻辑

**目标**：在 Snakemake 中提供一个与 PREPARE_GENOME + 基础预处理等价的工作流，使你可以在 Snakemake 中运行“nf-core 风格”的 Ribo-seq/RNA-seq 流程。

- 参考准备：
  - 实现 `prepare_reference_fasta_gtf` / `build_star_index` / `build_salmon_index` / `build_sortmerna_index` / `build_bbsplit_index` / `preprocess_transcript_fasta_gencode` rules；
  - 在 `config.yaml` 中允许两种模式：自动生成索引（默认）与使用已有索引（指定路径，跳过构建）；
  - 使 `snakemake --cores N reference_all` 可以单独跑参考准备。

- FASTQ 预处理与 QC：
  - 实现 `fqlint_raw` / `fastqc_raw` / `trimmer` / `fqlint_after_trimming` / `fastqc_trimmed`，复用 nf-core 推荐参数；
  - 抽象 nf-core 中的 `--trimmer` / `--extra_trimgalore_args` / `--extra_fastp_args` 等为 Snakemake config 字段。

### 阶段 P3（中优先级）：整合 nf-core 的 Ribo 特异模块与 TE 分析

**目标**：在 Snakemake 中提供与 nf-core/riboseq 功能上尽量对齐的 Ribo-specific 和 TE 分析能力。

- 实现 `ribowaltz` / `ribotish_quality` / `ribotish_predict` / `ribotricer_prepareorfs` / `ribotricer_detectorfs` rules：
  - 依赖 STAR BAM 和 annotation，输出 QC 与 ORF 结果；
  - 输出目录与 nf-core 命名保持尽可能一致，便于复用其 R 分析脚本或下游工具。

- 实现 `anota2seq_run`：
  - 输入 RNA-seq + Ribo-seq 的 counts/quant + contrasts 表；
  - 输出 TE 分类结果与可视化；
  - 可优先支持最常用的单因素/两组比较场景，其余复杂设计后续迭代。

- 将多种 QC / TE 结果串到一个 `rule all_nfcore_like` 中，方便单命令跑完 nf-core 风格主干。

### 阶段 P4（中低优先级）：统一 QC/报告与 R 分析 Snakemake 化

**目标**：减少手工运行 R 的步骤，使 Snakemake 从原始 FASTQ 一路驱动到所有主要 R 报告和图形。

- 将现有 R 脚本系统性包装为 Snakemake rules：
  - QC 模块（`RPFs_read_counts.R`、`Totals_read_counts.R`、`heatmaps.R`、`periodicity.R` 等）；
  - meta_plots（normalise/bin/plot）与单基因 profile；
  - GSEA / rrvgo / Interactive_scatters / codon_occupancy 等；
- 为每一类分析定义对应的 `rule all_*`，例如 `all_qc`, `all_meta_plots`, `all_gsea`，方便按模块触发。
- 引入 `multiqc` rule 整合 fastqc/trimmer/alignment/rRNA/ribowaltz 等 logs，与 nf-core 风格保持一致。

### 阶段 P5（低优先级）：性能优化、测试与文档完善

**目标**：在功能稳定后对 Snakemake 流程做系统优化和工程化改造。

- 性能与资源管理：
  - 为关键 rule 增加 `resources` 与 `threads`（统一从 `config.yaml` 读），优化集群/云环境调度；
  - 利用 `group` / `localrules` 对极短小步骤进行合并，减少调度开销；
  - 使用 `temp()`/`protected()` 控制中间文件生命周期，降低存储压力。

- 测试与 CI：
  - 利用 `testdata/GSE182201` 构建最小测试集，在 Snakemake 中添加 `profile` 或 `--config test=true` 的测试路径；
  - 如有需要，可使用 `pytest + snakemake` 或 `nf-test` 风格工具做规则级测试；
  - 对齐 nf-core 现有 GitHub Actions 流水线，在本仓库增加简单的 CI（lint + dry-run + 小数据跑通）。

- 文档：
  - 在本 README 基础上，补充 Snakemake 安装说明、典型运行命令、配置示例等；
  - 结合 Snakemake `--dag` / `--rulegraph` 输出，将流程图截图放入文档便于用户理解。

---

## 六、项目概述与安装配置指南

### 1. 项目概述

本项目在同一代码库中集成三套分析路径：

- 传统 Shell + Python + R 的 riboseq.sh_r_py 流程；
- nf-core/riboseq 的 Nextflow 实现；
- 基于 Snakemake 的新工作流 riboseq.smk，用于逐步接管上述两套流程的主要功能。

Snakemake 工作流的目标是：

- 用 config/config.yaml 和样本表统一配置；
- 通过模块化 rules 复用现有 Python/R 代码；
- 支持从小规模测试数据到生产级数据集的扩展；
- 在必要时复刻 nf-core/riboseq 的分析路径，便于结果对比与验证。

### 2. 环境准备与安装

#### 2.1 基础依赖

- 推荐使用 Conda 或 Mamba 管理环境；
- 需要安装 Snakemake（建议使用 `snakemake` 独立包或 `snakemake-minimal` 搭配 `mamba`）；
- 后续各 rule 将通过 workflow/envs 下的 YAML 文件声明工具依赖（如 fastqc、cutadapt、umi-tools 等）。

#### 2.2 创建 Snakemake 环境（示意）

可以按如下思路在本地或集群上创建专用环境：

- 使用 mamba/conda 创建基础环境并安装 Snakemake；
- 在运行时通过 `--use-conda` 让 Snakemake 自动为每个 rule 创建子环境；
- 如在 HPC 环境，可结合 `--profile` 使用集群配置文件。

### 3. 配置文件与目录结构

当前 Snakemake 骨架遵循如下结构：

- config/config.yaml：全局配置文件，包含路径、样本表位置、参考目录等关键信息；
- workflow/Snakefile：主入口文件，加载配置与核心规则；
- workflow/rules/core.smk：当前集中定义 `rule all` 等骨架规则，后续会按模块拆分；
- workflow/envs/：预留，用于后续添加各工具的 conda 环境定义；
- workflow/scripts/：预留，用于放置 Snakemake 封装用的 Python/R 脚本；
- reference/：用于存放小型通用参考与注释文件，并通过 README 记录大型公开数据的获取方式；
- testdata/：包含 GSE182201 子集 FASTQ、示例样本表与对比文件，用于小规模测试。

### 4. Riboseq 分析流程说明（Snakemake 视角）

从 Snakemake 的角度看，计划中的 Ribo-seq/RNA-seq 分析将分为三层：

- 参考准备层：下载/检查参考基因组与注释，构建各类索引；
- 读长处理与比对层：对原始 FASTQ 做 QC、接头去除、UMI 处理、rRNA/tRNA 去除、基因组/转录组比对或定量；
- 下游分析层：计数矩阵生成、差异分析（DESeq2 与 TE 分析）、周期性与 frame 分析、GSEA 与 meta plots 等。

当前仓库中的 Snakemake 骨架已经为上述三层预留了配置与模块划分，后续会在 P1/P2 阶段逐步将 riboseq.sh_r_py 与 nf-core/riboseq 的逻辑映射到具体 rules 中。

### 5. 使用示例

在完成基础环境配置后，可以在 riboseq.smk 目录下运行：

- `snakemake -n`：在 dry-run 模式下检查当前骨架是否能正常解析配置并构建 DAG；
- `snakemake -n --configfile config/config.yaml`：显式指定配置文件进行 dry-run；
- 后续在实现具体 rules 后，可增加目标规则，例如 `snakemake -n all_totals` 或 `snakemake -n all_rpf` 等。

### 6. 常见问题（FAQ）

**Q1：当前 Snakemake 骨架是否已经执行实际计算？**

A1：目前 `rule all` 仅作为占位符，不触发真实计算，主要用于验证配置与目录结构是否正确。后续阶段会在不改变整体结构的前提下逐步填充具体 rules。

**Q2：样本表与对比文件在哪里配置？**

A2：示例样本表与对比文件位于 testdata 目录中，并在 config/config.yaml 的 samples 字段中引用：

- `samples.sheet: "testdata/samplesheet_local.csv"`
- `samples.contrasts: "testdata/contrasts_local.csv"`

未来可根据需要将其替换为项目实际的样本表与对比设计文件，但保持相同的列结构有助于重用现有 R 与 Python 分析脚本。

**Q3：reference 目录下应该放哪些文件？**

A3：reference/ 仅用于存放体积较小、被多个项目或多个分析步骤复用的参考文件，例如 rRNA/tRNA 序列子集、region length 表等。大型 FASTA/GTF 及索引不直接纳入仓库，而是在 README 中记录获取方式，并在 config/config.yaml 中配置本地路径。

**Q4：如何与现有 riboseq.sh_r_py 与 nf-core/riboseq 保持一致性？**

A4：Snakemake 工作流会在规则命名、输入输出路径与参数设置上尽量贴近现有脚本和 nf-core 模块。通过测试数据对比关键输出（计数矩阵、reads 统计等），可以在每一阶段验证行为一致性。

### 7. 变更记录（与本次任务相关）

- 新增 config/config.yaml，集中管理路径、样本表与参考目录等基础配置；
- 新增 workflow/Snakefile 与 workflow/rules/core.smk，构建 Snakemake 工作流骨架并定义占位 `rule all`；
- 新增 reference/README.md，规范参考文件管理策略，并说明不在仓库中存储大型公开数据；
- 更新 testdata/samplesheet_local.csv，使 FASTQ 路径指向当前仓库内的 testdata/GSE182201 目录；
- 更新 testdata/contrasts_local.csv，增加 design_path 列以标记对比设计文件在仓库中的默认位置；

---

## 七、后续建议

- 若优先目标是尽快替换掉手工 Shell 流程，建议严格按照 P0 → P1 → P2 的顺序推进，先让 `riboseq.sh_r_py` 的功能在 Snakemake 中稳定落地，再逐步引入 nf-core/riboseq 的模块；
- 在 P1 完成后，就可以开始用 Snakemake 驱动现有 R 分析脚本，从而在同一 Snakemake DAG 中覆盖“预处理 + 比对 + 计数 + 差异分析 + meta plots”；
- 在 P2/P3 完成后，可以考虑逐步弃用 Nextflow 版本，只保留 Snakemake 作为统一入口，同时仍保留 nf-core 风格的参数与模块选择，方便将来与社区保持同步。

以上规划与实现为当前版本，后续可以根据测试结果和使用体验进一步细化和调整。
