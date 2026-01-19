# reference 目录说明

本目录用于存放与 riboseq.smk 工作流通用的、小体积参考与注释文件，例如：

- rRNA/tRNA 序列子集或索引元数据
- 区域长度表（region lengths）、转录本注释表等文本文件
- 与测试数据强相关的最小示例注释文件

## 不在仓库中存储的大型公开数据

以下类型的数据不直接存储在仓库中，只在此说明获取方式与本地路径配置示例：

- 全基因组 FASTA（例如人/小鼠参考基因组）
- 完整 GTF/GFF 注释文件
- 大体积的索引文件（STAR、Salmon、SortMeRNA、bbsplit 等）

请根据分析物种与版本，从对应的基因组数据库下载所需文件，并在 config/config.yaml 中通过 reference 相关字段配置本地路径。

## 推荐的本地目录结构

建议在本机或集群上采用类似结构组织参考数据：

- reference/
  - genome/
  - annotation/
  - rrna_trna/
  - indices/

其中 indices 下可按工具进一步分层（star、salmon、sortmerna、bbsplit 等），并在 Snakemake 配置中使用绝对路径或相对路径引用。

## 与 testdata 的关系

当前仓库中的示例注释文件（例如 testdata 目录下的 rRNA/tRNA 相关 FASTA 与列表）主要用于小规模测试。后续如将正式分析迁移到 riboseq.smk Snakemake 工作流中，可考虑将这些经常复用的小型注释文件移动或复制到 reference/ 目录，并在 config/config.yaml 中更新路径。

## 版本控制建议

- 对于小体积、经常复用的注释文件（例如 region length 表、简化注释表），建议直接纳入版本控制，并在文件名中体现物种与版本号。
- 对于大型 FASTA/GTF 及索引，仅在此 README 中记录获取方式、下载命令示例和期望的目录结构，在 config/config.yaml 中使用本地路径，并通过 Snakemake 规则检查文件是否存在。

