library(rtracklayer)
library(tidyverse)
readGFF("nextflow_nf_core/riboseq.smk/reference/Homo_sapiens.GRCh38.111_chr20.gtf") %>%
  select(gene_id) %>%
  distinct() %>%
  write_delim(file = "nextflow_nf_core/riboseq.smk/reference/GRCh38_chr20_gene.csv", 
              col_names = FALSE)
