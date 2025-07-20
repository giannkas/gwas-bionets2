library(dplyr)
magma <- readr::read_table("../results/psoriasis/magma/magma_output_50kb_window_snp-array_significant_snps/pso.genes.annot_only_gene_ids.tsv")
gprofiler2::gconvert(magma$GENE, numeric_ns="ENTREZGENE_ACC") %>% 
  mutate(GENE = as.numeric(input)) %>% 
  select(GENE, name) %>% 
  left_join(magma, by = "GENE") %>%
  rename(Gene = name) %>%
  readr::write_tsv("../results/psoriasis/magma/magma_output_50kb_window_snp-array_significant_snps/pso.genes.annot_gene_ids_names.tsv")
