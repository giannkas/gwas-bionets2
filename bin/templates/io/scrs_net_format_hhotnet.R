#!/usr/bin/env Rscript

# NETWORK
# SCORES
# index_gene.tsv, edge_list.tsv, gene_score.tsv

library(tidyverse)
library(igraph)

scores <- read_tsv("${SCORES}")
network <- read_tsv("${NETWORK}")

net_genes <- network %>%
  rename(gene1 = `Official Symbol Interactor A`, gene2 = `Official Symbol Interactor B`) %>%
  select(gene1, gene2) %>%
  mutate(gene1 = toupper(gene1), gene2 = toupper(gene2))

# Index-to-gene file
idx_gen <- scores %>%
  select(Gene) %>%
  mutate(Gene = toupper(Gene)) %>%
  distinct(Gene)

matched_genes <- net_genes %>%
  filter(gene1 %in% idx_gen\$Gene & gene2 %in% idx_gen\$Gene) %>%
  select(gene1, gene2)  %>%
  pivot_longer(cols = c(gene1, gene2), names_to = "gene_type", values_to = "genes") %>%
  distinct(genes) %>%
  arrange(genes)

idx_gen <- scores %>%
  select(Gene) %>%
  mutate(Gene = toupper(Gene)) %>%
  distinct(Gene) %>%
  filter(Gene %in% matched_genes\$genes) %>%
  mutate(GeneID = row_number()) %>%
  select(GeneID, Gene)

# Gene-to-score file
gen_scr <- scores %>%
  select(Gene, Pvalue) %>%
  mutate(Gene = toupper(Gene)) %>%
  filter(Gene %in% matched_genes\$genes) %>%
  distinct(Gene, .keep_all = TRUE)

# Edge list file
# Mapping between gene names and IDs
gene_map <- setNames(idx_gen\$GeneID, idx_gen\$Gene)

edg_lst <- network %>%
  rename(gene1 = `Official Symbol Interactor A`, gene2 = `Official Symbol Interactor B`) %>%
  select(gene1, gene2) %>%
  mutate(gene1 = gene_map[toupper(gene1)], gene2 = gene_map[toupper(gene2)]) %>%
  filter(!is.na(gene1) & !is.na(gene2))


write_tsv(idx_gen, "index_gene.tsv", col_names = FALSE)
write_tsv(gen_scr, "gene_score.tsv", col_names = FALSE)
write_tsv(edg_lst, "edge_list.tsv", col_names = FALSE)
