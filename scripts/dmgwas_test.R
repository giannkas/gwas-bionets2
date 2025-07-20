library(dmGWAS)
library(readr)
library(dplyr)

# scores <- read_tsv('../results/psoriasis/magma/magma_output_50kb_window_wes_ncbi38_tab2_dsl2/stable_consensus/data_splits/data_split_1/pso.scores.genes.out_converted_split_1.tsv') %>% 
#   select(Gene, Pvalue) %>%
#   mutate(Pvalue = ifelse(Pvalue == 1, 0.99999, Pvalue)) %>%
#   as.data.frame
scores <- read_tsv('results/AT4G20480_1/bash_scores/bash_output_10kb_window/AT4G20480.1_gene_pvalues.tsv') %>% 
  select(Gene, Pvalue) %>%
  mutate(Pvalue = ifelse(Pvalue == 1, 0.99999, Pvalue)) %>%
  as.data.frame
net <- read_tsv('data/genes/AraNet.txt') %>%
  rename(Interactor_A = `Official Symbol Interactor A`, 
         Interactor_B = `Official Symbol Interactor B`) %>%
  select(Interactor_A, Interactor_B) %>%
  filter(Interactor_A != Interactor_B) %>%
  as.data.frame

modules <- dms(net, scores, expr1 = NULL, expr2 = NULL, r = 0.1, d = 2)
top <- simpleChoose(modules)

tibble(Gene = names(V(top$subnetwork))) %>%
  write_tsv("selected_genes.dmgwas.txt")