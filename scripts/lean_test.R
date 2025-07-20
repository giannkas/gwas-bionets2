library(igraph)
library(LEANR)
library(dplyr)
library(tidyverse)
#set.seed(123456)
network <- read_tsv("../data/biogrid_ppi_interactorsAB.tsv")
gene_scores <- read_tsv("../results/psoriasis/magma/magma_output_50kb_window_wes_ncbi38_tab2_dsl2/stable_consensus/data_splits/data_split_1/pso.scores.genes.out_converted_split_1.tsv")

net_genes <- network %>%
  rename(gene1 = `Official Symbol Interactor A`, gene2 = `Official Symbol Interactor B`) %>%
  select(gene1, gene2) %>%
  mutate(gene1 = toupper(gene1), gene2 = toupper(gene2))

g_scrs <- gene_scores %>%
  select(Gene) %>%
  mutate(Gene = toupper(Gene)) %>%
  distinct(Gene)

matched_genes <- net_genes %>%
  filter(gene1 %in% g_scrs$Gene & gene2 %in% g_scrs$Gene) %>%
  select(gene1, gene2)  %>%
  pivot_longer(cols = c(gene1, gene2), names_to = "gene_type", values_to = "genes") %>%
  distinct(genes) %>%
  arrange(genes)

gen_scr <- gene_scores %>%
  select(Gene, Pvalue) %>%
  mutate(Gene = toupper(Gene)) %>%
  filter(Gene %in% matched_genes$genes) %>%
  distinct(Gene, .keep_all = TRUE)

net <- net_genes %>%
  filter(gene1 %in% g_scrs$Gene & gene2 %in% g_scrs$Gene) %>%
  graph_from_data_frame(directed = FALSE)

scores <- gen_scr[['Pvalue']]
names(scores) <- gen_scr[['Gene']]

results <- run.lean(scores, net, n_reps = 100000, add.scored.genes = TRUE, verbose = TRUE, ncores = 8)

significant_genes <- tibble(Gene = rownames(results$restab[results$restab[,'PLEAN']<=0.05,]))

if (nrow(significant_genes) == 0) {
  warning("No significant genes found. Creating empty output file.")
  tibble(gene = character(0)) %>%
    write_tsv('selected_genes.lean.txt')
} else {
  nh_genes <- significant_genes %>%
    mutate(neighbors = map(Gene, ~ {
      if (.x %in% names(results$nhs)) {
        results$nhs[[.x]]
      } else {
        character(0)
      }
    })) %>%
    unnest(neighbors)


  nh_genes_unique <- nh_genes %>% 
    rename(gene = neighbors) %>%
    distinct(gene) %>%
    write_tsv('selected_genes.lean.txt')

}

# Example from LEANR page

set.seed(123456)
# load network and CCM p-values
data(g2)
data(CCM.pvals)
data(gene.annots)
LEAN_results<-LEANR:::LEAN_results

LEAN_results<-lapply(names(CCM.pvals),function(ccm){
  run.lean(CCM.pvals[[ccm]], g2, n_reps = 10000, ncores = 8, verbose = TRUE)
})
names(LEAN_results)<-names(CCM.pvals)

# Extract significant local subnetworks
sign.genes<-lapply(LEAN_results,function(LEANres){
  rownames(LEANres$restab[LEANres$restab[,'PLEAN']<=0.05,])
})
names(sign.genes)<-names(CCM.pvals)

get.ls.info()


# modules <- dms(net, scores, expr1 = NULL, expr2 = NULL, r = 0.1, d = 2)
# top <- simpleChoose(modules)
# 
# tibble(Gene = names(V(top$subnetwork))) %>%
#   write_tsv("selected_genes.dmgwas.txt")