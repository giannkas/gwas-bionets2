library(igraph)
library(tidygraph)
library(magrittr)
library(ggrepel)
library(tidyverse)
library(cowplot)
suppressMessages(library(network))
suppressMessages(library(sna))
library(ggnetwork)
library(scatterpie)
library(ggupset)
library(RCy3)

theme_transparent <- theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
                           plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                           panel.grid.major = element_blank(), # get rid of major grid
                           panel.grid.minor = element_blank(), # get rid of minor grid
                           legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
                           legend.box.background = element_rect(fill = "transparent", color = NA) # get rid of legend panel bg
)

theme_set(theme_cowplot() + theme_transparent)
results <- '../results/psoriasis/magma/magma_output_50kb_window/stable_consensus/data_splits/'
data <- '../data/additional_knowledge/'

julio_psoriasis_genes <- read_tsv(paste0(data,'julio_psoriasis_genes.tsv'),
                                  col_types = 'cc') %>%
  select(gene)

ran_psoriasis_genes <- read_tsv(paste0(data,'ran_psoriasis_genes.tsv'),
                                col_types = 'cc') %>%
  select(gene)

tsoi_psoriasis_genes <- read_csv(paste0(data,'tsoi_psoriasis_genes.csv')) %>%
  select(`Closest Gene`)

il17_signaling_pathway <- read_tsv(paste0(data,'il17_signaling_pathway.tsv'),
                                   col_types = 'c') %>%
  select(gene)

il23_signaling_pathway <- read_tsv(paste0(data,'il23_signaling_pathway.tsv'),
                                   col_types = 'c') %>%
  select(gene)

tnf_signaling_pathway <- read_tsv(paste0(data,'tnf_signaling_pathway.tsv'),
                                  col_types = 'c') %>%
  select(gene)

ppi <- read_tsv('../data/biogrid_ppi_interactorsAB_match-magma-genes.tsv', col_types = 'cc') %>%
  select(`Official Symbol Interactor A`, `Official Symbol Interactor B`) %>%
  graph_from_data_frame(directed = FALSE) %>%
  as_tbl_graph

ppi_df <- as_tibble(ppi)

methods <- c('heinz','HotNet2', 'SigMod')
pathways <- c('IL17','IL23', 'TNF')
gen_list <- c('Julio', 'Ran2019', 'Tsoi2017')

method_palette <- c('heinz' = '#73471f',
                    'HotNet2' = '#377eb8', 'SigMod' = '#f4888a', 'All' = 'gray50', 
                    'Consensus' = 'black', 'MAGMA' = 'gray50')
pathways_palette <- c('IL17' = '#489ea3',
                      'IL23' = '#ff6f3c', 'TNF' = '#377eb8')
gen_list_palette <- c('Julio' = '#489ea3',
                      'Ran2019' = '#ff6f3c', 'Tsoi2017' = '#377eb8')

labs <- c('all_snps' = 'All', 'consensus' = 'Consensus',
          'heinz' = 'heinz', 'hotnet2' = 'HotNet2', 'sigmod' = 'SigMod')

# TOPOLOGY PANEL

plot_network <- function(net, selected_nodes) {
  
  graph <- activate(net, nodes) %>%
    filter(name %in% selected_nodes)
  class(graph) <- c('igraph','tbl_graph')
  
  ggnetwork(graph) %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = '#3282b8') +
    geom_nodes(color = '#233142') +
    ggtitle(' ') +
    theme_blank() + theme_transparent + 
    theme(title = element_text(hjust = 1, face = 'bold'))
}

# HEINZ genes

gene_lists_heinz <- list()

# Loop through the files and read them into the list
for (i in 1:5) {
  file_name <- paste0(results,"data_subset_", i, "/heinz/selected_genes.heinz.txt")
  gene_lists_heinz[[i]] <- read_tsv(file_name, col_types = 'c')$gene
}

intersection_heinz <- Reduce(intersect, gene_lists_heinz)
heinz_plt <- plot_network(ppi, intersection_heinz)

# HOTNET2 genes

gene_lists_hotnet2 <- list()

# Loop through the files and read them into the list
for (i in 1:5) {
  file_name <- paste0(results,"data_subset_", i, "/hotnet2/selected_genes.hotnet2.tsv")
  gene_lists_hotnet2[[i]] <- read_tsv(file_name, col_types = 'cc')$gene
}

intersection_hotnet2 <- Reduce(intersect, gene_lists_hotnet2)
hotnet2_plt <- plot_network(ppi, intersection_hotnet2)

# SIGMOD genes

gene_lists_sigmod <- list()

# Loop through the files and read them into the list
for (i in 1:5) {
  file_name <- paste0(results,"data_subset_", i, "/sigmod/selected_genes.sigmod.txt")
  gene_lists_sigmod[[i]] <- read_tsv(file_name, col_types = 'c')$gene
}

intersection_sigmod <- Reduce(intersect, gene_lists_sigmod)
sigmod_plt <- plot_network(ppi, intersection_sigmod)

# Topology panel

topology_panel <- plot_grid(heinz_plt, hotnet2_plt, sigmod_plt, 
                            labels = c('Heinz','HotNet2','SigMod'))

f1 <- plot_grid(topology_panel,
                nrow = 1,
                ncol = 2, 
                labels = 'C')
f1


### Consensus of intersections ###

# Heinz

g_heinz <- activate(ppi, nodes) %>%
  filter(name %in% intersection_heinz) %>%
  mutate(heinz = TRUE)
class(g_heinz) <- c('igraph','tbl_graph')

# HotNet2

g_hotnet2 <- activate(ppi, nodes) %>%
  filter(name %in% intersection_hotnet2) %>%
  mutate(hotnet2 = TRUE)
class(g_hotnet2) <- c('igraph','tbl_graph')

# SigMod

g_sigmod <- activate(ppi, nodes) %>%
  filter(name %in% intersection_sigmod) %>%
  mutate(sigmod = TRUE)
class(g_sigmod) <- c('igraph','tbl_graph')

# Consensus

options(repr.plot.width=25, repr.plot.height=15)
consensus <- graph_join(g_heinz, g_hotnet2, by = c('name')) %>%
  graph_join(g_sigmod, by = c('name')) %>%
  to_undirected %>%
  mutate(heinz = ifelse(is.na(heinz), FALSE, heinz),
         hotnet2 = ifelse(is.na(hotnet2), FALSE, hotnet2),
         sigmod = ifelse(is.na(sigmod), FALSE, sigmod),
         num_methods = rowSums(cbind(heinz, hotnet2, sigmod))) %>%
  filter(num_methods >= 2)

as_tibble(consensus) %>%
  rename(gene = name) %>%
  write_tsv(paste0(results,'../stable_consensus.tsv'))

consensus_gg <- consensus %>%
  ggnetwork %>%
  mutate(name = as.character(name))

nodes <- mutate(consensus_gg, heinz = as.numeric(heinz),
                HotNet2 = as.numeric(hotnet2),
                SigMod = as.numeric(sigmod),
                Julio = as.numeric(name %in% julio_psoriasis_genes$gene),
                Ran2019 = as.numeric(name %in% ran_psoriasis_genes$gene),
                Tsoi2017 =  as.numeric(name %in% tsoi_psoriasis_genes$`Closest Gene`),
                IL17 = as.numeric(name %in% il17_signaling_pathway$gene),
                IL23 = as.numeric(name %in% il23_signaling_pathway$gene),
                TNF =  as.numeric(name %in% tnf_signaling_pathway$gene)) %>%
  filter(xend == x & yend == y) %>% 
  select(x, y, name, methods, gen_list, pathways) %>%
  unique

edges <- filter(consensus_gg, xend != x | yend != y)


max_length_pathways <- max(length(nodes$name), length(il17_signaling_pathway$gene), 
                           length(il23_signaling_pathway$gene), length(tnf_signaling_pathway$gene))

max_length_genes <- max(length(nodes$name), length(julio_psoriasis_genes$gene), 
                        length(ran_psoriasis_genes$gene), length(tsoi_psoriasis_genes$`Closest Gene`))

nodes_padded_pathways <- c(nodes$name, rep(NA, max_length_pathways - length(nodes$name)))
il17_padded <- c(il17_signaling_pathway$gene, rep(NA, max_length_pathways - length(il17_signaling_pathway$gene)))
il23_padded <- c(il23_signaling_pathway$gene, rep(NA, max_length_pathways - length(il23_signaling_pathway$gene)))
tnf_padded <- c(tnf_signaling_pathway$gene, rep(NA, max_length_pathways - length(tnf_signaling_pathway$gene)))

nodes_padded_genes <- c(nodes$name, rep(NA, max_length_genes - length(nodes$name)))
julio_padded <- c(julio_psoriasis_genes$gene, rep(NA, max_length_genes - length(julio_psoriasis_genes$gene)))
ran_padded <- c(ran_psoriasis_genes$gene, rep(NA, max_length_genes - length(ran_psoriasis_genes$gene)))
tsoi_padded <- c(tsoi_psoriasis_genes$`Closest Gene`, rep(NA, max_length_genes - length(tsoi_psoriasis_genes$`Closest Gene`)))

upset_consensus_pathways <- data.frame(nodes_padded_pathways, 
                                       il17_padded, 
                                       il23_padded, 
                                       tnf_padded)

write.table(upset_consensus_pathways, paste0(results, "../stable_consensus_vs_pathways.txt"), 
            sep = "\t", row.names = FALSE, na = "", col.names = FALSE, quote = FALSE)

upset_consensus_genes <- data.frame(nodes_padded_genes, 
                                    julio_padded, 
                                    ran_padded, 
                                    tsoi_padded)

write.table(upset_consensus_genes, paste0(results, "../stable_consensus_vs_genelists.txt"), 
            sep = "\t", row.names = FALSE, na = "", col.names = FALSE, quote = FALSE)

#MAGMA is omitted
#common_magma <- !(nodes$name %in% magma$Gene)

common_julio_psoriasis_genes <- nodes$name %in% julio_psoriasis_genes$gene
common_ran_psoriasis_genes <- nodes$name %in% ran_psoriasis_genes$gene
common_tsoi_psoriasis_genes <- nodes$name %in% tsoi_psoriasis_genes$`Closest Gene`

common_il17_pathway <- nodes$name %in% il17_signaling_pathway$gene
common_il23_pathway <- nodes$name %in% il23_signaling_pathway$gene
common_tnf_pathway <- nodes$name %in% tnf_signaling_pathway$gene

# CONSENSUS PIES

consensus_pie_methods <- ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_scatterpie(data = nodes, aes(x = x, y = y, r=0.009), cols = methods) +
  geom_label_repel(data = nodes, 
                   aes(x = x, y = y, label = name), 
                   size = 4, point.padding = 1, max.time = 5, max.iter = 20000, 
                   max.overlaps = 30) +
  coord_fixed() +
  theme_blank() +
  labs(fill = 'Method') +
  scale_fill_manual(values = method_palette) +
  guides(color = "none") +
  theme(legend.position = 'left',
        legend.text = element_text(size = 20, vjust = .2),
        legend.title = element_text(size = 22, vjust = 1)) +
  theme_transparent

consensus_pie_methods

consensus_pie_gene_lists <- ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_scatterpie(data = nodes, aes(x = x, y = y), cols = gen_list) +
  geom_label_repel(data = nodes, 
                   aes(x = x, y = y, label = name), 
                   size = 3, point.padding = 1, max.time = 5, max.iter = 20000, 
                   max.overlaps = 20) +
  coord_fixed() +
  theme_blank() +
  labs(fill = 'Gen list') +
  scale_fill_manual(values = gen_list_palette) +
  guides(color = "none") +
  theme(legend.position = 'left',
        legend.text = element_text(size = 20, vjust = .2),
        legend.title = element_text(size = 22, vjust = 1)) +
  theme_transparent

consensus_pie_gene_lists

consensus_pie_pathways <- ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_scatterpie(data = nodes, aes(x = x, y = y), cols = pathways) +
  geom_label_repel(data = nodes, 
                   aes(x = x, y = y, label = name), 
                   size = 3, point.padding = 1, max.time = 5, max.iter = 20000, 
                   max.overlaps = 20) +
  coord_fixed() +
  theme_blank() +
  labs(fill = 'Signaling pathways') +
  scale_fill_manual(values = pathways_palette) +
  guides(color = "none") +
  theme(legend.position = 'left',
        legend.text = element_text(size = 20, vjust = .2),
        legend.title = element_text(size = 22, vjust = 1)) +
  theme_transparent

consensus_pie_pathways

# CONSENSUS NAMES

# MAGMA omitted

# Gene lists

consensus_names_julio = ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_point(data = filter(nodes, !common_julio_psoriasis_genes), 
             aes(x = x, y = y), 
             size = 3, color = 'black') +
  geom_label_repel(data = filter(nodes, !common_julio_psoriasis_genes), 
                   aes(x = x, y = y, label = name), 
                   size = 3, fill = 'white', color = 'black') +
  geom_point(data = filter(nodes, common_julio_psoriasis_genes), 
             aes(x = x, y = y), 
             size = 3, fill = 'white', color = '#c51b7d') +
  geom_label_repel(data = filter(nodes, common_julio_psoriasis_genes), 
                   aes(x = x, y = y, label = name), 
                   size = 3, fill = 'white', color = '#c51b7d') +
  coord_fixed() +
  theme_blank() +
  labs(title="Genes in common with consensus from Julio's list") +
  theme(text = element_text(size = 40), legend.position = 'none', 
        plot.title = element_text(color = '#c51b7d', size = 20, hjust = 0.5))

consensus_names_julio

consensus_names_ran = ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_point(data = filter(nodes, !common_ran_psoriasis_genes), 
             aes(x = x, y = y), 
             size = 3, color = 'black') +
  geom_label_repel(data = filter(nodes, !common_ran_psoriasis_genes), 
                   aes(x = x, y = y, label = name), 
                   size = 3, fill = 'white', color = 'black') +
  geom_point(data = filter(nodes, common_ran_psoriasis_genes), 
             aes(x = x, y = y), 
             size = 3, fill = 'white', color = '#1f78b4') +
  geom_label_repel(data = filter(nodes, common_ran_psoriasis_genes), 
                   aes(x = x, y = y, label = name), 
                   size = 3, fill = 'white', color = '#1f78b4') +
  coord_fixed() +
  theme_blank() +
  labs(title="Genes in common with consensus from Ran2019's list") +
  theme(text = element_text(size = 40), legend.position = 'none', 
        plot.title = element_text(color = '#1f78b4', size = 20, hjust = 0.5))

consensus_names_ran

consensus_names_tsoi = ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_point(data = filter(nodes, !common_tsoi_psoriasis_genes), 
             aes(x = x, y = y), 
             size = 3, color = 'black') +
  geom_label_repel(data = filter(nodes, !common_tsoi_psoriasis_genes), 
                   aes(x = x, y = y, label = name), 
                   size = 3, fill = 'white', color = 'black') +
  geom_point(data = filter(nodes, common_tsoi_psoriasis_genes), 
             aes(x = x, y = y), 
             size = 3, fill = 'white', color = '#33a02c') +
  geom_label_repel(data = filter(nodes, common_tsoi_psoriasis_genes), 
                   aes(x = x, y = y, label = name), 
                   size = 3, fill = 'white', color = '#33a02c') +
  coord_fixed() +
  theme_blank() +
  labs(title="Genes in common with consensus from Tsoi2017's list") +
  theme(text = element_text(size = 40), legend.position = 'none', 
        plot.title = element_text(color = '#33a02c', size = 20, hjust = 0.5))

consensus_names_tsoi

# Pathways

consensus_names_il17 = ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_point(data = filter(nodes, !common_il17_pathway), 
             aes(x = x, y = y), 
             size = 3, color = 'black') +
  geom_label_repel(data = filter(nodes, !common_il17_pathway), 
                   aes(x = x, y = y, label = name), 
                   size = 3, fill = 'white', color = 'black') +
  geom_point(data = filter(nodes, common_il17_pathway), 
             aes(x = x, y = y), 
             size = 3, fill = 'white', color = '#7fbc41') +
  geom_label_repel(data = filter(nodes, common_il17_pathway), 
                   aes(x = x, y = y, label = name), 
                   size = 3, fill = 'white', color = '#7fbc41') +
  coord_fixed() +
  theme_blank() +
  labs(title="Genes in common with consensus from IL17 signaling pathway") +
  theme(text = element_text(size = 40), legend.position = 'none', 
        plot.title = element_text(color = '#7fbc41', size = 20, hjust = 0.5))

consensus_names_il17

consensus_names_il23 = ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_point(data = filter(nodes, !common_il23_pathway), 
             aes(x = x, y = y), 
             size = 3, color = 'black') +
  geom_label_repel(data = filter(nodes, !common_il23_pathway), 
                   aes(x = x, y = y, label = name), 
                   size = 3, fill = 'white', color = 'black') +
  geom_point(data = filter(nodes, common_il23_pathway), 
             aes(x = x, y = y), 
             size = 3, fill = 'white', color = '#8c510a') +
  geom_label_repel(data = filter(nodes, common_il23_pathway), 
                   aes(x = x, y = y, label = name), 
                   size = 3, fill = 'white', color = '#8c510a') +
  coord_fixed() +
  theme_blank() +
  labs(title="Genes in common with consensus from IL23 signaling pathway") +
  theme(text = element_text(size = 40), legend.position = 'none', 
        plot.title = element_text(color = '#8c510a', size = 20, hjust = 0.5))

consensus_names_il23

consensus_names_tnf = ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_point(data = filter(nodes, !common_tnf_pathway), 
             aes(x = x, y = y), 
             size = 3, color = 'black') +
  geom_label_repel(data = filter(nodes, !common_tnf_pathway), 
                   aes(x = x, y = y, label = name), 
                   size = 3, fill = 'white', color = 'black') +
  geom_point(data = filter(nodes, common_tnf_pathway), 
             aes(x = x, y = y), 
             size = 3, fill = 'white', color = '#762a83') +
  geom_label_repel(data = filter(nodes, common_tnf_pathway), 
                   aes(x = x, y = y, label = name), 
                   size = 3, fill = 'white', color = '#762a83') +
  coord_fixed() +
  theme_blank() +
  labs(title="Genes in common with consensus from TNF signaling pathway") +
  theme(text = element_text(size = 40), legend.position = 'none', 
        plot.title = element_text(color = '#762a83', size = 20, hjust = 0.5))

consensus_names_tnf

# MAGMA omitted
#ggsave(paste0(results,'../figures/consensus_names_magma.pdf'), consensus_names_magma, width=15, height=20, bg = "transparent")
ggsave(paste0(results,'../figures/stable_consensus_pie_methods.pdf'), consensus_pie_methods, width=15, height=20, bg = "transparent")

consensus_names_gene_lists <- list(consensus_names_julio, consensus_names_ran, consensus_names_tsoi, consensus_pie_gene_lists)

pdf(paste0(results,'../figures/stable_consensus_names_gene_lists.pdf'), width=15, height=15, bg = "transparent")
for (plot in consensus_names_gene_lists) {
  print(plot)
}
dev.off()

consensus_names_pathways <- list(consensus_names_il23, consensus_names_il17, consensus_names_tnf, consensus_pie_pathways)
pdf(paste0(results,'../figures/stable_consensus_names_pathways.pdf'), width=15, height=15, bg = "transparent")
for (plot in consensus_names_pathways) {
  print(plot)
}
dev.off()
