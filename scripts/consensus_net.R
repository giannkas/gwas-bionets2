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
library(ggnewscale)

theme_transparent <- theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
                           plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                           panel.grid.major = element_blank(), # get rid of major grid
                           panel.grid.minor = element_blank(), # get rid of minor grid
                           legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
                           legend.box.background = element_rect(fill = "transparent", color = NA) # get rid of legend panel bg
)

theme_set(theme_cowplot() + theme_transparent)
results <- '../results/psoriasis/magma/magma_output_50kb_window_snp-array/'
results_gwas <- '../results/psoriasis/magma/magma_output_50kb_window_snp-array_significant_snps/'
results_dsl2 <- '../results/psoriasis/magma/magma_output_50kb_window_snp-array_dsl2/'
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

ppi <- read_tsv('../data/biogrid_ppi_interactorsAB.tsv', col_types = 'cc') %>%
  select(`Official Symbol Interactor A`, `Official Symbol Interactor B`) %>%
  graph_from_data_frame(directed = FALSE) %>%
  as_tbl_graph

ppi_df <- as_tibble(ppi)

# methods <- c('Heinz','HotNet2', 'SigMod', 'GWAS', 'dmGWAS', 'LEAN')
methods <- c('Heinz','HotNet2', 'SigMod', 'cGWAS', 'MAGMABH', 'dmGWAS', 'LEAN', 'MAGMABHcGWAS')
pathways <- c('IL17','IL23', 'TNF')
gen_list <- c('Julio', 'Ran2019', 'Tsoi2017')

method_palette <- c('Heinz' = '#73471f',
                    'HotNet2' = '#377eb8', 'SigMod' = '#f4888a', 'All' = 'gray50',
                    'Consensus' = 'black', 'MAGMABH' = '#abdda4', 'cGWAS' = '#fdae61', 'dmGWAS' = '#80cdc1', 'LEAN' = '#4d9221')

genesrc_palette <- c('MAGMABH' = '#abdda4', 'cGWAS' = '#fdae61', 'MAGMABH and cGWAS' = '#d3c783', 'Methods only' = '#ffffff')
#genesrc_palette <- c('GWAS' = '#e2c7c5', 'Methods' = '#ffffff')

# method_palette <- c('Heinz' = '#73471f',
#                     'HotNet2' = '#377eb8', 'SigMod' = '#f4888a', 'All' = 'gray50',
#                     'Consensus' = 'black', 'MAGMA' = 'gray50', 'dmGWAS' = '#80cdc1', 'LEAN' = '#4d9221')
#'GWAS' = '#e2c7c5'
pathways_palette <- c('IL17' = '#489ea3',
                      'IL23' = '#ff6f3c', 'TNF' = '#377eb8')
gen_list_palette <- c('Julio' = '#fb9a99',
                      'Ran2019' = '#b15928', 'Tsoi2017' = '#6a3d9a')


labs <- c('all_snps' = 'All', 'consensus' = 'Consensus',
          'heinz' = 'heinz', 'hotnet2' = 'HotNet2', 'sigmod' = 'SigMod', 'magmabh' = 'MAGMABH', 'cgwas' = 'cGWAS', 'dmgwas' = 'dmGWAS', 'lean' = 'LEAN',
          'magmabh_cgwas' = 'MAGMABHcGWAS')

# TOPOLOGY PANEL

plot_network <- function(net, selected_nodes) {
  
  graph <- activate(net, nodes) %>%
    filter(name %in% selected_nodes)
  class(graph) <- c('igraph','tbl_graph')
  
  ggnetwork(graph) %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = '#3282b8') +
    geom_nodes(color = '#233142', size = 5) +
    ggtitle(' ') +
    theme_blank() + theme_transparent + 
    theme(title = element_text(hjust = 1, face = 'bold'))
}

# inter_hesi <- intersect(heinz, sigmod)
# inter_heho <- intersect(heinz, hotnet2)
# inter_hosi <- intersect(hotnet2, sigmod)
# combined_intersections <- unique(c(inter_hesi, inter_heho, inter_hosi))
# heinz
heinz <- read_tsv(paste0(results,'heinz/selected_genes.heinz.txt'), col_types = 'c')$gene
heinz_plt <- plot_network(ppi, heinz) 
heinz_plt_square <- heinz_plt + coord_fixed(ratio = 1)# Add coord_fixed(ratio = 1) to make the plot square

# hotnet2
hotnet2 <- read_tsv(paste0(results,'hotnet2/selected_genes.hotnet2.tsv'), col_types = 'cc')$gene
hotnet2_plt <- plot_network(ppi, hotnet2)
hotnet2_plt_square <- hotnet2_plt + coord_fixed(ratio = 1)

# sigmod
sigmod <- read_tsv(paste0(results,'sigmod/selected_genes.sigmod.txt'), col_types = 'c')$gene
sigmod_plt <- plot_network(ppi, sigmod)
sigmod_plt_square <- sigmod_plt + coord_fixed(ratio = 1)

# dmgwas
dmgwas <- read_tsv(paste0(results_dsl2,'dmgwas/selected_genes.dmgwas.txt'), col_types = 'c')$gene
dmgwas_plt <- plot_network(ppi, dmgwas)
dmgwas_plt_square <- dmgwas_plt + coord_fixed(ratio = 1)

lean <- read_tsv(paste0(results_dsl2,'lean/selected_genes.lean.txt'), col_types = 'c')$gene
lean_plt <- plot_network(ppi, lean)
lean_plt_square <- lean_plt + coord_fixed(ratio = 1)

# Classical GWAS
cgwas_genes <- read_tsv(paste0(results_gwas,'pso.genes.annot_gene_ids_names.tsv'), col_types = 'c')$Gene
cgwas_plt <- plot_network(ppi, cgwas_genes)
cgwas_plt_square <- plot_network(ppi, cgwas_genes) + coord_fixed(ratio = 1)

# Magma output after Benjamini Hochberg correction
magmabh <- read_tsv(paste0(results,'pso.scores.genes.out_converted_bhcorrection.tsv'), 
                    col_types = 'iiciiiiiidnn')$Gene
magmabh_plt <- plot_network(ppi, magmabh)
magmabh_plt_square <- plot_network(ppi, magmabh) + coord_fixed(ratio = 1)

magmabh_cgwas <- union(magmabh, cgwas_genes)
magmabh_cgwas_plt <- plot_network(ppi, magmabh_cgwas)
magmabh_cgwas_plt_square <- plot_network(ppi, magmabh_cgwas) + coord_fixed(ratio = 1)

topology_panel <- plot_grid(heinz_plt, hotnet2_plt, sigmod_plt, dmgwas_plt, lean_plt, magmabh_plt, cgwas_plt, magmabh_cgwas_plt, 
                            labels = c('\nHeinz','\nHotNet2','\nSigMod', 'dmGWAS', 'LEAN', 'MAGMABH', '\ncGWAS', 'MAGMABH+cGWAS'),
                            label_size = 30,
                            nrow = 4,
                            ncol = 2)

f1 <- plot_grid(lean_plt_square,
                labels = c('LEAN'))
#  labels = c("A"),
#  label_size = 60,
#  ncol = 2)

#options(repr.plot.width=15, repr.plot.height=15)
f1

ggsave(paste0(results_dsl2,'figures/networks_foreach_method.pdf'), topology_panel, width=30, height=60, limitsize = FALSE, bg = "transparent")

  ### Consensus ###

# Heinz

g_heinz <- activate(ppi, nodes) %>%
  filter(name %in% heinz) %>%
  mutate(heinz = TRUE)
class(g_heinz) <- c('igraph','tbl_graph')

# HotNet2

g_hotnet2 <- activate(ppi, nodes) %>%
  filter(name %in% hotnet2) %>%
  mutate(hotnet2 = TRUE)
class(g_hotnet2) <- c('igraph','tbl_graph')

# SigMod

g_sigmod <- activate(ppi, nodes) %>%
  filter(name %in% sigmod) %>%
  mutate(sigmod = TRUE)
class(g_sigmod) <- c('igraph','tbl_graph')

# dmGWAS

g_dmgwas <- activate(ppi, nodes) %>%
  filter(name %in% dmgwas) %>%
  mutate(dmgwas = TRUE)
class(g_dmgwas) <- c('igraph','tbl_graph')

# LEAN

g_lean <- activate(ppi, nodes) %>%
  filter(name %in% lean) %>%
  mutate(lean = TRUE)
class(g_lean) <- c('igraph','tbl_graph')

# cGWAS

g_cgwas <- activate(ppi, nodes) %>%
  filter(name %in% cgwas_genes) %>%
  mutate(cgwas_genes = TRUE)
class(g_cgwas) <- c('igraph','tbl_graph')

# MAGMA_BH

g_magmabh <- activate(ppi, nodes) %>%
  filter(name %in% magmabh) %>%
  mutate(magmabh = TRUE)
class(g_magmabh) <- c('igraph','tbl_graph')

# MAGMA_BH and cGWAS

g_magmabh_cgwas <- activate(ppi, nodes) %>%
  filter(name %in% magmabh_cgwas) %>%
  mutate(magmabh_cgwas = TRUE)
class(g_magmabh_cgwas) <- c('igraph','tbl_graph')


# Consensus

options(repr.plot.width=25, repr.plot.height=15)
consensus <- graph_join(g_heinz, g_hotnet2, by = c('name')) %>%
  graph_join(g_sigmod, by = c('name')) %>%
  graph_join(g_dmgwas, by = c('name')) %>%
  graph_join(g_lean, by = c('name')) %>%
  graph_join(g_magmabh_cgwas, by = c('name')) %>%
  graph_join(g_cgwas, by = c('name')) %>%
  graph_join(g_magmabh, by = c('name')) %>%
  to_undirected %>%
  mutate(heinz = ifelse(is.na(heinz), FALSE, heinz),
         hotnet2 = ifelse(is.na(hotnet2), FALSE, hotnet2),
         sigmod = ifelse(is.na(sigmod), FALSE, sigmod),
         dmgwas = ifelse(is.na(dmgwas), FALSE, dmgwas),
         lean = ifelse(is.na(lean), FALSE, lean),
         magmabh_cgwas = ifelse(is.na(magmabh_cgwas), FALSE, magmabh_cgwas),
         cgwas_genes = ifelse(is.na(cgwas_genes), FALSE, cgwas_genes),
         magmabh = ifelse(is.na(magmabh), FALSE, magmabh),
         num_methods = rowSums(cbind(heinz, hotnet2, sigmod, dmgwas, lean))) %>%
  filter(num_methods >= 3)

as_tibble(consensus) %>%
  rename(gene = name) %>%
  write_tsv(paste0(results_dsl2,'consensus.tsv'))

consensus_gg <- consensus %>%
  ggnetwork %>%
  mutate(name = as.character(name))


nodes <- mutate(consensus_gg, Heinz = as.numeric(heinz),
                HotNet2 = as.numeric(hotnet2),
                SigMod = as.numeric(sigmod),
                dmGWAS = as.numeric(dmgwas),
                LEAN = as.numeric(lean),
                MAGMABHcGWAS = as.numeric(magmabh_cgwas),
                cGWAS = as.numeric(cgwas_genes),
                MAGMABH = as.numeric(magmabh),
                Julio = as.numeric(name %in% julio_psoriasis_genes$gene),
                Ran2019 = as.numeric(name %in% ran_psoriasis_genes$gene),
                Tsoi2017 =  as.numeric(name %in% tsoi_psoriasis_genes$`Closest Gene`),
                IL17 = as.numeric(name %in% il17_signaling_pathway$gene),
                IL23 = as.numeric(name %in% il23_signaling_pathway$gene),
                TNF =  as.numeric(name %in% tnf_signaling_pathway$gene)) %>%
  filter(xend == x & yend == y) %>% 
  select(x, y, name, all_of(methods), all_of(gen_list), all_of(pathways)) %>%
  unique

#genesrc_palette <- c('MAGMABH' = '#abdda4', 'cGWAS' = '#fdae61', 'MAGMABH and cGWAS' = '#d3c783', 'Methods only' = '#ffffff')

nodes <- nodes %>%
  mutate(label_group = case_when(
                        MAGMABHcGWAS == 1 ~ 'MAGMABH and cGWAS',
                        MAGMABH == 1 ~ 'MAGMABH',
                        cGWAS == 1 ~ 'cGWAS',
                        TRUE ~ 'Methods only')
  )

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

write.table(upset_consensus_pathways, paste0(results, "consensus_vs_pathways.txt"), 
            sep = "\t", row.names = FALSE, na = "", col.names = FALSE, quote = FALSE)

upset_consensus_genes <- data.frame(nodes_padded_genes, 
                                    julio_padded, 
                                    ran_padded, 
                                    tsoi_padded)

write.table(upset_consensus_genes, paste0(results, "consensus_vs_genelists.txt"), 
            sep = "\t", row.names = FALSE, na = "", col.names = FALSE, quote = FALSE)

#highlight <- consensus_gg$num_methods == 3
common_magma <- !(nodes$name %in% magma$Gene)

common_julio_psoriasis_genes <- nodes$name %in% julio_psoriasis_genes$gene
common_ran_psoriasis_genes <- nodes$name %in% ran_psoriasis_genes$gene
common_tsoi_psoriasis_genes <- nodes$name %in% tsoi_psoriasis_genes$`Closest Gene`
common_all_gen_lists <- common_julio_psoriasis_genes | common_ran_psoriasis_genes | common_tsoi_psoriasis_genes

common_il17_pathway <- nodes$name %in% il17_signaling_pathway$gene
common_il23_pathway <- nodes$name %in% il23_signaling_pathway$gene
common_tnf_pathway <- nodes$name %in% tnf_signaling_pathway$gene
common_all_pathways <- common_il17_pathway | common_il23_pathway | common_tnf_pathway

# CONSENSUS PIES

# methods <- c('Heinz','HotNet2', 'SigMod', 'cGWAS', 'MAGMABH', 'dmGWAS', 'LEAN', 'MAGMABHcGWAS')

consensus_pie_methods <- ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_scatterpie(data = nodes, aes(x = x, y = y, r=0.009), cols = setdiff(methods, c('MAGMABH', 'cGWAS', 'MAGMABHcGWAS'))) +
  scale_fill_manual(name = "Methods", values = method_palette) +
  ggnewscale::new_scale_fill() +
  geom_label_repel(data = nodes, 
                   aes(x = x, y = y, label = name, fill = label_group),
                   size = 3, point.padding = 1, max.time = 5, max.iter = 20000, 
                   max.overlaps = 30) +
  scale_fill_manual(name = "Gene source", values = genesrc_palette) +
  coord_fixed() +
  theme_blank() +
  guides(color = "none", fill = guide_legend(
    override.aes = aes(label = "  ")
  )) +
  theme(legend.position = 'left',
        legend.text = element_text(size = 20, vjust = .2),
        legend.title = element_text(size = 22, vjust = 1)) +
  theme_transparent

consensus_pie_methods

consensus_pie_gene_lists <- ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_point(data = filter(nodes, !common_all_gen_lists), 
             aes(x = x, y = y), 
             size = 1, color = '#393e46') +
  geom_scatterpie(data = nodes, aes(x = x, y = y, r=0.011), cols = gen_list) +
  geom_label_repel(data = filter(nodes, common_all_gen_lists), 
                   aes(x = x, y = y, label = name), 
                   size = 6, point.padding = 1, max.time = 5, max.iter = 20000, 
                   max.overlaps = 30) +
  coord_fixed() +
  theme_blank() +
  labs(fill = 'Gen list') +
  scale_fill_manual(values = gen_list_palette) +
  guides(color = "none") +
  theme(legend.position = 'none',
        legend.text = element_text(size = 20, vjust = .2),
        legend.title = element_text(size = 22, vjust = 1)) +
  theme_transparent

consensus_pie_gene_lists

consensus_pie_pathways <- ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_point(data = filter(nodes, !common_all_pathways), 
             aes(x = x, y = y), 
             size = 1, color = '#393e46') +
  geom_scatterpie(data = nodes, aes(x = x, y = y, r=0.011), cols = pathways) +
  geom_label_repel(data = filter(nodes, common_all_pathways), 
                   aes(x = x, y = y, label = name), 
                   size = 6, point.padding = 1, max.time = 5, max.iter = 20000, 
                   max.overlaps = 50) +
  coord_fixed() +
  theme_blank() +
  labs(fill = 'Signaling pathways') +
  scale_fill_manual(values = pathways_palette) +
  guides(color = "none") +
  theme(legend.position = 'none',
        legend.text = element_text(size = 20, vjust = .2),
        legend.title = element_text(size = 22, vjust = 1)) +
  theme_transparent

consensus_pie_pathways

# CONSENSUS NETWORK

consensus_network = ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#3282b8', linewidth = 1) +
  geom_point(data = nodes, 
             aes(x = x, y = y), 
             size = 5, color = '#233142') +
  # geom_text_repel(data = nodes, 
  #                  aes(x = x, y = y, label = name), 
  #                  size = 4, color = 'black', max.overlaps = 30) +
  coord_fixed() +
  theme_blank() 
# +
#   labs(title="Genes selected by the methods") +
#   theme(text = element_text(size = 40), legend.position = 'none',
#         plot.title = element_text(color = '#233142', size = 20, hjust = 0.5))

consensus_network



consensus_names_magma = ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_point(data = filter(nodes, !common_magma), 
             aes(x = x, y = y), 
             size = 3, color = 'black') +
  geom_label_repel(data = filter(nodes, !common_magma), 
                   aes(x = x, y = y, label = name), 
                   size = 3, fill = 'white', color = 'black') +
  geom_point(data = filter(nodes, common_magma), 
             aes(x = x, y = y), 
             size = 3, fill = 'white', color = 'brown') +
  geom_label_repel(data = filter(nodes, common_magma), 
                   aes(x = x, y = y, label = name), 
                   size = 3, fill = 'white', color = 'brown') +
  coord_fixed() +
  theme_blank() +
  labs(title="Genes selected by the methods but not by MAGMA") +
  theme(text = element_text(size = 40), legend.position = 'none',
        plot.title = element_text(color = 'brown', size = 20, hjust = 0.5))

consensus_names_magma

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

ggsave(paste0(results,'figures/consensus_names_magma.pdf'), consensus_names_magma, width=15, height=20, bg = "transparent")
ggsave(paste0(results_dsl2,'figures/consensus_pie_>=3methods_outof5_magmabh_cgwas.pdf'), consensus_pie_methods, width=14, height=10, bg = "transparent")
ggsave(paste0(results_dsl2,'figures/consensus_pie_>=3methods_outof5_gene_lists.pdf'), consensus_pie_gene_lists, width=10, height=10, bg = "transparent")
ggsave(paste0(results_dsl2,'figures/consensus_pie_>=3methods_outof5_pathways.pdf'), consensus_pie_pathways, width=10, height=10, bg = "transparent")
ggsave(paste0(results,'figures/consensus_network_methods.pdf'), consensus_pie_methods, width=15, height=15, limitsize = FALSE, bg = "transparent")

consensus_names_gene_lists <- list(consensus_names_julio, consensus_names_ran, consensus_names_tsoi, consensus_pie_gene_lists)

pdf(paste0(results,'figures/consensus_names_gene_lists.pdf'), width=15, height=15, bg = "transparent")
for (plot in consensus_names_gene_lists) {
  print(plot)
}
dev.off()

consensus_names_pathways <- list(consensus_names_il23, consensus_names_il17, consensus_names_tnf, consensus_pie_pathways)
pdf(paste0(results,'figures/consensus_names_pathways.pdf'), width=15, height=15, bg = "transparent")
for (plot in consensus_names_pathways) {
  print(plot)
}
dev.off()

#ggsave(paste0(results,'figures/consensus_names.pdf'), consensus_names_magma, consensus_names_psoriasis, width=15, height=20, bg = "transparent")

print(length(V(ppi)))
print(nrow(magmabh))
print(summary(magmabh))
print(vertex_attr_names(ppi))
