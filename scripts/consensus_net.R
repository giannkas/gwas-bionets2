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

# To install RCy3, you need BiocManager.
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("RCy3")

theme_transparent <- theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
                           plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                           panel.grid.major = element_blank(), # get rid of major grid
                           panel.grid.minor = element_blank(), # get rid of minor grid
                           legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
                           legend.box.background = element_rect(fill = "transparent", color = NA) # get rid of legend panel bg
)

theme_set(theme_cowplot() + theme_transparent)
results <- '../results/'

ppi <- read_tsv('../data/network.tsv', col_types = 'cc') %>%
  select(`Official Symbol Interactor A`, `Official Symbol Interactor B`) %>%
  graph_from_data_frame(directed = FALSE) %>%
  as_tbl_graph

ppi_df <- as_tibble(ppi)

methods <- c('Heinz','HotNet2', 'SigMod', 'cGWAS', 'MAGMABH', 'dmGWAS', 'LEAN', 'MAGMABHcGWAS')

method_palette <- c('Heinz' = '#73471f',
                    'HotNet2' = '#377eb8', 'SigMod' = '#f4888a', 'All' = 'gray50',
                    'Consensus' = 'black', 'MAGMABH' = '#abdda4', 'cGWAS' = '#fdae61', 'dmGWAS' = '#80cdc1', 'LEAN' = '#4d9221')

genesrc_palette <- c('MAGMABH' = '#abdda4', 'cGWAS' = '#fdae61', 'MAGMABH and cGWAS' = '#d3c783', 'Methods only' = '#ffffff')

labs <- c('all_snps' = 'All', 'consensus' = 'Consensus',
          'heinz' = 'heinz', 'hotnet2' = 'HotNet2', 'sigmod' = 'SigMod', 'magmabh' = 'MAGMABH', 'cgwas' = 'cGWAS', 'dmgwas' = 'dmGWAS', 'lean' = 'LEAN',
          'magmabh_and_cgwas' = 'MAGMABHcGWAS')

# TOPOLOGY PANEL

plot_network <- function(net, selected_nodes) {
  
  graph <- activate(net, nodes) %>%
    filter(name %in% selected_nodes)
  class(graph) <- c('igraph','tbl_graph')
  
  ggnetwork(graph) %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = '#3282b8') +
    geom_nodes(color = '#233142', size = 4) +
    #geom_label(aes(label = name),
    #                 size = 2) +
    ggtitle(' ') +
    theme_blank() + theme_transparent + 
    theme(title = element_text(hjust = 1, face = 'bold'))
}

# heinz
heinz <- read_tsv(paste0(results,'heinz/selected_genes_split_1.heinz.txt'), col_types = 'c')$gene
heinz_plt <- plot_network(ppi, heinz) 
heinz_plt_square <- heinz_plt + coord_fixed(ratio = 1)# Add coord_fixed(ratio = 1) to make the plot square

# hotnet2
hotnet2 <- read_tsv(paste0(results,'hotnet2/selected_genes_split_1.hotnet2.tsv'), col_types = 'cc')$gene
hotnet2_plt <- plot_network(ppi, hotnet2)
hotnet2_plt_square <- hotnet2_plt + coord_fixed(ratio = 1)

# sigmod
sigmod <- read_tsv(paste0(results,'sigmod/selected_genes_split_1.sigmod.txt'), col_types = 'c')$gene
sigmod_plt <- plot_network(ppi, sigmod)
sigmod_plt_square <- sigmod_plt + coord_fixed(ratio = 1)

# dmgwas
dmgwas <- read_tsv(paste0(results,'dmgwas/selected_genes_split_1.dmgwas.txt'), col_types = 'c')$gene
dmgwas_plt <- plot_network(ppi, dmgwas)
dmgwas_plt_square <- dmgwas_plt + coord_fixed(ratio = 1)

# lean
lean <- read_tsv(paste0(results,'lean/selected_genes_split_1.lean.txt'), col_types = 'c')$gene
lean_plt <- plot_network(ppi, lean)
lean_plt_square <- lean_plt + coord_fixed(ratio = 1)

# Classical GWAS
cgwas_genes <- read_tsv(paste0(results,'pso.genes.annot_gene_ids_names.tsv'), col_types = 'c')$Gene %>% unique()
cgwas_plt <- plot_network(ppi, cgwas_genes)
cgwas_plt_square <- plot_network(ppi, cgwas_genes) + coord_fixed(ratio = 1)

# Magma output after Benjamini Hochberg correction
magmabh <- read_tsv(paste0(results,'pso.scores.genes.out_converted_bhcorrection.tsv'), 
                    col_types = 'iiciiiiiidnn')$Gene %>% unique()
magmabh_plt <- plot_network(ppi, magmabh)
magmabh_plt_square <- plot_network(ppi, magmabh) + coord_fixed(ratio = 1)

# Magma output together with classical gwas annotation output
magmabh_cgwas <- union(magmabh, cgwas_genes)
magmabh_cgwas_plt <- plot_network(ppi, magmabh_cgwas)
magmabh_cgwas_plt_square <- plot_network(ppi, magmabh_cgwas) + coord_fixed(ratio = 1)

# Magma output intersected with classical gwas annotation output
magmabh_and_cgwas <- intersect(magmabh, cgwas_genes)
magmabh_and_cgwas_plt <- plot_network(ppi, magmabh_and_cgwas)
magmabh_and_cgwas_plt_square <- plot_network(ppi, magmabh_and_cgwas) + coord_fixed(ratio = 1)

topology_panel <- plot_grid(heinz_plt, hotnet2_plt, sigmod_plt, dmgwas_plt, lean_plt, magmabh_plt, cgwas_plt, magmabh_cgwas_plt, 
                            labels = c('\nHeinz','\nHotNet2','\nSigMod', 'dmGWAS', 'LEAN', 'MAGMABH', '\ncGWAS', 'MAGMABH+cGWAS'),
                            label_size = 30,
                            nrow = 4,
                            ncol = 2)

# networks for each of the methods

ggsave(paste0(results,'figures/networks_foreach_method.pdf'), topology_panel, width=30, height=60, limitsize = FALSE, bg = "transparent")

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

# MAGMA_BH + cGWAS (union)

g_magmabh_cgwas <- activate(ppi, nodes) %>%
  filter(name %in% magmabh_cgwas) %>%
  mutate(magmabh_cgwas = TRUE)
class(g_magmabh_cgwas) <- c('igraph','tbl_graph')

# MAGMA_BH and cGWAS (intersection)

g_magmabh_and_cgwas <- activate(ppi, nodes) %>%
  filter(name %in% magmabh_and_cgwas) %>%
  mutate(magmabh_and_cgwas = TRUE)
class(g_magmabh_and_cgwas) <- c('igraph','tbl_graph')


# Consensus

options(repr.plot.width=25, repr.plot.height=15)
consensus <- graph_join(g_heinz, g_hotnet2, by = c('name')) %>%
  graph_join(g_sigmod, by = c('name')) %>%
  graph_join(g_dmgwas, by = c('name')) %>%
  graph_join(g_lean, by = c('name')) %>%
  graph_join(g_magmabh_and_cgwas, by = c('name')) %>%
  graph_join(g_cgwas, by = c('name')) %>%
  graph_join(g_magmabh, by = c('name')) %>%
  to_undirected %>%
  mutate(heinz = ifelse(is.na(heinz), FALSE, heinz),
         hotnet2 = ifelse(is.na(hotnet2), FALSE, hotnet2),
         sigmod = ifelse(is.na(sigmod), FALSE, sigmod),
         dmgwas = ifelse(is.na(dmgwas), FALSE, dmgwas),
         lean = ifelse(is.na(lean), FALSE, lean),
         magmabh_and_cgwas = ifelse(is.na(magmabh_and_cgwas), FALSE, magmabh_and_cgwas),
         cgwas_genes = ifelse(is.na(cgwas_genes), FALSE, cgwas_genes),
         magmabh = ifelse(is.na(magmabh), FALSE, magmabh),
         num_methods = rowSums(cbind(heinz, hotnet2, sigmod, dmgwas, lean))) %>%
  filter(num_methods >= 3)

as_tibble(consensus) %>%
  rename(gene = name) %>%
  write_tsv(paste0(results,'consensus.tsv'))

consensus_gg <- consensus %>%
  ggnetwork %>%
  mutate(name = as.character(name))


nodes <- mutate(consensus_gg, Heinz = as.numeric(heinz),
                HotNet2 = as.numeric(hotnet2),
                SigMod = as.numeric(sigmod),
                dmGWAS = as.numeric(dmgwas),
                LEAN = as.numeric(lean),
                cGWAS = as.numeric(cgwas_genes),
                MAGMABH = as.numeric(magmabh),
                MAGMABHcGWAS = as.numeric(magmabh_and_cgwas),
                ) %>%
  filter(xend == x & yend == y) %>% 
  select(x, y, name, all_of(methods)) %>%
  unique

nodes <- nodes %>%
  mutate(label_group = case_when(
                        MAGMABHcGWAS == 1 ~ 'MAGMABH and cGWAS',
                        MAGMABH == 1 ~ 'MAGMABH',
                        cGWAS == 1 ~ 'cGWAS',
                        TRUE ~ 'Methods only')
  )

edges <- filter(consensus_gg, xend != x | yend != y)


# CONSENSUS PIES

consensus_pie_methods <- ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_scatterpie(data = nodes, aes(x = x, y = y, r=0.007), cols = setdiff(methods, c('MAGMABH', 'cGWAS', 'MAGMABHcGWAS'))) +
  scale_fill_manual(name = "Methods", values = method_palette) +
  ggnewscale::new_scale_fill() +
  geom_label_repel(data = nodes, 
                   aes(x = x, y = y, label = name, fill = label_group),
                   size = 3, point.padding = 1, max.time = 5, max.iter = 20000, 
                   max.overlaps = 50) +
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


ggsave(paste0(results,'figures/consensus_pie_>=3methods_outof5_magmabh_cgwas.pdf'), consensus_pie_methods, width=14, height=10, bg = "transparent")


