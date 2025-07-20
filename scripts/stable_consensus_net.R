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

# To install tidyverse, you may need to install these libraries in Ubuntu:
# apt install libcurl4-openssl-dev libfontconfig1-dev libxml2-dev libharfbuzz-dev libfribidi-dev
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
#results <- '../results/psoriasis/magma/magma_output_50kb_window_wes_ncbi38_tab2gwen/stable_consensus/data_splits/'
#results <- '../results/psoriasis/magma/magma_output_50kb_window_wes/stable_consensus/data_splits/'
results <- '../results/psoriasis/magma/magma_output_50kb_window_snp-array/stable_consensus/data_splits/'
results_dsl2 <- '../results/psoriasis/magma/magma_output_50kb_window_snp-array_dsl2/stable_consensus/data_splits/'
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

methods <- c('Heinz','HotNet2', 'SigMod', 'dmGWAS', 'LEAN')
pathways <- c('IL17','IL23', 'TNF')
gen_list <- c('Julio', 'Ran2019', 'Tsoi2017')

method_palette <- c('Heinz' = '#73471f',
                    'HotNet2' = '#377eb8', 'SigMod' = '#f4888a', 'All' = 'gray50', 
                    'Consensus' = 'black', 'MAGMA' = 'gray50', 'GWAS' = '#e2c7c5', 'dmGWAS' = '#80cdc1', 'LEAN' = '#4d9221')
pathways_palette <- c('IL17' = '#489ea3',
                      'IL23' = '#ff6f3c', 'TNF' = '#377eb8')
gen_list_palette <- c('Julio' = '#fb9a99',
                      'Ran2019' = '#b15928', 'Tsoi2017' = '#6a3d9a')

labs <- c('all_snps' = 'All', 'consensus' = 'Consensus',
          'heinz' = 'heinz', 'hotnet2' = 'HotNet2', 'sigmod' = 'SigMod', 'gwas' = 'GWAS', 'dmgwas' = 'dmGWAS', 'lean' = 'LEAN' )

# TOPOLOGY PANEL

# color some nodes according to the 5 modules found in the gene set of 11 solutions
node_colors_gene_set11 <- c(
  "CDC37" = "#c51b7d",
  "RPS6KA4" = "#c51b7d",
  "DDR1" = "#c51b7d",
  "FKBPL" = "#c51b7d",
  "TUBB" = "#c51b7d",
  "UBD" = "#ff0000",
  "PSMA6" = "#ff0000",
  "PSMB8" = "#ff0000",
  "PSMB9" = "#ff0000",
  "MRPS18B" = "#00ff00",
  "CS" = "#00ff00",
  "MRPL24" = "#00ff00",
  "MRPL4" = "#00ff00",
  "IRF1" = "#00ff00",
  "TRIM40" = "#00ffef",
  "TRIM65" = "#00ffef",
  "IFIH1" = "#00ffef",
  "USP3" = "#00ffef",
  "SNAI1" = "#00ffef",
  "RING1" = "#00ffef",
  "CSNK2B" = "#00ffef",
  "BRD2" = "#00ffef",
  "ABT1" = "#00ffef",
  "NFKBIA" = "#00ffef",
  "REL" = "#00ffef",
  "NFKB1" = "#00ffef",
  "ETS1" = "#00ffef",
  "TNIP1" = "#00ffef",
  "TNFAIP3" = "#00ffef",
  "TNF" = "#00ffef",
  "RNF114" = "#00ffef",
  "SPATA2" = "#00ffef",
  "HLA-A" = "#fdae61",
  "TAP1" = "#fdae61",
  "TAP2" = "#fdae61",
  "TAPBP" = "#fdae61",
  "CD8A" = "#fdae61",
  "HLA-G" = "#fdae61",
  "ERBB3" = "#fdae61",
  "PA2G4" = "#fdae61",
  "HSPA1A" = "#fdae61",
  "HSPA1L" = "#fdae61",
  "VPS52" = "#fdae61",
  "USP8" = "#fdae61",
  "RNF41" = "#fdae61",
  "CLEC16A" = "#fdae61",
  "BAG6" = "#1a9850",
  "NCR3" = "#1a9850",
  "KEAP1" = "#A38658",
  "RAD50" = "#35978f",
  "MDC1" = "#35978f",
  "RBM17" = "#35978f",
  "RNF5" = "#878787",
  "ABHD16A" = "#878787", 
  "FBF1" = "#d6604d",
  "TRIM27" = "#d6604d",
  "RPS26" = "#abd9e9",
  "RPS18" = "#abd9e9",
  "MRPL38" = "#543005"
)

node_names_gene_set8 <- c(
  "STAMBPL1",
  "GTSE1",
  "XPO1",
  "MAPRE1",
  "TUBB",
  "HTRA4",
  "DDR1",
  "FKBPL",
  "NTRK1",
  "RPS6KA4",
  "CDC37",
  "MAP3K14",
  "TRAF3IP2",
  "CHUK",
  "CEP170",
  "SUPT16H",
  "ABT1",
  "KRI1",
  "KRR1",
  "BRD2",
  "CSNK2B",
  "NFKBIA",
  "REL",
  "NFKB1",
  "TRIM27",
  "LYST",
  "RING1",
  "FBF1",
  "ETS1",
  "TNIP1",
  "SPATA2",
  "RIPK3",
  "TNF",
  "TNFAIP3",
  "SNAI1",
  "RNF114",
  "AGO2",
  "TNRC6A",
  "CNOT8",
  "UNK",
  "TRIM65",
  "IFIH1",
  "USP3",
  "TRIM40",
  
  "UBD",
  "PSMG2",
  "PSMB9",
  "PSMB8",
  "PSMA6",
  "CDK3",
  "CDK2",
  "CDC7",
  "CLSPN",
  "DBF4",
  "MSH5",
  "MSH4",
  "ID2",
  "GATA4",
  "JARID2",
  "MTF2",
  
  "TAP2",
  "TAPBP",
  "TAP1",
  "HLA−G",
  "CD8A",
  "HLA−A",
  "CTSB",
  "RHBDD1",
  "RNF5",
  "ABHD16A",
  "S100A10",
  "ERBB3",
  "HSPA1A",
  "HSPA1L",
  "PA2G4",
  "RNF41",
  "USP8",
  "CLEC16A",
  "VPS52",
  "AAMP"
)

node_names_gene_set5 <- c(
  "CS",
  "VARS2",
  "MRPS18B",
  "MRPL24",
  "MRPL4",
  "MRPL38",
  "DDX28",
  "PPIF",
  "XRN2",
  "MRPL51",
  "MRPL52",
  "MRPL14",
  
  "HMGXB4",
  "BPTF",
  "MAZ",
  "KEAP1",
  "STK4",
  "IRF1",
  "IRF9",
  "STAT2"
)

node_names_gene_mitochondrial <- c(
  "CS",
  "MRPS18B",
  "MRPL24",
  "MRPL4",
  "MRPL38",
  "DDX28"
)

node_names_genes_consensus <- c(
  "NFKB1",
  "ETS1",
  "REL",
  "TNFAIP3",
  "ID2",
  "RING1",
  "RIPK3",
  "LSM2",
  "TOPORS",
  "VPS52",
  "RAD50",
  "ABHD16A",
  "CSNK2B",
  "TRIM39",
  "TRIM17",
  "TRIM27",
  "BAG6",
  "RPS26",
  "TNIP1",
  "COG6",
  "CDC37",
  "MRPL24",
  "KRI1",
  "CDK2",
  "ERBB3",
  "PAAF1",
  "PPM1B",
  "FAS",
  "MDC1",
  "HSPA1A",
  "RNF41",
  "KEAP1",
  "CTSB",
  "TNF",
  "HLA-G",
  "GRB7",
  "DIABLO",
  "NFKBIA",
  "RPS6KA4",
  "IRF1",
  "TNRC6A",
  "PSMA6",
  "CD8A",
  "MAP3K14",
  "HLA-A",
  "BRD2",
  "PA2G4",
  "UBE2K",
  "USP8",
  "HLA-C",
  "NOS2",
  "LTA",
  "TAP1",
  "TAPBP",
  "POLD2",
  "TRAF3IP2",
  "UBD",
  "LTB",
  "SNAI1",
  "CCHCR1",
  "SPATA2",
  "MRPL4",
  "PTPN2",
  "TRIM15",
  "SLBP",
  "PPP1R10",
  "FBXO45",
  "TUBB",
  "RNF5",
  "IFIH1",
  "POLD1",
  "FKBPL",
  "TYK2",
  "PSMB5",
  "PSMB9",
  "ACIN1",
  "TRIM31",
  "PSMB8",
  "RNF114",
  "PSMG2",
  "PRCC",
  "DDX39B",
  "CS",
  "MRPL38",
  "MRPS18B",
  "DDX28",
  "HSPA1L",
  "RBM17",
  "RPS18",
  "TRIM65",
  "NCR3",
  "FBF1",
  "PPP1R11",
  "PGBD1",
  "UBLCP1",
  "TAP2",
  "ZKSCAN4",
  "TRIM40",
  "ABT1",
  "DDR1",
  "PPP1R18",
  "LY6G6F",
  "HLA-DRA",
  "HLA-B",
  "HLA-DMB"
)

plot_network <- function(net, selected_nodes, node_colors = NULL, label_genes = NULL) {
  
  graph <- activate(net, nodes) %>%
    filter(name %in% selected_nodes)
  
  graph <- graph %>%
    mutate(color = ifelse(name %in% names(node_colors), 
                          '#233142',
                          #node_colors[name], 
                          '#233142'))
  
  # Ensure the graph is of class igraph and tbl_graph
  class(graph) <- c('igraph','tbl_graph')
  # graph <- as.igraph(graph)
  
  # Get the ggnetwork data
  network_data <- ggnetwork(graph)
  
  # Remove duplicate rows for labeling
  label_data <- network_data %>% filter(name %in% label_genes) %>% distinct(name, x, y)
  
  # Start plotting
  p <- ggplot(network_data, aes(x = x, y = y))
  
  # Add edges only if xend and yend are present
  if("xend" %in% colnames(network_data) && "yend" %in% colnames(network_data)) {
    p <- p + geom_edges(aes(xend = xend, yend = yend), color = '#3282b8')
  }
  
  # Add nodes, labels, and other styling
  p <- p +
    geom_nodes(aes(color = I(color)), size = 5) +
    geom_label_repel(
      data = label_data, 
      aes(label = name), 
      size = 3,
      box.padding = 0.2,
      point.padding = 0.3,
      segment.color = 'grey50'
    ) +
    ggtitle(' ') +
    theme_blank() + 
    theme_transparent + 
    theme(title = element_text(hjust = 1, face = 'bold'), legend.position="none")
  
  return(p)
}

# compilation of genes selected by all the methods across the different datasets.
heinz_genes <- c()
hotnet2_genes <- c()
sigmod_genes <- c()
dmgwas_genes <- c()
lean_genes <- c()

nsplits <- 5 # Number of splits applied to the data.

for (j in 1:nsplits) {
  # Read genes selected by Heinz
  file_name <- paste0(results, "data_split_", j, "/heinz/selected_genes_split_", j, ".heinz.txt")
  heinz_genes <- unique(c(heinz_genes, read_tsv(file_name, col_types = 'c')$gene))
  
  # Read genes selected by HotNet2
  file_name <- paste0(results, "data_split_", j, "/hotnet2/selected_genes_split_", j, ".hotnet2.tsv")
  hotnet2_genes <- unique(c(hotnet2_genes, read_tsv(file_name, col_types = 'cc')$gene))
  
  # Read genes selected by SigMod
  file_name <- paste0(results, "data_split_", j, "/sigmod/selected_genes_split_", j, ".sigmod.txt")
  sigmod_genes <- unique(c(sigmod_genes, read_tsv(file_name, col_types = 'c')$gene))
  
  # Read genes selected by DMGWAS
  file_name <- paste0(results_dsl2, "data_split_", j, "/dmgwas/selected_genes_split_", j, ".dmgwas.txt")
  dmgwas_genes <- unique(c(dmgwas_genes, read_tsv(file_name, col_types = 'c')$gene))
  
  # Read genes selected by Lean
  file_name <- paste0(results_dsl2, "data_split_", j, "/lean/selected_genes_split_", j, ".lean.txt")
  lean_genes <- unique(c(lean_genes, read_tsv(file_name, col_types = 'c')$gene))
}


# powersets

gene_lists_subsets <- c()
gene_lists_src <- list()
nsplits <- 5 # Number of splits applied to the data.
nmethods <- 5 # Number of biological network methods used

for (j in 1:nsplits) {
  file_name <- paste0(results,"data_split_", j, "/heinz/selected_genes_split_", j, ".heinz.txt")
  gene_lists_src[[length(gene_lists_src) + 1]] <- read_tsv(file_name, col_types = 'c')$gene

  file_name <- paste0(results,"data_split_", j, "/hotnet2/selected_genes_split_", j, ".hotnet2.tsv")
  gene_lists_src[[length(gene_lists_src) + 1]] <- read_tsv(file_name, col_types = 'cc')$gene

  file_name <- paste0(results,"data_split_", j, "/sigmod/selected_genes_split_", j, ".sigmod.txt")
  gene_lists_src[[length(gene_lists_src) + 1]] <- read_tsv(file_name, col_types = 'c')$gene
  
  file_name <- paste0(results_dsl2,"data_split_", j, "/dmgwas/selected_genes_split_", j, ".dmgwas.txt")
  gene_lists_src[[length(gene_lists_src) + 1]] <- read_tsv(file_name, col_types = 'c')$gene

  file_name <- paste0(results_dsl2,"data_split_", j, "/lean/selected_genes_split_", j, ".lean.txt")
  gene_lists_src[[length(gene_lists_src) + 1]] <- read_tsv(file_name, col_types = 'c')$gene
}

for (gene_list in gene_lists_src) {
  for (gene in gene_list) {
    if (gene %in% names(gene_lists_subsets)) {
      gene_lists_subsets[gene] <- gene_lists_subsets[gene] + 1
    } else {
      gene_lists_subsets[gene] <- 1
    }
  }
}

# Calculate the number of genes appearing in at least 1, 2, 3, ..., nmethods * nsplits lists
max_occurrences <- nmethods * nsplits
occurrence_counts <- sapply(1:max_occurrences, function(x) sum(gene_lists_subsets >= x))

# Create a data frame for plotting
occurrence_df <- data.frame(Occurrences = 1:max_occurrences, Frequency = occurrence_counts)

# Plot the data
ggplot(occurrence_df, aes(x = Occurrences, y = Frequency)) +
  geom_bar(stat = "identity", fill = "#016784") +
  geom_text(aes(label = Frequency), vjust = -0.5) +
  geom_line(stat = "identity") +
  labs(
    title = "Number of genes appearing in at least n solutions", 
    x = "Solutions", 
    y = "Number of genes") +
  theme_minimal()

# Initialize an empty list to store the plots
gene_plots <- vector("list", max_occurrences)

# Loop through the process n times
for (i in seq(1, max_occurrences)) {
  # Genes in at least i solutions
  gene_set <- names(gene_lists_subsets[gene_lists_subsets >= i])
  # if (i == 8)
  #   gene_plots[[i]] <- plot_network(ppi, gene_set, node_colors_gene_set11, node_names_gene_set8)
  # else if (i == 5)
  #   gene_plots[[i]] <- plot_network(ppi, gene_set, node_colors_gene_set11, node_names_gene_set5)
  # else
  gene_plots[[i]] <- plot_network(ppi, gene_set)
}

# Create the topology panel
topology_panel <- plot_grid(plotlist = gene_plots, 
                            labels = paste0("Set", 1:max_occurrences),
                            label_size = 28)

# Create the final plot
f1 <- plot_grid(topology_panel, 
                align = "v",        # Aligns vertically
                axis = "tblr")         # Aligns the top and bottom axes
f1

ggsave(paste0(results_dsl2,'../figures/stable_consensus_1+_solutions_5methods.pdf'), f1, 
       width=50, height=50, limitsize = FALSE, bg = "transparent")

# Genes in at least 7 solutions

gene_set7 <- names(gene_lists_subsets[gene_lists_subsets >= 7])
gene_set7_plt <- plot_network(ppi, gene_set7)

# Genes in at least 11 solutions

gene_set11 <- names(gene_lists_subsets[gene_lists_subsets >= 11])
gene_set11_plt <- plot_network(ppi, gene_set11)

gene_set12 <- names(gene_lists_subsets[gene_lists_subsets >= 12])
gene_set15 <- names(gene_lists_subsets[gene_lists_subsets >= 15])
gene_set15_plt <- plot_network(ppi, gene_set11)

# Genes in at least 13 solutions

gene_set13 <- names(gene_lists_subsets[gene_lists_subsets >= 13])
gene_set13_plt <- plot_network(ppi, gene_set13)

# Genes in at least 16 solutions

gene_set16 <- names(gene_lists_subsets[gene_lists_subsets >= 16])
gene_set16_plt <- plot_network(ppi, gene_set16)

# Topology panel

topology_panel <- plot_grid(gene_set16_plt,
                        labels = c("Set16"),
                        label_size = 36)

#ggsave(paste0(results,'../figures/gene_selection8_snp-array.pdf'), topology_panel, 
#       width=25, height=25, limitsize = FALSE, bg = "transparent")

f1 <- plot_grid(topology_panel,
                nrow = 2,
                ncol = 3,
                labels = 'C')
f1

### Consensus of gene sets ###

# Gene set 7

gset7 <- activate(ppi, nodes) %>%
  filter(name %in% gene_set7) %>%
  mutate(gset7 = TRUE)
class(gset7) <- c('igraph','tbl_graph')

# Gene set 11

gset11 <- activate(ppi, nodes) %>%
  filter(name %in% gene_set11) %>%
  mutate(gset11 = TRUE)
class(gset11) <- c('igraph','tbl_graph')

# Gene set 12

gset12 <- activate(ppi, nodes) %>%
  filter(name %in% gene_set12) %>%
  mutate(gset12 = TRUE)
class(gset12) <- c('igraph','tbl_graph')

# Gene set 13

gset13 <- activate(ppi, nodes) %>%
  filter(name %in% gene_set13) %>%
  mutate(gset13 = TRUE)
class(gset13) <- c('igraph','tbl_graph')

# Gene set 14

gset14 <- activate(ppi, nodes) %>%
  filter(name %in% gene_set14) %>%
  mutate(gset14 = TRUE)
class(gset14) <- c('igraph','tbl_graph')

# Gene set 15

gset15 <- activate(ppi, nodes) %>%
  filter(name %in% gene_set15) %>%
  mutate(gset15 = TRUE)
class(gset15) <- c('igraph','tbl_graph')

# Gene set 16

gset16 <- activate(ppi, nodes) %>%
  filter(name %in% gene_set16) %>%
  mutate(gset16 = TRUE)
class(gset16) <- c('igraph','tbl_graph')

# Consensus

options(repr.plot.width=25, repr.plot.height=15)

consensus_gg <- function(gset) {
  
  gset <- gset %>% to_undirected

  as_tibble(gset) %>%
  rename(gene = name) %>%
  write_tsv(paste0(results_dsl2,'../stable_consensus.tsv'))
  
  gset %>%
  ggnetwork %>%
  mutate(name = as.character(name))
}

consensus <- consensus_gg(gset7)

nodes <- mutate(consensus,
                Heinz = as.numeric(name %in% heinz_genes),
                HotNet2 = as.numeric(name %in% hotnet2_genes),
                SigMod = as.numeric(name %in% sigmod_genes),
                dmGWAS = as.numeric(name %in% dmgwas_genes),
                LEAN = as.numeric(name %in% lean_genes),
                Julio = as.numeric(name %in% julio_psoriasis_genes$gene),
                Ran2019 = as.numeric(name %in% ran_psoriasis_genes$gene),
                Tsoi2017 =  as.numeric(name %in% tsoi_psoriasis_genes$`Closest Gene`),
                IL17 = as.numeric(name %in% il17_signaling_pathway$gene),
                IL23 = as.numeric(name %in% il23_signaling_pathway$gene),
                TNF =  as.numeric(name %in% tnf_signaling_pathway$gene)) %>%
  filter(xend == x & yend == y) %>% 
  select(x, y, name, all_of(gen_list), all_of(pathways), all_of(methods)) %>%
  unique

edges <- filter(consensus, xend != x | yend != y)


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
common_all_gen_lists <- common_julio_psoriasis_genes | common_ran_psoriasis_genes | common_tsoi_psoriasis_genes

common_il17_pathway <- nodes$name %in% il17_signaling_pathway$gene
common_il23_pathway <- nodes$name %in% il23_signaling_pathway$gene
common_tnf_pathway <- nodes$name %in% tnf_signaling_pathway$gene
common_all_pathways <- common_il17_pathway | common_il23_pathway | common_tnf_pathway

# CONSENSUS PIES

consensus_pie_gene_lists <- ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_point(data = filter(nodes, !common_all_gen_lists), 
             aes(x = x, y = y), 
             size = 1, color = '#393e46') +
  geom_scatterpie(data = nodes, aes(x = x, y = y, r=0.009), cols = gen_list) +
  geom_label_repel(data = filter(nodes, common_all_gen_lists), 
                   aes(x = x, y = y, label = name), 
                   size = 3, point.padding = 1, max.time = 5, max.iter = 20000, 
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
  geom_scatterpie(data = nodes, aes(x = x, y = y, r=0.009), cols = pathways) +
  geom_label_repel(data = filter(nodes, common_all_pathways), 
                   aes(x = x, y = y, label = name), 
                   size = 3, point.padding = 1, max.time = 5, max.iter = 20000, 
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

consensus_pie_methods <- ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_scatterpie(data = nodes, aes(x = x, y = y, r=0.009), cols = methods) +
  geom_label_repel(data = nodes, 
                   aes(x = x, y = y, label = name), 
                   size = 2, point.padding = 1, max.time = 5, max.iter = 20000, 
                   max.overlaps = 50) +
  coord_fixed() +
  theme_blank() +
  labs(fill = 'Methods') +
  scale_fill_manual(values = method_palette) +
  guides(color = "none") +
  theme(legend.position = 'left',
        legend.text = element_text(size = 20, vjust = .2),
        legend.title = element_text(size = 22, vjust = 1)) +
  theme_transparent

consensus_pie_methods

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
ggsave(paste0(results,'../figures/stable_consensus_pie_gene_lists_wes_ncbi38_gset11.pdf'), consensus_pie_gene_lists, width=20, height=18, bg = "transparent")
ggsave(paste0(results_dsl2,'../figures/stable_consensus_gene_lists_gset7.pdf'), consensus_pie_gene_lists, width=10, height=10, bg = "transparent")
ggsave(paste0(results_dsl2,'../figures/stable_consensus_pathways_gset7.pdf'), consensus_pie_pathways, width=10, height=10, bg = "transparent")
ggsave(paste0(results,'../figures/stable_consensus_pie_pathways_wes_ncbi38_gset11.pdf'), consensus_pie_pathways, width=20, height=18, bg = "transparent")
ggsave(paste0(results_dsl2,'../figures/stable_consensus_methods_gset7.pdf'), consensus_pie_methods, width=15, height=15, bg = "transparent")

consensus_names_gene_lists <- list(consensus_names_julio, consensus_names_ran, consensus_names_tsoi, consensus_pie_gene_lists)

pdf(paste0(results,'../figures/stable_consensus_names_gene_lists.pdf'), width=15, height=15, bg = "transparent")
for (plot in consensus_names_gene_lists) {
  print(plot)
}
dev.off()

consensus_names_pathways <- list(consensus_names_il23, consensus_names_il17, consensus_names_tnf, consensus_pie_pathways)
pdf(paste0(results,'../figures/gene_sets_15solutions/stable_consensus_names_pathways_in14solutions.pdf'), width=15, height=15, bg = "transparent")
for (plot in consensus_names_pathways) {
  print(plot)
}
dev.off()
