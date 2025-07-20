library(hypeR)

#gene_list <- readLines("../data/c7_all_v2024-1_symbols_unique.gmt")

# geneset_data <- readLines("../data/c7.all.v2024.1.Hs.symbols_nourls_uppercase.gmt")
pathway_rpgd <- readLines("../results/psoriasis/enrichment/magma_output_50kb_window_snp-array/gene_comparison_mitochondria_ribosome/mitochondrial_ribosomal_genes_rpgd.txt")
pathway_reactome <- readLines("../results/psoriasis/enrichment/magma_output_50kb_window_snp-array/gene_comparison_mitochondria_ribosome/mitochondrial_ribosomal_participating_molecules_[R-HSA-5368287]_reactome.tsv")
pathway_gwas_catalog <- readLines("../results/psoriasis/enrichment/magma_output_50kb_window_snp-array/gene_comparison_mitochondria_ribosome/mitochondrial_ribosomal_genes_gwas_catalog.txt")
genesets <- list()
genesets[["mito_ribo_genes_rpgd"]] <- pathway_rpgd
genesets[["mito_ribo_genes_reactome"]] <- pathway_reactome
genesets[["mito_ribo_genes_gwas_catalog"]] <- pathway_gwas_catalog

for (line in geneset_data) {
  words <- strsplit(line, "\t")[[1]]
  # first word is the gene set name
  geneset_name <- words[1]
  # remaining words are the gene symbols
  gene_symbols <- words[-1]
  # remove any empty gene symbols
  gene_symbols <- gene_symbols[gene_symbols != ""]
  # Assign the gene symbols to the gene set name in the list
  genesets[[geneset_name]] <- gene_symbols
}

g5_snparray_dsl1_11th <- c("MRPL4", "MRPS18B", "MRPL38", "MRPL24", "CS")

g12_snparray_dsl1_5th <- c("MRPL4", "MRPL51", "DDX28", "MRPL52", "XRN2", "MRPS18B",
    "PPIF", "MRPL38", "MRPL24", "MRPL14", "CS", "VARS2")

g8_snparray_dsl1_5th <- c("BPTF", "IRF1", "MAZ", "STAT2", "KEAP1", "IRF9", "HMGXB4", 
    "STK4")

g18_snparray_dsl1_6th <- c("GATA4", "CDC7", "PSMA6", "PSMB5", "DBF4", "PSMB9", "CLSPN",
    "PSMA8", "MSH4", "MTF2", "MSH5", "UBD", "CDK3", "ID2", "CDK2", "PSMG2", 
    "TIMELESS", "JARID2")

signatures <- list(
  G12_snparray_5th = g12_snparray_5th, 
  G8_snparray_5th = g8_snparray_5th, 
  G18_snparray_6th = g18_snparray_6th)

str(signatures)
print(genesets)

population <- readLines("../data/biogrid_tab2_interactorsAB_glemoine_merged_uniq_common_snp-array_magma_genes.tsv")

hyp_obj5 <- hypeR(g5_snparray_dsl1_11th, genesets, test = "hypergeometric", background = population)
hyp_obj12 <- hypeR(g12_snparray_dsl1_5th, genesets, test = "hypergeometric", background = population)
hyp_obj8 <- hypeR(g8_snparray_dsl1_5th, genesets, test = "hypergeometric", background = population)
hyp_obj18 <- hypeR(g18_snparray_dsl1_6th, genesets, test = "hypergeometric", background = population)

# hyp_dots(hyp_obj)
# hyp_show(hyp_obj)
# hyp_emap(hyp_obj)
# hyp_hmap(hyp_obj)

# Generate markdown report
hyp_to_rmd(hyp_obj5,
           show_emaps = FALSE,
           file_path="report_hype_g5_snparray_11th.rmd",
           title="Enrichement analysis for the set of 5 genes in the 11th solution (SNP-array, DSL1)",
           author="Giann Karlo Aguirre-Samboní",
           versioning = TRUE,
)

hyp_to_rmd

hyp_to_rmd(hyp_obj12,
           show_emaps = FALSE,
           file_path="report_hype_g12_snparray_5th.rmd",
           title="Enrichement analysis for the set of 12 genes in the 5th solution (SNP-array)",
           author="Giann Karlo Aguirre-Samboní",
           versioning = TRUE,
)

hyp_to_rmd(hyp_obj8,
           show_emaps = FALSE,
           file_path="report_hype_g8_snparray_5th.rmd",
           title="Enrichement analysis for the set of 8 genes in the 5th solution (SNP-array)",
           author="Giann Karlo Aguirre-Samboní",
           versioning = TRUE,
)

hyp_to_rmd(hyp_obj18,
           show_emaps = FALSE,
           file_path="report_hype_g18_snparray_6th.rmd",
           title="Enrichement analysis for the set of 18 genes in the 6th solution (SNP-array)",
           author="Giann Karlo Aguirre-Samboní",
           versioning = TRUE,
)

