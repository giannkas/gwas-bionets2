#!/usr/bin/env nextflow

process lean {

  publishDir params.out, overwrite: true, mode: "copy"

    input:
      path SCORES
      path EDGELIST
      val REPS
      val SPLIT

    output:
      path "selected_genes${SPLIT}.lean.txt"

    script:

    """
      #!/usr/bin/env Rscript

      library(igraph)
      library(LEANR)
      library(dplyr)
      library(tidyverse)

      # read edgelist
      network <- read_tsv("${EDGELIST}")

      # gene scores
      gene_scores <- read_tsv("${SCORES}")    

      net_genes <- network %>%
        rename(gene1 = `Official Symbol Interactor A`, gene2 = `Official Symbol Interactor B`) %>%
        select(gene1, gene2) %>%
        mutate(gene1 = toupper(gene1), gene2 = toupper(gene2))

      g_scrs <- gene_scores %>%
        select(Gene) %>%
        mutate(Gene = toupper(Gene)) %>%
        distinct(Gene)

      matched_genes <- net_genes %>%
        filter(gene1 %in% g_scrs\$Gene & gene2 %in% g_scrs\$Gene) %>%
        select(gene1, gene2)  %>%
        pivot_longer(cols = c(gene1, gene2), names_to = "gene_type", values_to = "genes") %>%
        distinct(genes) %>%
        arrange(genes)

      gen_scr <- gene_scores %>%
        select(Gene, Pvalue) %>%
        mutate(Gene = toupper(Gene)) %>%
        filter(Gene %in% matched_genes\$genes) %>%
        distinct(Gene, .keep_all = TRUE)

      net <- net_genes %>%
        filter(gene1 %in% g_scrs\$Gene & gene2 %in% g_scrs\$Gene) %>%
        graph_from_data_frame(directed = FALSE)

      scores <- gen_scr[['Pvalue']]
      names(scores) <- gen_scr[['Gene']]

      results <- run.lean(scores, net, n_reps = ${REPS}, add.scored.genes = TRUE)

      significant_genes <- tibble(Gene = rownames(results\$restab[results\$restab[,'PLEAN']<=0.05,]))
      
      if (nrow(significant_genes) == 0) {
        warning("No significant genes found. Creating empty output file.")
        tibble(gene = character(0)) %>%
          write_tsv('selected_genes${SPLIT}.lean.txt')
      } else {
        nh_genes <- significant_genes %>%
          mutate(neighbors = map(Gene, ~ {
            if (.x %in% names(results\$nhs)) {
              results\$nhs[[.x]]
            } else {
              character(0)
            }
          })) %>%
          unnest(neighbors)


        nh_genes_unique <- nh_genes %>% 
          rename(gene = neighbors) %>%
          distinct(gene) %>%
          write_tsv('selected_genes${SPLIT}.lean.txt')

      }
    """

}

params.out = '.'
params.reps = 100000
params.network = ''
params.scores = ''

workflow {

  def net = file(params.network)
  def scores = file(params.scores)
  def reps = params.reps
  def split = params.i > 0 && params.d_samp != 0 ? "_split_${params.i}" : params.i > 0 && params.d_samp == 0 ? "_chunk_${params.i}" : ""

  lean(scores, net, reps, split)
}