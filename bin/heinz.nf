#!/usr/bin/env nextflow

process heinz {

  publishDir params.out, overwrite: true, mode: "copy"

  input:
    path MAGMA
    path TAB2
    val FDR
    val SPLIT

  output:
    path "selected_genes${SPLIT}.heinz.txt", emit: genes_heinz
      
  """
  #!/usr/bin/env Rscript

  library(tidyverse)
  library(igraph)
  library(BioNet)

  scores <- read_tsv('${MAGMA}')
  net <- read_tsv("${TAB2}") %>%
      rename(gene1 = `Official Symbol Interactor A`, 
             gene2 = `Official Symbol Interactor B`) %>%
      filter(gene1 %in% scores\$Gene & gene2 %in% scores\$Gene) %>%
      select(gene1, gene2) %>%
      graph_from_data_frame(directed = FALSE)
  scores <- filter(scores, Gene %in% names(V(net)))

  # search subnetworks
  pvals <- scores\$Pvalue
  names(pvals) <- scores\$Gene
  fb <- fitBumModel(pvals, plot = FALSE)
  scores <- scoreNodes(net, fb, fdr = ${FDR})

  if (sum(scores > 0)) {
      selected <- runFastHeinz(net, scores)    
      tibble(gene = names(V(selected))) %>% 
          write_tsv('selected_genes${SPLIT}.heinz.txt')
  } else {
      write_tsv(tibble(gene = character()), 'selected_genes${SPLIT}.heinz.txt')
  }
  """

}

params.out = '.'
params.fdr = 0.1
params.i = 0
params.d_samp = 1
params.tab2 = ''
params.scores = ''


workflow {
  
  def tab2 = file(params.tab2)
  def scores = file(params.scores)
  def split = params.i > 0 && params.d_samp != 0 ? "_split_${params.i}" : params.i > 0 && params.d_samp == 0 ? "_chunk_${params.i}" : ""

  heinz(scores, tab2, params.fdr, split)
}
