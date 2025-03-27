#!/usr/bin/env nextflow

process make_network {

  input:
    path tab2

  output:
    path 'node_index.tsv' into node_index
    path 'edge_list.tsv' into edge_list

  script:
    template 'io/tab2_2hotnet.R'

}

process sparse_scores {

  publishDir params.out, overwrite: true, mode: 'copy'

  input:
    path magma
    val lfdr_cutoff
    val split

  output:
    path "scored_genes${split}.sparse.txt" into sparse_scores
    path "lfdr_plot${split}.pdf"

  """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(twilight)
    library(cowplot)

    theme_set(theme_cowplot())

    scores <- read_tsv('${SCORES}')

    lfdr <- twilight(scores\$Pvalue, B=1000)
    lfdr <- tibble(Gene = scores\$Gene[as.numeric(rownames(lfdr\$result))],
                magma_p = scores\$Pvalue[as.numeric(rownames(lfdr\$result))],
                lfdr = lfdr\$result\$fdr)

    ggplot(lfdr, aes(x = magma_p, y = 1 - lfdr)) +
     geom_line() +
     geom_vline(xintercept = ${CUTOFF}, color = 'red') +
     labs(x = 'P-value', y = '1 - lFDR')
    ggsave('lfdr_plot${SPLIT}.pdf', width=7, height=6)

    lfdr %>%
     mutate(Pvalue = ifelse(magma_p < ${CUTOFF}, magma_p, 1)) %>%
     write_tsv('scored_genes${SPLIT}.sparse.txt')
  """

}

process magma2hotnet {

  input:
    path sparse_scores

  output:
    path 'scores.ht' into scores

  script:
    template 'io/magma2hotnet.R'

}


process make_h5_network {

  input:
    path hotnet2
    path node_index
    path edge_list
    val beta

  output:
    path "ppin_ppr_${beta}.h5" into h5
    path 'permuted' into permutations

  script:

    """
      python2 ${HOTNET2}/makeNetworkFiles.py \
        --edgelist_file ${EDGE_LIST} \
        --gene_index_file ${NODE_IDX} \
        --network_name ppin \
        --prefix ppin \
        --beta ${BETA} \
        --cores -1 \
        --num_permutations ${network_permutations} \
        --output_dir .
    """

 }

process make_heat_data {

  input:
    path hotnet2
    path scores

  output:
    path 'heat.json' into heat

  script:
    """
    python2 ${HOTNET2}/makeHeatFile.py \
      scores \
      --heat_file ${SCORES} \
      --output_file heat.json \
      --name gwas
    """
}

process hotnet2 {

  input:
    path hotnet2
    path heat
    path network
    path permutations
    val beta

  output:
    path 'consensus/subnetworks.tsv' into subnetworks

  script:
    """
    python2 ${HOTNET2}/HotNet2.py \
      --network_files ${NETWORK} \
      --permuted_network_path ${PERMS}/ppin_ppr_${BETA}_##NUM##.h5 \
      --heat_files ${HEAT} \
      --network_permutations ${network_permutations} \
      --heat_permutations ${heat_permutations} \
      --num_cores -1 \
      --output_directory .
    """

}

process process_output {

  publishDir params.out, overwrite: true, mode: 'copy'

  input:
    path subnetworks
    val split

  output:
    path "selected_genes${split}.hotnet2.tsv" into genes_hotnet2

  script:
    """
    #!/usr/bin/env Rscript

     library(tidyverse)

     read_tsv('${SUBNETWORKS}', col_types = 'cc', comment = '#', col_names = F) %>%
         select(X1) %>%
         mutate(cluster = 1:n()) %>%
         separate_rows(X1, sep = ' ') %>%
         rename(gene = X1) %>%
         write_tsv('selected_genes${SPLIT}.hotnet2.tsv')
    """

}

// Parameters
params.out = "."
params.hotnet2_path = null
params.i = 0
params.d_samp = 1
params.lfdr_cutoff = 0.05

def split = (params.i > 0 && params.d_samp != 0) ? "_split_${params.i}" : 
            (params.i > 0 && params.d_samp == 0) ? "_chunk_${params.i}" : ""

def img_hotnet2 = '/gwas-bionets/hotnet2'
def HOTNET2 = params.hotnet2_path ? file(params.hotnet2_path) : img_hotnet2

def network_permutations = 100
def heat_permutations = 1000
def beta = 0.4

workflow {
  // input files
  def tab2 = file(params.tab2)
  def magma = file(params.scores)

  // processes
  node_index, edge_list = make_network(tab2)
  sparse_scores = sparse_scores(magma, params.lfdr_cutoff, split)
  scores = magma2hotnet(sparse_scores)

  h5, permutations = make_h5_network(HOTNET2, node_index, edge_list, beta)
  heat = make_heat_data(HOTNET2, scores)

  subnetworks = hotnet2(HOTNET2, heat, h5, permutations, beta)
  genes_hotnet2 = process_output(subnetworks, split)
}