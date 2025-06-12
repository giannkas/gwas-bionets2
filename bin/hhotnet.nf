#!/usr/bin/env nextflow

process data_formatting {

  publishDir params.out, overwrite: true, mode: 'copy'

  input:
    path NETWORK
    path SCORES

  output:
    path 'index_gene.tsv', emit: index_gene
    path 'gene_score.tsv', emit: gene_score
    path 'edge_list.tsv', emit: edge_list

  script:
    template 'io/scrs_net_format_hhotnet.R'

}


process similarity_matrix {

  publishDir params.out, overwrite: true, mode: 'copy'

  input:
    path HHOTNET
    path NETWORK
    val BETA

  output:
    path "similarity_matrix.h5", emit: sim_mtx
    path "beta.txt", emit: beta_file

  script:
    """
      python ${HHOTNET}/src/construct_similarity_matrix.py \
        -b   ${BETA} \
        -i   ${NETWORK} \
        -o   similarity_matrix.h5 \
        -bof beta.txt
    """
}

process permuting_network {

  input:
    tuple val(I), path(NETWORK), path(HHOTNET)

  output:
    path "edge_list_${I}.tsv", emit: edge_list_$I

  script:
    """
      python ${HHOTNET}/src/permute_network.py \
        -i ${NETWORK} \
        -s "${I}" \
        -c \
        -o edge_list_${I}.tsv
    """

}

process find_permutation_bins {

  input:
    path NETWORK
    path SCORES
    path IDXGENE
    path HHOTNET

  output:
    path "score_bins.tsv", emit: score_bins

  script:
    """
      python ${HHOTNET}/src/find_permutation_bins.py \
        -gsf ${SCORES} \
        -igf ${IDXGENE} \
        -elf ${NETWORK} \
        -ms 5 \
        -o score_bins.tsv
    """

}

process permuting_scores {

  input:
    tuple val(PERMS), path(BINS), path(SCORES), path(HHOTNET)

  output:
    path "scores_${PERMS}.tsv", emit: scores

  script:
    """
      python ${HHOTNET}/src/permute_scores.py \
        -i ${SCORES} \
        -bf ${BINS} \
        -s "${PERMS}" \
        -o scores_${PERMS}.tsv
    """

}

process construct_hierarchies {

  input:
    tuple val(PERMS), path(SIMMATRIX), path(IDXGENE), path(SCORES), path(HHOTNET)

  output:
    tuple path("hierarchy_edge_list_${PERMS}.tsv"), path("hierarchy_index_gene_${PERMS}.tsv")

  script:
    """
      python ${HHOTNET}/src/construct_hierarchy.py \
        -smf ${SIMMATRIX} \
        -igf ${IDXGENE} \
        -gsf ${SCORES} \
        -helf hierarchy_edge_list_${PERMS}.tsv \
        -higf hierarchy_index_gene_${PERMS}.tsv
    """

}

process processing_hierarchies {

  publishDir params.out, overwrite: true, mode: 'copy'
  
  input:
    tuple path(HEDGE_LIST_0, stageAs: "?/*"), path(HINDEX_GENE_0, stageAs: "?/*"), path(HEDGE_LISTS), path(HINDEX_GENES), val(NETWORKNAME), val(SCORESNAME), path(HHOTNET), val(SPLIT)

  output:
    path "clusters_${NETWORKNAME}_${SCORESNAME}${SPLIT}.hhotnet.tsv", emit: clusters
    path "sizes_${NETWORKNAME}_${SCORESNAME}${SPLIT}.hhotnet.pdf", emit: sizes

  script:
    """
      python ${HHOTNET}/src/process_hierarchies.py \
        -oelf ${HEDGE_LIST_0} \
        -oigf ${HINDEX_GENE_0} \
        -pelf ${HEDGE_LISTS} \
        -pigf ${HINDEX_GENES} \
        -lsb 10 \
        -cf clusters_${NETWORKNAME}_${SCORESNAME}${SPLIT}.hhotnet.tsv \
        -pl ${NETWORKNAME} gene_score \
        -pf sizes_${NETWORKNAME}_${SCORESNAME}${SPLIT}.hhotnet.pdf
    """

}

process performing_consensus {

  publishDir params.out, overwrite: true, mode: 'copy'

  input:
    path CLUSTERS
    path IDXGENE
    path NETWORK
    val NETWORKNAME
    val SCORESNAME
    path HHOTNET
    val SPLIT

  output:
    path "selected_genes${SPLIT}.hhotnet.tsv"
    path "consensus_edges${SPLIT}.hhotnet.tsv"

  script:
    """
      python ${HHOTNET}/src/perform_consensus.py \
        -cf ${CLUSTERS} \
        -igf ${IDXGENE} \
        -elf ${NETWORK} \
        -n ${NETWORKNAME} \
        -s ${SCORESNAME} \
        -t 1 \
        -cnf selected_genes${SPLIT}.hhotnet.tsv \
        -cef consensus_edges${SPLIT}.hhotnet.tsv
    """
}


// Parameters
params.out = "."
params.hhotnet_path = null
params.i = 0
params.d_samp = 1
params.connectivity = 1
params.beta = 0
params.permutations = 5

workflow {
  // input files
  def network = file(params.network)
  def scores = file(params.scores)

  def split = (params.i > 0 && params.d_samp != 0) ? "_split_${params.i}" : 
            (params.i > 0 && params.d_samp == 0) ? "_chunk_${params.i}" : ""

  def img_hhotnet = "/gwas-bionets2/hhotnet"
  def HHOTNET = (params.hhotnet_path != null && params.hhotnet_path != "/default/path") ? params.hhotnet_path : img_hhotnet

  // def network_permutations = 100
  // def heat_permutations = 1000
  def network_permutations = params.permutations
  def beta = params.beta
  def conn = params.connectivity
  def network_name = network.getBaseName()
  def scores_name = scores.getBaseName()

  // processes
  data_formatting(network, scores)
  // folders(network, scores)
  similarity_matrix(HHOTNET, data_formatting.out.edge_list, beta)
  net_nperms = Channel.of(1..4)

  net_perms = data_formatting.out.edge_list
    .combine(net_nperms)
    .map { edge_list, nperm -> [nperm, edge_list, HHOTNET] }

  permuted_nets = permuting_network(net_perms)


  // permutation bins

  find_permutation_bins(
    data_formatting.out.edge_list,
    data_formatting.out.gene_score,
    data_formatting.out.index_gene,
    HHOTNET)

  // permutation scores

  chl1_nperms = Channel.of(1..network_permutations)

  scores_perms = find_permutation_bins.out.score_bins
    .combine(data_formatting.out.gene_score)
    .combine(chl1_nperms)
    .map { score_bins, gene_score, nperm -> [nperm, score_bins, gene_score, HHOTNET] }

  permuting_scores(scores_perms)

  // hierarchies construction

  chl0_nperms = Channel.of(0..network_permutations)

  hierarchies_cns = permuting_scores.out.scores
    .combine(similarity_matrix.out.sim_mtx)
    .combine(data_formatting.out.index_gene)
    .combine(chl0_nperms)
    .map { nscores, sim_mtx, index_gene, nperm -> [nperm, sim_mtx, index_gene, nscores, HHOTNET] }

  cnstd_hierarchies = construct_hierarchies(hierarchies_cns)

  //From here on, it has not been tested

  // hierarchies processing

  initial_hierarchies = cnstd_hierarchies.first()

  head_hedge_list = initial_hierarchies.map { hedge_list, hindex_gene -> hedge_list }
  head_hindex_gene = initial_hierarchies.map { hedge_list, hindex_gene -> hindex_gene }

  tail_hedgel_hindexg = cnstd_hierarchies.filter { hedge_list, hindex_gene ->
    hedge_list != head_hedge_list && hindex_gene != head_hindex_gene
  }

  tail_hedge_list = tail_hedgel_hindexg.map { hedge_list, hindex_gene ->
    hedge_list }

  tail_hindex_gene = tail_hedgel_hindexg.map { hedge_list, hindex_gene ->
    hindex_gene }

  hierarchies_procs = head_hedge_list
    .combine(head_hindex_gene)
    .combine(tail_hedge_list)
    .combine(tail_hindex_gene)
    .map { head_hel, head_hig, tail_hel, tail_hig -> [head_hel, head_hig, tail_hel, tail_hig, network_name, scores_name, HHOTNET, split] }

  processing_hierarchies(hierarchies_procs)

  // performing consensus

  performing_consensus(processing_hierarchies.out.clusters,
    data_formatting.out.index_gene,
    data_formatting.out.edge_list,
    network_name,
    scores_name,
    HHOTNET,
    split)

}