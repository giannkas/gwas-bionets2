#!/usr/bin/env nextflow

process sigmod {

  publishDir params.out, overwrite: true, mode: 'copy'

  input:
    path MAGMA_OUT
    path TAB2
    path sigmod_path
    val lambdamax
    val nmax
    val maxjump
    val split

  output:
    path "selected_genes${split}.sigmod.txt" into genes_sigmod

  script:
    template 'discovery/run_sigmod.R'
}

params.out = '.'
params.lambdamax = 1
params.nmax = 300
params.maxjump = 10
params.sigmod_path = null
params.i = 0
params.d_samp = 1


def img_sigmod = '/gwas-bionets/sigmod'
def sigmod_path = params.sigmod_path ? file(params.sigmod_path) : img_sigmod

// Split suffix
def split = (params.i > 0 && params.d_samp != 0) ? "_split_${params.i}" : 
            (params.i > 0 && params.d_samp == 0) ? "_chunk_${params.i}" : ""
}

workflow {
  // input files
  def magma = file(params.scores)
  def tab2 = file(params.tab2)

  // run sigmod process
  genes_sigmod = sigmod(magma, tab2, sigmod_path, params.lambdamax, params.nmax, params.maxjump, split)
}
