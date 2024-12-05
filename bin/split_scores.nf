#!/usr/bin/env nextflow

/////////////////////////////////////
//  CHUNK PREPARATION
/////////////////////////////////////
process make_chunks {

  publishDir { params.out + (params.k > 1 ? '/scores_chunks/scores_chunk_' + task.index : '') }, overwrite: true, mode: "copy"

  input:
    path scores_file
    val filename
    val K
    val I from 1..K

  output:
    path { filename + (K > 1 ? "_chunk_" + I : '') + ".tsv" } into chunks

  script:
      
  """
    #!/usr/bin/env Rscript

      library(tidyverse)
      scores <- read_tsv('${SCORES}')

      #define number of data frames to split into
      n <- ${K}
      if (n > 1) {
        chunk_suffix <- paste0("_chunk_", ${I})
      } else {
        chunk_suffix <- ""
      }
      nrows <- nrow(scores)

      #split data frame into n equal-sized data frames
      #chunks <- split(scores, factor(sort(rank(row.names(scores))%%n)))

      #loop through each chunk and save as tsv file
      #for (i in seq_along(chunks)) {

      # Calculate the range of rows to exclude in this chunk
      exclude_start <- ceiling((${I} - 1) * nrows / n) + 1
      exclude_end <- ceiling(${I} * nrows / n)

      # Create the chunk by excluding the specified rows
      chunk <- scores[-(exclude_start:exclude_end), ]

      #create filename based on chunk number
      outfile <- paste0("${filename}", chunk_suffix, ".tsv")
      
      #write.table(chunks[[i]], file = outfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
      #write the chunk to a tsv file with row.names=FALSE to avoid saving row names
      write.table(chunk, file = outfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      #}
  """  

}

// Parameters
params.out = '.'
params.k = 1
params.scores = null

// Input Validation
if (!params.scores) {
    log.error "Error: --scores parameter must be specified."
    exit 1
}

// Extract filename and basename from the scores file
scores_file = file(params.scores)
filename = scores_file.baseName

// Workflow
workflow {
  val K = params.k
  def chunks = make_chunks(scores_file, filename, K)
}
