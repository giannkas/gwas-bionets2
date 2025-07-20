library(tidyverse)
#readr::read_table("/cluster/CBIO/data1/glemoine/ukb_gwas-tools/data/network/biogrid_ppi.tsv") %>% 
ppi <- read_tsv('../data/biogrid_ppi.tsv')
results <- '../results/psoriasis/magma/magma_output_50kb_window/stable_consensus/data_splits/'
#define number of data frames to split into
n <- 5
#split data frame into n equal-sized data frames
chunks <- split(ppi, factor(sort(rank(row.names(ppi))%%n)))

# Loop through each chunk and save as tsv file
for (i in 1:length(chunks)) {
  # Create filename based on chunk number
  filename <- paste0(results, "biogrid_chunk_", i, ".tsv")
  
  # Write the chunk to a tsv file with row.names=FALSE to avoid saving row names
  write.table(chunks[[i]], file = filename, sep = "\t", row.names = FALSE)
}
