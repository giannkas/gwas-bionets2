library(UpSetR)
library(ComplexUpset)
library(ggplot2)
library(tidyverse)
# example of list input (list of named vectors)
#listInput <- list(one = c(1, 2, 3, 5, 7, 8, 11, 12, 13), two = c(1, 2, 4, 5, 
#                                                                 10), three = c(1, 5, 6, 7, 8, 9, 10, 12, 13))
data_upset <- read_tsv("../results/psoriasis/magma/magma_output_50kb_window_snp-array_dsl2/selected_genes_heinz-bhcorrection-hotnet2-sigmod-dmgwas-lean-cgwas.txt")
#data_upset <- read_tsv("../results/psoriasis/magma/magma_output_50kb_window_snp-array/consensus_vs_pathways.txt")
#data_upset <- read_tsv("../results/psoriasis/magma/magma_output_50kb_window_snp-array/consensus_vs_genelists.txt")
#data_upset <- read_tsv("../results/psoriasis/magma/magma_output_50kb_window_wes_ncbi38_tab2gwen/selected_genes_heinz-bhcorrection-hotnet2-sigmod.txt")
# example of expression input
#expressionInput <- c(one = 2, two = 1, three = 2, `one&two` = 1, `one&three` = 4, 
#                     `two&three` = 1, `one&two&three` = 2)

listInput <- as.list(data_upset)
names(listInput) <- c("Heinz","MAGMABH", "HotNet2", "SigMod", "dmGWAS", "LEAN", "cGWAS")
#names(listInput) <- c("Consensus","IL17", "IL23", "TNF")
#names(listInput) <- c("Consensus","Julio", "Ran2019", "Tsoi2017")


ComplexUpset::upset(fromList(listInput), c("Heinz", "HotNet2", "SigMod"), 
                    sort_intersections_by=c('degree', 'cardinality'),
                    set_sizes=FALSE,
                    encode_sets=FALSE,
                    queries=list(upset_query(set=c('Heinz'), color='#73471f'),
                                 upset_query(set=c('HotNet2'), color='#377eb8'),
                                 upset_query(set=c('SigMod'), color='#f4888a')
                                 ),
                    # set_sizes = (upset_set_size()
                    #              + geom_text(aes(label=..count..), hjust=1.1, stat='count')
                    #              + expand_limits(y=800)
                    #              + theme(axis.text.x=element_text(angle=90))
                    #              ),
                    min_degree=1,
                    # queries=list(
                    #           upset_query(
                    #             intersect='BHcorrection',
                    #             color='#016784ff',
                    #             fill='#016784ff',
                    #             only_components=c('intersections_matrix', 'Intersection size')
                    #           )
                    #         ),
                    base_annotations=list(
                      'Size'=(
                        intersection_size(
                          mode='exclusive_intersection',
                          mapping=aes(fill=exclusive_intersection),
                          size=0,
                          text=list(check_overlap=TRUE)
                        )
                        +
                        scale_fill_manual(
                          values = c("Heinz" = "#73471f", "Heinz-HotNet2" = "#233142",
                                     "Heinz-HotNet2-SigMod" = "#233142",
                                     "Heinz-SigMod" = "#233142",
                                     "HotNet2-SigMod" = "#233142",
                                     "HotNet2" = "#377eb8", "SigMod" = "#f4888a")
                        )
                        # + scale_fill_venn_mix(
                        #   data=fromList(listInput),
                        #   guide='none',
                        #   colors=c('Heinz'='red', 'HotNet2'='blue', 'SigMod'='green3')
                        # )
                        # add some space on the top of the bars
                        + scale_y_continuous(expand=expansion(mult=c(0, 0.1)))
                        + theme(
                          # hide grid lines
                          panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),
                          # show axis lines
                          axis.line=element_line(colour='black')
                          )
                        )
                      ),
                    stripes=upset_stripes(
                      geom=geom_segment(size=12),  # make the stripes larger
                      colors=c('grey95', 'white')
                    )
                  )
  
UpSetR::upset(fromList(listInput), nsets = 7, order.by = "freq", text.scale = c(3,3,3,3,3,3)) 
