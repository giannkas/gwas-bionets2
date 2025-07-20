from bioinfokit import analys, visuz
import pandas as pd
# load dataset as pandas dataframe
#df = analys.get_data('mhat').data
df = pd.read_csv('results/psoriasis/magma/magma_output_50kb_window_snp-array/pso.scores.genes.out_converted.tsv', sep="\t")
