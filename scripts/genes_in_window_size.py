# python3 genes_in_window_size.py --geneloc additional_tools/magma/genome_data/NCBI37.3/NCBI37.3.gene.loc 
# --geneset ../results/psoriasis/enrichment/magma_output_50kb_window_snp-array/module_5genes_11thsolution.txt --window 10000

import pandas as pd
import argparse

def find_genes(query_gene, df_geneloc, window):
  same_chr = df_geneloc[df_geneloc['chrom'] == query_gene['chrom']]

  # remove the same query genes
  same_chr = same_chr[same_chr['gene_name'] != query_gene['gene_name']]

  if query_gene['strand'] == '+':
    upstream = same_chr[same_chr['end'] <= query_gene['start']]
    downstream = same_chr[same_chr['start'] >= query_gene['end']]

    upstream = upstream.assign(dist=query_gene['start'] - upstream['end'])
    downstream = downstream.assign(dist=downstream['start'] - query_gene['end'])

  else:  # strand is -
    upstream = same_chr[same_chr['start'] >= query_gene['end']]
    downstream = same_chr[same_chr['end'] <= query_gene['start']]

    upstream = upstream.assign(dist=upstream['start'] - query_gene['end'])
    downstream = downstream.assign(dist=query_gene['start'] - downstream['end'])
  
  # # apply distance filter (window)
  upstream = upstream[upstream['dist'] <= window]
  downstream = downstream[downstream['dist'] <= window]

  # # Sort and take closest if any
  upstream = upstream.sort_values('dist')
  downstream = downstream.sort_values('dist')

  upstream = upstream.assign(gene_ref=query_gene['gene_name'])
  downstream = downstream.assign(gene_ref=query_gene['gene_name'])

  return upstream, downstream

def parse_args():
  parser = argparse.ArgumentParser(description="Finding genes within a given window size of base pairs for each gene in a geneset.")
  parser.add_argument('--geneset', type=str, required=True, help='Input geneset path')
  parser.add_argument('--first', action='store_true', help='Flag to indicate whether this is the very first entry in the database, and REPLACE an old one')
  parser.add_argument('--geneloc', type=str, required=True, help='Input gene position reference path')
  parser.add_argument('--window', type=float, default=10000, required=False, help='Input window size in base pairs')
  parser.add_argument('--verbose', action='store_true', help='Increase output verbosity')
  parser.add_argument('--out', type=str, required=False, help='Output file path')

  args = parser.parse_args()

  if args.verbose: print(f"Reading geneset file: {args.geneset}")

  df_geneset = pd.read_csv(args.geneset, sep="\t", names=['gene_name'])

  if args.verbose: print(f"Reading geneloc file: {args.geneloc}")

  if args.out is None:
    args.out = f"up_down_genes_{int(args.window)}bp.tsv"

  columns = ['gene_id', 'chrom', 'start', 'end', 'strand', 'gene_name']
  df_geneloc = pd.read_csv(args.geneloc, sep="\t", names=columns)

  # print(df_geneset)
  # print(df_geneloc)
  # print(args.window)

  query_genes = df_geneloc[df_geneloc['gene_name'].isin(df_geneset['gene_name'])]

  # print(str(query_genes['chrom']) == str(df_geneloc['chrom']))
  # print(df_geneloc['chrom'])

  #find_genes(query_genes[query_genes['gene_name'] == "MRPL4"], df_geneloc, args.window)

  results = pd.DataFrame()

  for _, query_gene in query_genes.iterrows():
    upstream, downstream = find_genes(query_gene, df_geneloc, args.window)

    streams = pd.concat([upstream, downstream], ignore_index=True)

    results = pd.concat([results, streams], ignore_index=True)
    

  results.to_csv(args.out, sep='\t', index=False)


if __name__ == "__main__":
  pd.set_option('display.max_rows', 350)
  args = parse_args()

