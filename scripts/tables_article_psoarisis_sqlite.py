import pandas as pd
import sqlite3
from sys import stdin
import argparse

def pathways_pvalue(df, db_file, first):
    
  # Drop the 'Old P-value' and 'Old Adjusted P-value' columns
  df = df.drop(columns=['Old P-value', 'Old Adjusted P-value'])
  
  # Round numerical columns to 5 decimal places
  df['P-value'] = df['P-value'].apply(lambda x: f"{x:.5e}")
  df['Adjusted P-value'] = df['Adjusted P-value'].apply(lambda x: f"{x:.5e}")
  df['Odds Ratio'] = df['Odds Ratio'].apply(lambda x: f"{x:.5e}")
  
  # Expand the genes column
  records = []
  for _, row in df.iterrows():
      genes = row['Genes'].split(';')
      for gene in genes:
          records.append([row['Source'], row['Module name'], row['Term'], row['Overlap'], row['P-value'], row['Adjusted P-value'], row['Odds Ratio'], gene])
  
  # Create a new DataFrame
  df_expanded = pd.DataFrame(records, columns=['source', 'module_name', 'pathway', 'overlap', 'p-value', 'adjusted_p-value', 'odds_ratio', 'gene'])

  # Save to SQLite database
  conn = sqlite3.connect(db_file)
  if first:
    df_expanded.to_sql('pathway_pvalue', conn, if_exists='replace', index=False)
  else:
    df_expanded.to_sql('pathway_pvalue', conn, if_exists='append', index=False)
  conn.close()
  
  print(f"SQLite database saved to {db_file}")


def check_gene_in_method(gene, df_methods):
  return ["Yes" if gene in df_methods[col].values else "No" for col in df_methods.columns]

def count_gene_occurrences(gene, gene_lists_src):
  heinz_count = sum(gene in gene_lists_src[i] for i in range(0, 15, 3))
  hotnet2_count = sum(gene in gene_lists_src[i] for i in range(1, 15, 3))
  sigmod_count = sum(gene in gene_lists_src[i] for i in range(2, 15, 3))
    
  return heinz_count, hotnet2_count, sigmod_count

def load_gene_lists(results_path):
  nsplits = 5
  nmethods = 3
  gene_lists_src = []
  
  for j in range(1, nsplits + 1):
    heinz_file = f"{results_path}/data_split_{j}/heinz/selected_genes_split_{j}.heinz.txt"
    hotnet2_file = f"{results_path}/data_split_{j}/hotnet2/selected_genes_split_{j}.hotnet2.tsv"
    sigmod_file = f"{results_path}/data_split_{j}/sigmod/selected_genes_split_{j}.sigmod.txt"
    
    for file in [heinz_file, hotnet2_file, sigmod_file]:
      try:
        genes = pd.read_csv(file, sep="\t", usecols=["gene"]).squeeze().tolist()
        gene_lists_src.append(genes)
      except FileNotFoundError:
        gene_lists_src.append([])
  
  return gene_lists_src

def gene_pathways(df_pathways, df_scores_dict, df_methods, db_file, gene_lists_src, first):
  
  records = []
  for _, row in df_pathways.iterrows():
    source = row['Source']
    module = row['Module name']
    pathway = row['Term']
    genes = row['Genes'].split(';')
    
    for gene in genes:
      bh_value = df_scores_dict.get(gene, '-')
      heinz, magma, hotnet, sigmod = check_gene_in_method(gene, df_methods)
      heinz_5, hotnet2_5, sigmod_5 = count_gene_occurrences(gene, gene_lists_src)

      records.append([source, module, pathway, gene, bh_value, magma, heinz, hotnet, sigmod, heinz_5, hotnet2_5, sigmod_5])
  
  df_expanded = pd.DataFrame(records, columns=["source", "module_name", "pathway", "gene", "bh-value", "magmabh", "heinz", "hotnet2", "sigmod", "heinz_5", "hotnet2_5", "sigmod_5"])
  
  conn = sqlite3.connect(db_file)
  if first:
    df_expanded.to_sql('gene_pathways', conn, if_exists='replace', index=False)
  else:
    df_expanded.to_sql('gene_pathways', conn, if_exists='append', index=False)
  conn.close()
  
  print(f"SQLite database saved to {db_file}")


def main():
  parser = argparse.ArgumentParser(description="Creating or adding more entries to pathways pvalue table in pathways and genes database.")
  parser.add_argument('--pathways', type=str, required=True, help='Input pathways path')
  parser.add_argument('--first', action='store_true', help='Flag to indicate whether this is the very first entry in the database, and REPLACE an old one')
  parser.add_argument('--module', type=str, required=True, help='Input module name or label')
  parser.add_argument('--source', type=str, required=True, help='Input source database name or label')
  parser.add_argument('--alpha', type=float, default=0.05, required=False, help='Input alpha or threshold value to exclude pathways')
  parser.add_argument('--scores', type=str, required=False, help='Input scores path')
  parser.add_argument('--method_sel', type=str, required=False, help='Input method gene selection path')
  parser.add_argument('--results_path', type=str, required=False, help='Path to results folder where all different method\'s outputs are located')
  parser.add_argument('--database', type=str, required=False, help='If a database already exists, then do not create a new one.')
  parser.add_argument('--output', type=str, default="pathways_and_genes", required=False, help='Output database path')
  parser.add_argument('--verbose', action='store_true', help='Increase output verbosity')

  args = parser.parse_args()

  if args.verbose: print(f"Reading input file: {args.pathways}")

  df_pathways = pd.read_csv(args.pathways, sep="\t")

  if args.verbose: print(f"Excluding pathways that do not meet the established threshold: {args.alpha}")
  
  df_pathways = df_pathways[df_pathways['Adjusted P-value'] < args.alpha]

  if args.verbose: print(f"Adding {args.source} and {args.module} columns to the table")

  df_pathways.insert(0, 'Source', args.source)
  df_pathways.insert(1, 'Module name', args.module)

  if args.scores is None and args.method_sel is None:
    if args.verbose: print(f"Creating or adding entries to patway_pvalue table")
    pathways_pvalue(df_pathways, args.output if args.database is None else args.database, args.first)

  elif args.scores is not None and args.method_sel is not None and args.results_path is not None:
    if args.verbose: print(f"Reading input file: {args.scores}")
    df_scores = pd.read_csv(args.scores, sep="\t")
    df_scores_dict = df_scores.set_index(df_scores.columns[2])[df_scores.columns[-1]].to_dict()

    if args.verbose: print(f"Reading input file: {args.method_sel}")
    df_methods = pd.read_csv(args.method_sel, sep="\t", header=None)

    if args.verbose: print(f"Processing files in: {args.results_path}")
    gene_lists_src = load_gene_lists(args.results_path)
    gene_pathways(df_pathways, df_scores_dict, df_methods, args.output if args.database is None else args.database, gene_lists_src, args.first)



# Example usage
if __name__ == '__main__':
  main()