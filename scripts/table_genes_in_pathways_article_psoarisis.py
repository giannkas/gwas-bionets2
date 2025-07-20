import argparse
import pandas as pd

# script to format Latex tables for the Psoriasis article

def split_chunks_char_limit(text, max_len=22):
	words = text.split()
	result = []
	i = 0
	while i < len(words):
		triple = ' '.join(words[i:i+3])
		double = ' '.join(words[i:i+2])
		if len(triple) <= max_len and i + 2 <= len(words):
			result.append(triple)
			i += 3
		elif len(double) <= max_len and i + 1 <= len(words):
			result.append(double)
			i += 2
		else:
			result.append(words[i])
			i += 1
	return f'\\makecell{{{" \\\\ ".join(result)}}}'

def load_data(qngenes = 1):
	# Load (A): KEGG pathways
	df_a = pd.read_csv("../results/psoriasis/enrichment/magma_output_50kb_window_snp-array_dsl2/reactome_pathways_2024_table_net_105genes_consensus.txt", sep="\t")
	df_a = df_a[df_a['Adjusted P-value'] < 0.05]  # Filter by significance
	# Filter entries where 'Adjusted P-value' is less than 0.05, 
	# and Overlap is greater or equal than 0.05 with at least half of genes in the query being in the overlapping genes.
	df_a = df_a[(df_a['Overlap'].apply(eval) >= 0.05) | (df_a['Genes'].apply(lambda x: len(x.split(';'))) / qngenes >= 0.5)]
	df_a['Term'] = df_a['Term'].apply(split_chunks_char_limit)
	
	# Load (C): Gene BH values
	df_c = pd.read_csv("../results/psoriasis/magma/magma_output_50kb_window_snp-array/pso.scores.genes.out_converted_bhcorrection_all.tsv", sep="\t")
	df_c_dict = df_c.set_index(df_c.columns[2])[df_c.columns[-1]].to_dict()
	
	# Load (B): Gene categorization
	df_b = pd.read_csv("../results/psoriasis/magma/magma_output_50kb_window_snp-array_dsl2/selected_genes_heinz-bhcorrection-hotnet2-sigmod-dmgwas-lean-cgwas.txt", sep="\t", header=None)
	
	return df_a, df_c_dict, df_b

def check_gene_in_b(gene, df_b):
	return ["Yes" if gene in df_b[col].values else "No" for col in df_b.columns]

def generate_latex_table(df_a, df_c_dict, df_b):
	latex_table = """\\begin{table*}[ht!]
\\noindent\\makebox[\\textwidth]{%
\\resizebox{1.5\\textwidth}{!}{%
\\begin{tabular}{C{4cm}C{2cm}C{2cm}C{1.5cm}C{2cm}C{1.5cm}C{1.5cm}C{1.8cm}C{1.5cm}C{1.5cm}}
\\hline
\\textbf{Pathway} & \\textbf{Genes} & \\textbf{BH-value} & \\textbf{Heinz} & \\textbf{MagmaBH} & \\textbf{HotNet2} & \\textbf{SigMod} & \\textbf{dmGWAS} & \\textbf{LEAN} & \\textbf{cGWAS} \\\\ \\hline
"""
	
	for _, row in df_a.iterrows():
		pathway = row['Term']
		genes = row['Genes'].split(';')
		first_entry = True
		
		for gene in genes:
			bh_value = df_c_dict.get(gene, '-')
			heinz, magma, hotnet, sigmod, dmgwas, lean, cgwas = check_gene_in_b(gene, df_b)
			
			if first_entry:
				latex_table += f"\\multirow{{{len(genes)}}}{{*}}{{{pathway}}} " if len(genes) > 1 else f"{pathway} "
				first_entry = False
			else:
				latex_table += " "
			
			latex_table += f"& {gene} & {bh_value} & {heinz} & {magma} & {hotnet} & {sigmod} & {dmgwas} & {lean} & {cgwas}\\\\ \\cmidrule(lr){{2-10}}\n"
		latex_table = latex_table.rstrip(" \\cmidrule(lr){2-10}\n") + " \\\\ \\hline\n"
	
	latex_table += "\\end{tabular}\n}\n}\n\\end{table*}"
	
	return latex_table

def main():
	parser = argparse.ArgumentParser(description="Creating genes in pathways LaTeX table.")
	parser.add_argument('--qngenes', type=int, required=True, help='number of genes in the query')
	args = parser.parse_args()


	df_a, df_c_dict, df_b = load_data(args.qngenes)
	latex_code = generate_latex_table(df_a, df_c_dict, df_b)
	
	with open("genes_pathways.tex", "w") as f:
		f.write(latex_code)
	
	print("Table genes in pathways generated in genes_pathways.tex")

if __name__ == "__main__":
	main()
