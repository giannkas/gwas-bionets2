import argparse
import pandas as pd

# def split_chunks_char_limit(text, max_len=22):
# 	words = text.split()
# 	result = []
# 	i = 0
# 	while i < len(words):
# 		triple = ' '.join(words[i:i+3])
# 		double = ' '.join(words[i:i+2])
# 		if len(triple) <= max_len and i + 2 <= len(words):
# 			result.append(triple)
# 			i += 3
# 		elif len(double) <= max_len and i + 1 <= len(words):
# 			result.append(double)
# 			i += 2
# 		else:
# 			result.append(words[i])
# 			i += 1
# 	return f'\\makecell{{{" \\\\ ".join(result)}}}'

def convert_tsv_to_latex(tsv_file = "", qngenes = 1, output_tex = ""):
	# Load the data
	df = pd.read_csv(tsv_file, sep="\t")
	
	# Filter entries where 'Adjusted P-value' is less than 0.05, 
	# and Overlap is greater or equal than 0.05 with at least half of genes in the query being in the overlapping genes.
	df = df[df['Adjusted P-value'] < 0.05]
	df = df[(df['Overlap'].apply(eval) >= 0.05) | (df['Genes'].apply(lambda x: len(x.split(';'))) / qngenes >= 0.5)]
	#df['Term'] = df['Term'].apply(split_chunks_char_limit)


	# Drop the 'Old P-value' and 'Old Adjusted P-value' columns
	df = df.drop(columns=['Old P-value', 'Old Adjusted P-value'])

	# Round numerical columns to 5 decimal places
	df['P-value'] = df['P-value'].apply(lambda x: f"{x:.5e}")
	df['Adjusted P-value'] = df['Adjusted P-value'].apply(lambda x: f"{x:.5e}")
	df['Odds Ratio'] = df['Odds Ratio'].apply(lambda x: f"{x:.5e}")
	
	# Generate LaTeX table header
	latex_table = """\\begin{table}[ht!]
\\begin{tabular}{C{3cm}C{2cm}C{2cm}C{2cm}C{2cm}C{2cm}}
\\hline
\\textbf{Pathway} & \\textbf{Overlap} & \\textbf{P-value} & \\textbf{Adjusted P-value} & \\textbf{Odds Ratio} & \\textbf{Genes} \\\\ \\hline
"""
		
	# Iterate through the rows to populate LaTeX table
	for _, row in df.iterrows():
		pathway = row['Term']
		overlap = row['Overlap']
		p_value = row['P-value']
		adj_p_value = row['Adjusted P-value']
		odds_ratio = row['Odds Ratio']
		genes = row['Genes'].split(';')
		
		first_entry = True
		for gene in genes:
			if len(genes) >= 1 and first_entry:
				latex_table += f"{pathway} & {overlap} & {p_value} & {adj_p_value} & {odds_ratio} & \\makecell[c]{{{gene} \\\\"
				first_entry = False
			else:
				latex_table += f" {gene} \\\\ "
		latex_table = latex_table.rstrip(" \\\\") + "} \\\\ \\hline\n"
	
	# Close the table
	latex_table += "\\end{tabular}\n\\end{table}"
	
	# Save to output file
	with open(output_tex, "w") as f:
		f.write(latex_table)
	
	print(f"Table pathways pvalues generated in {output_tex}")

def main():
	parser = argparse.ArgumentParser(description="Creating pathways pvalue LaTeX table.")
	parser.add_argument('--qngenes', type=int, required=True, help='number of genes in the query')
	args = parser.parse_args()
	convert_tsv_to_latex("../results/psoriasis/enrichment/magma_output_50kb_window_snp-array_dsl2/reactome_pathways_2024_table_net_105genes_consensus.txt", args.qngenes, "pathways_pvalues.tex")

if __name__ == "__main__":
	main()