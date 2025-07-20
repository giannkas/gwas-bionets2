library(ggplot2)

# Read data
data_mr <- read.table("../results/psoriasis/enrichment/magma_output_50kb_window_snp-array/gene_comparison_mitochondria_ribosome/mitochondrial_ribosomal_common_snp_array_magma_genes_pvalue_sorted.txt", header=FALSE, stringsAsFactors=FALSE)
data_m <- read.table("../results/psoriasis/enrichment/magma_output_50kb_window_snp-array/gene_comparison_mitochondria_ribosome/mitochondrial_common_snp_array_magma_genes_pvalue_sorted_excl_ribosomal.txt", header=FALSE, stringsAsFactors=FALSE)
data_r <- read.table("../results/psoriasis/enrichment/magma_output_50kb_window_snp-array/gene_comparison_mitochondria_ribosome/ribosomal_common_snp_array_magma_genes_pvalue_sorted.txt", header=FALSE, stringsAsFactors=FALSE)

colnames(data_mr) <- c("Gene", "P_value")
colnames(data_m) <- c("Gene", "P_value")
colnames(data_r) <- c("Gene", "P_value")

# Plot histogram
qqnorm(data_mr$P_value, main="Q-Q Plot of P-values - Mitochondrial-ribosomal genes", pch=19, col="blue")
qqline(data_mr$P_value, col="red", lwd=2)  # Adds reference line

qqnorm(data_m$P_value, main="Q-Q Plot of P-values - Mitochondrial genes", pch=19, col="blue")
qqline(data_m$P_value, col="red", lwd=2)  # Adds reference line
# remove those genes in mitochondrial-ribosomal to only examine mitochondrial genes

qqnorm(data_r$P_value, main="Q-Q Plot of P-values - Ribosomal genes", pch=19, col="blue")
qqline(data_r$P_value, col="red", lwd=2)  # Adds reference line

# Histograms 

ggplot(data_mr, aes(x = P_value)) +
  geom_histogram(binwidth = 0.02, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Ribosomal Gene P-values",
       x = "P-value", y = "Frequency") +
  theme_minimal()

ggplot(data_m, aes(x = P_value)) +
  geom_histogram(binwidth = 0.02, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Mitochondrial Gene P-values",
       x = "P-value", y = "Frequency") +
  theme_minimal()

ggplot(data_r, aes(x = P_value)) +
  geom_histogram(binwidth = 0.02, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Ribosomal Gene P-values",
       x = "P-value", y = "Frequency") +
  theme_minimal()

# Distribution of P-values by gene

ggplot(data_mr, aes(x = reorder(Gene, P_value), y = -log(P_value))) +
  geom_col(fill = "skyblue") +
  geom_hline(yintercept = -log(0.05), color = "red", linetype = "dashed") +
  labs(title = "Distribution of P-values by Mitochondrial ribosomal gene", x = "Gene", y = "P-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_m, aes(x = reorder(Gene, P_value), y = -log(P_value))) +
  geom_col(fill = "skyblue") +
  geom_hline(yintercept = -log(0.05), color = "red", linetype = "dashed") +
  labs(title = "Distribution of P-values by Mitochondrial gene", x = "Gene", y = "P-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data_r, aes(x = reorder(Gene, P_value), y = -log(P_value))) +
  geom_col(fill = "skyblue") +
  geom_hline(yintercept = -log(0.05), color = "red", linetype = "dashed") +
  labs(title = "Distribution of P-values by Ribosomal gene", x = "Gene", y = "P-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggplot(data_mr, aes(x = P_value)) +
#   geom_density(fill = "purple", alpha = 0.5) +
#   labs(title = "Density Distribution of Ribosomal Gene P-values",
#        x = "P-value", y = "Density") +
#   theme_minimal()
