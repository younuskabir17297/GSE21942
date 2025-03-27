#Install packages
packages <- c("ggplot2", "ggrepel", "RColorBrewer", "umap", "Rtsne", 
              "tidyverse", "reshape2", "GEOquery", "limma", "pheatmap", 
              "org.Hs.eg.db", "clusterProfiler", "hgu133plus2.db")

new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "pheatmap", "org.Hs.eg.db", 
                       "clusterProfiler", "hgu133plus2.db"))

# Load necessary libraries
library(GEOquery)
library(limma)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(tidyverse)
library(umap)
library(Rtsne)
library(hgu133plus2.db)
library(reshape2)

# Download and extract data
gse <- getGEO("GSE21942", GSEMatrix = TRUE)[[1]]

# Extract expression matrix and metadata
expr_data <- exprs(gse)
metadata <- pData(gse) %>% 
  select(title, geo_accession, `disease state:ch1`) %>% 
  rename(diagnosis = `disease state:ch1`) %>%
  mutate(diagnosis = as.factor(diagnosis))

# 🔹 (OPTIMIZATION) - Set rownames for better referencing
rownames(metadata) <- metadata$geo_accession

# 🔹 (OPTIMIZATION) - Improved filtering of low-expression genes
keep <- rowSums(expr_data > 5) > (0.5 * ncol(expr_data)) 
expr_filtered <- expr_data[keep, ]
expr_filtered <- expr_filtered[!apply(expr_filtered, 1, var) == 0, ]  # Remove zero variance genes
expr_filtered <- log2(expr_filtered + 1)  # Log transformation

# Identify disease status column
metadata$Group <- ifelse(metadata$diagnosis == "healthy", "Control", "MS")
metadata$Group <- factor(metadata$Group, levels = c("Control", "MS"))

# Define design matrix for differential expression analysis
design <- model.matrix(~ 0 + metadata$Group)
colnames(design) <- levels(metadata$Group)

# Perform differential expression analysis
fit <- lmFit(expr_filtered, design)
contrast_matrix <- makeContrasts(MS_vs_Control = MS - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
deg_results <- topTable(fit2, coef = "MS_vs_Control", adjust.method = "fdr", number = Inf)

# 🔹 (OPTIMIZATION) - Keep only necessary columns
deg_results <- deg_results[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
write.csv(deg_results, "DEGs_MS_vs_Control.csv", row.names = TRUE)

# 🔹 (OPTIMIZATION) - Batch effect correction
expr_filtered <- removeBatchEffect(expr_filtered, batch = metadata$geo_accession)

# Gene annotation and enrichment analysis
gene_symbols <- mapIds(hgu133plus2.db, keys = rownames(deg_results), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
gene_symbols <- na.omit(gene_symbols)
entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
entrez_ids <- na.omit(entrez_ids)

# GO Enrichment Analysis
ego <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",  
                pAdjustMethod = "fdr", pvalueCutoff = 0.05, qvalueCutoff = 0.05)



#Figure 1 :
# KEGG Enrichment Analysis

ego_kegg <- enrichKEGG(gene = entrez_ids, organism = "hsa", pvalueCutoff = 0.05)
barplot(ego_kegg, showCategory = 10, title = "Top 10 Enriched KEGG Pathways")
write.csv(as.data.frame(ego_kegg), "KEGG_Enrichment_Results.csv", row.names = FALSE)

# Figure 2
#Heatmap of Top 50 Most Variable Genes
top_genes <- names(sort(apply(expr_filtered, 1, sd), decreasing = TRUE)[1:50])
pheatmap(expr_filtered[top_genes, ], 
         annotation_col = metadata[, "diagnosis", drop = FALSE],
         show_rownames = FALSE, 
         main = "Heatmap of Top 50 Variable Genes")



#Figure 3 
#Perplexity adjustment for t-SNE
perplexity_value <- max(5, min(30, floor(ncol(expr_filtered) / 3)))
tsne_res <- Rtsne(t(expr_filtered), perplexity = perplexity_value, verbose = TRUE, max_iter = 500)
tsne_df <- as.data.frame(tsne_res$Y) %>% mutate(diagnosis = metadata$diagnosis)
ggplot(tsne_df, aes(x = V1, y = V2, color = diagnosis)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = paste("t-SNE (Perplexity =", perplexity_value, ")")) +
  theme_minimal()


#Figure 4
# UMAP analysis
umap_res <- umap(t(expr_filtered))
umap_df <- as.data.frame(umap_res$layout) %>% mutate(diagnosis = metadata$diagnosis)
ggplot(umap_df, aes(x = V1, y = V2, color = diagnosis)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "UMAP: MS vs. Controls") +
  theme_minimal()

# Figure 5 
# Volcano Plot with Gene Labels

# Define threshold for significance
deg_results$threshold <- ifelse(deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1, "Significant", "Not Significant")

# Select top genes to label
top_labels <- deg_results %>%
  filter(threshold == "Significant") %>%
  arrange(adj.P.Val) %>%
  head(20)  # Select top 20 most significant genes

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("red", "black")) +
  geom_text_repel(data = top_labels, aes(label = rownames(top_labels)), size = 3, max.overlaps = 10) +
  theme_minimal() +
  labs(title = "Volcano Plot: MS vs. Control", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value")


#Figure 6 
#Boxplot


# Select top differentially expressed genes safely
top_upregulated <- head(rownames(deg_results[deg_results$logFC > 1, ]), 10)
top_downregulated <- head(rownames(deg_results[deg_results$logFC < -1, ]), 10)
top_genes_to_plot <- c(top_upregulated, top_downregulated)

# Ensure genes exist in the expression dataset
top_genes_to_plot <- intersect(top_genes_to_plot, rownames(expr_filtered))
expr_plot <- expr_filtered[top_genes_to_plot, , drop = FALSE]

# Reshape data for ggplot
expr_plot_melt <- melt(as.matrix(expr_plot), varnames = c("Gene", "Sample"), value.name = "Expression")

# Ensure correct grouping by matching sample names
expr_plot_melt$Group <- metadata$Group[match(expr_plot_melt$Sample, rownames(metadata))]

# Remove any unmatched (NA) samples
expr_plot_melt <- na.omit(expr_plot_melt)

# Plot with ggplot, ensuring color is correctly mapped
ggplot(expr_plot_melt, aes(x = Gene, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, color = "black") +  # Black box outlines
  geom_jitter(aes(color = Group), width = 0.2, size = 1.5, alpha = 0.7) +  # Add jitter for visibility
  coord_flip() +
  theme_minimal() +
  labs(title = "Expression of Top DEGs", x = "Gene", y = "Expression Level") +
  scale_fill_manual(values = c("Control" = "red", "MS" = "blue")) +  # Fill for boxes
  scale_color_manual(values = c("Control" = "red", "MS" = "blue"))  # Color for jitter points


#Figure 7 
#Violin Plot


ggplot(expr_plot_melt, aes(x = Gene, y = Expression, fill = Group)) +
  geom_violin(alpha = 0.6, color = "black", trim = TRUE) +  # Violin plot to show distributions
  geom_jitter(aes(color = Group), width = 0.3, size = 1.5, alpha = 0.8) +  # Add jitter for better visibility
  coord_flip() +  # Flip x and y axes
  theme_minimal() +
  labs(title = "Violin Plot: Expression of Top DEGs", x = "Gene", y = "Expression Level") +
  scale_fill_manual(values = c("Control" = "red", "MS" = "blue")) +  # Fill colors for violin
  scale_color_manual(values = c("Control" = "red", "MS" = "blue")) +  # Color for jitter points
  theme(legend.position = "top")  # Move legend to the top


# Save DEGs (Differentially Expressed Genes) results to CSV
write.csv(deg_results, "DEGs_MS_vs_Control.csv", row.names = TRUE)

# Check the current working directory
getwd()

# Set the working directory (if needed)
setwd("D:/Stata/pcaaaaaaaaaaa")
write.csv(deg_results, "DEGs_MS_vs_Control.csv", row.names = TRUE)


# Save KEGG Enrichment Analysis results to CSV
write.csv(as.data.frame(ego_kegg), "KEGG_Enrichment_Results.csv", row.names = FALSE)

# Save GO Enrichment Analysis results to CSV
getwd()
setwd("D:/Stata/pcaaaaaaaaaaa")
write.csv(as.data.frame(ego), "GO_Enrichment_Results.csv", row.names = FALSE)
