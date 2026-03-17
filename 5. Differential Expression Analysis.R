library(DESeq2)
library(readr)

##########################################
#read featureCounts text file
##########################################
#load gene count data from featureCounts output
fc <- read.table(
  "/Users/qianzhao/Desktop/academics/thesis/data/aim 3/gene_counts.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE
)

#set gene IDs as row names
rownames(fc) <- fc$Geneid
#extract count columns only
countData <- fc[, 8:ncol(fc)]
#clean sample names
colnames(countData) <- sub("_Aligned.sortedByCoord.out.bam$", "", colnames(countData))
sample_names <- colnames(countData)

#define sample conditions: A1-7 are KP, A8-14 are KPC
condition <- ifelse(
  as.integer(sub("^A", "", sample_names)) <= 7,
  "KP",
  "KPC"
)

#create sample dataframe for DESeq2
colData <- data.frame(
  sample = sample_names,
  condition = factor(condition, levels = c("KP", "KPC"))  #KP is the reference
)
rownames(colData) <- colData$sample

##########################################
#run DESeq2 
##########################################
#create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(
  countData = round(countData),  #use integer counts for DESeq2's negative binomial distribution
  colData   = colData,
  design    = ~ condition
)

#run DESeq2
dds <- DESeq(dds)

#extract results for KPC vs KP (log2FC > 0 means higher in KPC)
res <- results(dds, contrast = c("condition", "KPC", "KP"))
#order genes by adjusted p-value, most significant first
resOrdered <- res[order(res$padj), ]
##########################################
#PCA plot using top 1000 most variable genes
##########################################
#apply variance stabilizing transformation for visualization
vsd <- vst(dds, blind = FALSE)

#perform PCA on top 1000 most variable genes
pcaData <- plotPCA(vsd, intgroup = "condition", ntop = 1000, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#libraries for plot visualization
library(ggplot2)
library(ggrepel)
library(dplyr)

#PCA plot
ggplot(pcaData, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3) +
  xlab(paste0("PC1: ", round(percentVar[1], 1), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2], 1), "% variance")) +
  theme_classic(base_size = 14)

#quantify clustering quality with silhouette score
library(cluster)
df <- data.frame(
  PC1 = pcaData$PC1,
  PC2 = pcaData$PC2,
  condition = pcaData$condition
)

X <- as.matrix(df[, c("PC1","PC2")])
d <- dist(X, method = "euclidean") #euclidean distance between samples in PCA space
sil <- silhouette(as.integer(df$condition), d)
mean(sil[, "sil_width"]) #closer to 1 = better clustering
##########################################
#heatmap of top 1000 most significant genes
##########################################
library(pheatmap)
library(matrixStats)

ann <- as.data.frame(colData(dds)[, "condition", drop = FALSE])
sample_names <- rownames(ann)
#order samples numerically A1 ... A14
sample_order <- sample_names[order(as.numeric(sub("^A", "", sample_names)))]

#select top 1000 genes by adjusted p-value
res <- res[order(res$padj), ]
top <- rownames(res)[1:1000]
#extract expression data for these 1000 genes
mat_ordered <- assay(vsd)[top, sample_order]
ann_ordered <- ann[sample_order, , drop = FALSE]
ann_colors <- list(condition = c(KP = "salmon", KPC = "mediumturquoise"))

#heatmap with row scaling with z-scores
pheatmap(
  mat_ordered,
  scale = "row",
  annotation_col = ann_ordered,
  annotation_colors = ann_colors,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  treeheight_row = 0,
  show_rownames = FALSE
)

#quantify clustering based on silhouette score
X <- t(mat_ordered)
d <- dist(X, method = "euclidean") #euclidean distance
sil <- silhouette(as.integer(ann_ordered$condition), d)
mean(sil[, "sil_width"])
##########################################
#MA plot
##########################################
library(ggrepel)
library(dplyr)
library(AnnotationDbi)
library(org.Mm.eg.db)

df <- as.data.frame(resOrdered)

#extract ENSEMBL ids
ens_ids <- sub("\\..*", "", rownames(resOrdered))
#map ENSEMBL ids to gene symbols
symbols <- mapIds(
  org.Mm.eg.db,
  keys = ens_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
df$symbol <- as.character(symbols)  

#set threshold for mean expression
padj_cutoff <- 0.05
mean_cutoff <- 10
log2FoldChange_cutoff <- 1

#categorize genes based on expression level, significance, fold change
df <- df %>%
  mutate(
    baseMean_plot = pmax(baseMean, 1),
    status = case_when(
        #upregulated in KPC = significant, log2FC >=1, mean expression >= 1e2
        baseMean_plot >= mean_cutoff &
        !is.na(padj) & padj < padj_cutoff &
        log2FoldChange >= log2FoldChange_cutoff  ~ "Upregulated in KPC",
        #downregulated in KPC = significant, log2FC <= -1, mean expression >= 1e2
        baseMean_plot >= mean_cutoff &
        !is.na(padj) & padj < padj_cutoff &
        log2FoldChange <= -1*log2FoldChange_cutoff ~ "Downregulated in KPC",
      TRUE ~ "Filtered Out"
    )
  )

#retrieve genes that passed thresholds
deg_genes <- df %>%
  filter(status %in% c("Upregulated in KPC", "Downregulated in KPC")) %>%
  arrange(padj)

#label with cleaned ensembl ids
deg_genes$ensembl_clean <- sub("\\..*$", "", rownames(deg_genes))

#add label column for plotting
deg_genes$label <- ifelse(is.na(deg_genes$symbol) | deg_genes$symbol == "",
                          deg_genes$ensembl_clean,
                          deg_genes$symbol)

#create dataframe for genes that passed thresholds
deg_genes$has_symbol <- !(is.na(deg_genes$symbol) | deg_genes$symbol == "")
label_df <- deg_genes %>%
  filter(has_symbol)
label_df$symbol_italic <- paste0("italic('", label_df$symbol, "')")

#MA plot showing log2FC vs mean expression
p <- ggplot(df, aes(x = baseMean_plot, y = log2FoldChange)) +
  geom_point(aes(color = status), alpha = 0.6, size = 0.7) +
  geom_point(data = label_df, aes(color = status), size = 2) +
  geom_text_repel(data = label_df, aes(label = symbol),
                  size = 4, fontface = "italic", max.overlaps = Inf) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_hline(yintercept = c(-1*log2FoldChange_cutoff, log2FoldChange_cutoff), linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = mean_cutoff, linetype = "dashed", linewidth = 0.4) +
  scale_x_log10() +
  scale_color_manual(
    values = c(
      "Upregulated in KPC"   = "#D62728",
      "Downregulated in KPC" = "#1F77B4",
      "Filtered Out"    = "grey80"
    ) 
  ) +
  labs(x = "Mean expression (all samples, log10 scale)",
       y = "log2 fold change (KPC vs KP)", color = NULL) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

print(p)
##########################################
#bar plot of differentially expressed genes that satisfied thresholds
##########################################
#add direction column to deg_genes for coloring purposes
deg_genes$direction <- ifelse(deg_genes$log2FoldChange > 0, "Positive", "Negative")

#order genes by log2FC
deg_genes$label <- factor(
  deg_genes$label,
  levels = deg_genes$label[order(deg_genes$log2FoldChange)]
)
deg_genes$lfc_label <- sprintf("%.2f", deg_genes$log2FoldChange)

#format plot
ymin <- min(deg_genes$log2FoldChange)
ymax <- max(deg_genes$log2FoldChange)
pad  <- 0.25 * max(abs(c(ymin, ymax)))

#bar plot of log2 fold changes
ggplot(deg_genes, aes(x = label, y = log2FoldChange, fill = direction)) +
  geom_col() +
  geom_text(
    aes(label = lfc_label),
    hjust = ifelse(deg_genes$log2FoldChange >= 0, -0.15, 1.15),
    size = 4
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(ymin - pad, ymax + pad)) +
  scale_fill_manual(
    values = c("Negative" = "#1F77B4",
               "Positive" = "#D62728")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(5.5, 15, 5.5, 5.5),
    axis.text.y = element_text(face = "italic")  # Make y-axis gene names italic
  ) +
  labs(
    title = "Differentially Expressed Genes",
    x = "Gene Name",
    y = "log2 Fold Change",
    fill = "Direction"
  )

##########################################
#heatmap of genes that satisfied thresholds
##########################################
#build matrix with ENSEMBL rownames
mat <- assay(vsd)
mat_ens <- sub("\\..*$", "", rownames(mat))
rownames(mat) <- mat_ens

#subset the matrix to include only DEGs
mat_deg <- mat[deg_genes$ensembl_clean, sample_order, drop = FALSE]

#set rownames to gene symbols
rownames(mat_deg) <- deg_genes$label

#make the gene names italics
italic_labels <- paste0("italic('", rownames(mat_deg), "')")

#plot heatmap
pheatmap(
  mat_deg,
  scale = "row",
  annotation_col = ann_ordered,
  annotation_colors = ann_colors,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  treeheight_row = 0,
  show_rownames = TRUE,
  labels_row = parse(text = italic_labels),
  fontsize_row = 8  # Adjust font size for readability
)

##########################################
#PCA of genes that satisfied thresholds
##########################################
#perform PCA on selected genes
pca_res <- prcomp(t(mat_deg), center = TRUE, scale. = FALSE)

#calculate percent variance explained by each PC component
percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

#create dataframe for plotting
pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$sample <- rownames(pca_df)

#add condition information
pca_df$condition <- ifelse(
  as.integer(sub("^A", "", pca_df$sample)) <= 7,
  "KP",
  "KPC"
)

#plot PCA results
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(min.segment.length = Inf) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = c("KP" = "salmon", "KPC" = "mediumturquoise")) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

##########################################
#pathway analysis
##########################################
library(clusterProfiler)
library(enrichplot)
library(stringr)

#prepare gene list
df$ensembl_clean <- sub("\\..*$", "", rownames(resOrdered))

#for genes with multiple transcripts, keep the one with the largest abs val of Wald statistic
gene_df <- df |>
  dplyr::filter(!is.na(stat)) |>
  dplyr::group_by(ensembl_clean) |>
  dplyr::summarise(stat = stat[which.max(abs(stat))], .groups="drop")

#create named vector of Wald statistic
geneList <- gene_df$stat
names(geneList) <- gene_df$ensembl_clean
geneList <- sort(geneList, decreasing = TRUE)

#map ENSEMBL ids to ENTREZ ids
map <- bitr(names(geneList),
            fromType = "ENSEMBL",
            toType   = "ENTREZID",
            OrgDb    = org.Mm.eg.db)

#keep only mappable genes
geneList_entrez <- geneList[map$ENSEMBL]
names(geneList_entrez) <- map$ENTREZID

#collapse duplicate ENTREZ (keep max |stat|)
tmp <- data.frame(ENTREZID = names(geneList_entrez),
                  stat    = as.numeric(geneList_entrez),
                  stringsAsFactors = FALSE)

tmp <- tmp |>
  dplyr::group_by(ENTREZID) |>
  dplyr::summarise(stat = stat[which.max(abs(stat))], .groups="drop")

geneList_entrez <- tmp$stat
names(geneList_entrez) <- tmp$ENTREZID

#ensure numeric and sorted
geneList_entrez <- as.numeric(geneList_entrez)
names(geneList_entrez) <- tmp$ENTREZID
geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)

#GSEA: GO
gsea_go_bp <- gseGO(
  geneList     = geneList_entrez,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP", #biological process ontology
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

#create dotplot of top GO results
dotplot(gsea_go_bp, showCategory = 20) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 50)
  )

#GSEA: KEGG
gsea_kegg <- gseKEGG(
  geneList     = geneList_entrez,
  organism     = "mmu",
  minGSSize    = 10,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

#create dotplot of top KEGG results
dotplot(gsea_kegg, showCategory = 20) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 50)
  )

#GSEA: Reactome
library(ReactomePA)

gsea_reactome <- gsePathway(
  geneList     = geneList_entrez,
  organism     = "mouse",
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

#create dotplot of top Reactome results
dotplot(gsea_reactome, showCategory = 10) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 50)
  )

