#### Single-cell RNA analysis of Breast Cancer and save workplace####

save.image("breast.Rdata")
savehistory("breast.Rhistory")
loadhistory("breast.Rdata")
load("breast.Rdata")

#### Set seed for reproducibility####
set.seed(12345)

####Load packages####
install.packages("Seurat")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("celldex")

BiocManager::install("celldex", update = TRUE)
install.packages("pheatmap")
install.packages("viridis")

install.packages("devtools")
library(devtools)
devtools::install_github("immunogenomics/harmony")
install.packages("harmony")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("monocle", "clusterProfiler", "org.Hs.eg.db", "enrichplot"))

install.packages("immgen")

library(monocle)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(celldex)
library(SingleR)
library(celldex)
library(tidyverse)
library(Seurat)
library(pheatmap)
library(viridis)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(harmony)
library(immgen)

####Load data ####

# Define paths
path_normal <- "E:/Breast_Cancer/Breast-Cancer/normal.txt"
path_tumor <- "E:/Breast_Cancer/Breast-Cancer/tumor.txt"

# Load data using read.table (assuming tab-delimited gene expression matrices)
normal.data <- read.table(path_normal, header=TRUE, row.names=1, sep="\t")
tumor.data <- read.table(path_tumor, header=TRUE, row.names=1, sep="\t")

# Create Seurat objects
normal <- CreateSeuratObject(counts = normal.data)
tumor <- CreateSeuratObject(counts = tumor.data)

#### Pre-processing (Quality check)####

## Check droplets

# Visualize number of features (genes) and UMIs per cell
par(mfrow=c(2,2))
hist(normal$nFeature_RNA, breaks=50, main="Normal - No. of Genes per Cell", xlab="# Genes", col="lightblue")
hist(tumor$nFeature_RNA, breaks=50, main="Tumor - No. of Genes per Cell", xlab="# Genes", col="salmon")
hist(normal$nCount_RNA, breaks=50, main="Normal - No. of UMIs per Cell", xlab="# UMIs", col="lightblue")
hist(tumor$nCount_RNA, breaks=50, main="Tumor - No. of UMIs per Cell", xlab="# UMIs", col="salmon")

# Check for mouse rRNA genes
rrna_genes <- c("Rn45s", "Rn5.8s", "Rn18s", "Rn28s")

# Calculate the percentage of rRNA for each dataset
normal[["percent.rrna"]] <- PercentageFeatureSet(normal, pattern = paste(rrna_genes, collapse="|"))
tumor[["percent.rrna"]] <- PercentageFeatureSet(tumor, pattern = paste(rrna_genes, collapse="|"))

# Visualization
VlnPlot(normal, features = c("percent.rrna"), ncol = 1)
VlnPlot(tumor, features = c("percent.rrna"), ncol = 1)

# Identify mitochondrial genes
mito.genes <- grep("^mt-", rownames(normal), value = TRUE)

# Calculate the percentage of mitochondrial transcripts for the normal dataset
normal$percent.mt <- PercentageFeatureSet(normal, features = mito.genes)

# Repeat the process for the tumor dataset
mito.genes.tumor <- grep("^mt-", rownames(tumor), value = TRUE)
tumor$percent.mt <- PercentageFeatureSet(tumor, features = mito.genes.tumor)

# To check the percentage 
head(normal@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt")])
head(tumor@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt")])

# For the "normal" dataset
VlnPlot(normal, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

# For the "tumor" dataset
VlnPlot(tumor, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

# Now filter assigned to a new variable called filtered_normal and tumor
filtered_normal <- subset(normal, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)

filtered_tumor <- subset(tumor, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)

# Number of cells before and after filtering for the normal dataset
print(paste("Normal before filtering:", dim(normal)[2]))
print(paste("Normal after filtering:", dim(filtered_normal)[2]))

# Number of cells before and after filtering for the tumor dataset
print(paste("Tumor before filtering:", dim(tumor)[2]))
print(paste("Tumor after filtering:", dim(filtered_tumor)[2]))

#### Merge the data sets####

seurat.combined <- merge(x = filtered_normal, y = filtered_tumor, add.cell.ids = c("normal", "tumor"), project = "combined")


#### Normalize the Data####
seurat.combined <- NormalizeData(seurat.combined)


####Identify Variable Genes####

seurat.combined <- FindVariableFeatures(seurat.combined)


####Scale the Data####

seurat.combined <- ScaleData(seurat.combined, features = VariableFeatures(object = seurat.combined))


####Perform PCA on Individual Datasets####

seurat.combined <- RunPCA(seurat.combined, features = VariableFeatures(object = seurat.combined))


#Find the dimentionality using elbowplot
# Elbow plot for filtered_normal
ElbowPlot(seurat.combined) + ggtitle("Elbow Plot - Combined")

#### Run Harmony for batch correction####

pca_embeddings <- Embeddings(seurat.combined, reduction = "pca")

##check if "orig.ident" exists
"orig.ident" %in% colnames(seurat.combined@meta.data)

colnames(seurat.combined@meta.data)

## Use pca_embeddings in the HarmonyMatrix
harmony_embeddings <- HarmonyMatrix(
  data_mat = pca_embeddings,
  meta_data = seurat.combined@meta.data,
  vars_use = "orig.ident",
  max.iter.cluster = 100
)

##Store the Harmonized Embeddings in Seurat Object

seurat.combined[["harmony"]] <- CreateDimReducObject(embeddings = harmony_embeddings, key = "harmony_")

####Run UMAP on the harmonized data####
#Reductions(seurat.combined)


seurat.combined <- RunUMAP(seurat.combined, features = VariableFeatures(seurat.combined))
DimPlot(seurat.combined, group.by = "orig.ident", reduction = "umap")

##Change label
seurat.combined@meta.data$orig.ident[seurat.combined@meta.data$orig.ident == "pymt"] <- "tumor"
seurat.combined@meta.data$orig.ident[seurat.combined@meta.data$orig.ident == "wt"] <- "normal"
DimPlot(seurat.combined, group.by = "orig.ident", reduction = "umap")

#### Perform Clustering####
##Find neighbors
seurat.combined <- FindNeighbors(seurat.combined, dims = 1:20)
##Find clusters
seurat.combined <- FindClusters(seurat.combined, resolution = 0.5)
DimPlot(seurat.combined, group.by = "seurat_clusters", reduction = "umap")

##################################################################################
#### Cell annotation####

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SingleR")

library(SingleR)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scRNAseq")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TabulaMurisData")

library(ExperimentHub)

exp_hub <- ExperimentHub()

query(exp_hub, "TabulaMurisData")

exp_hub[["EH1617"]]

##Subset based on spleen

spleen_ref <- exp_hub[["EH1617"]]

spleen_ref <- spleen_ref[, spleen_ref$tissue == "Spleen"]

spleen_ref <- spleen_ref[, !is.na(spleen_ref$cell_ontology_class)]

##Install scuttle


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(scuttle)

##Log normal spleen_ref
spleen_ref <- logNormCounts(spleen_ref)

seurat_as_sce <- Seurat::as.SingleCellExperiment(seurat.combined)


singleR_results <- SingleR(test = as.SingleCellExperiment(seurat.combined), ref = spleen_ref, labels = spleen_ref$cell_ontology_class)

seurat.combined$singlr_label <- singleR_results$labels

DimPlot(seurat.combined, reduction = "umap", group.by = "singlr_label", label = TRUE)


head(singleR_results)

seurat.combined$SingleR <- singleR_results$labels

table(seurat.combined$SingleR)


manual_colors <- c("B cell" = "#E41A1C", 
                   "dendritic cell" = "#377EB8", 
                   "macrophage" = "#4DAF4A", 
                   "natural killer cell" = "#984EA3", 
                   "T cell" = "#FF7F00")

DimPlot(seurat.combined, group.by = "SingleR", reduction = "umap", cols = manual_colors)


manual_colors <- c("B cell" = "#E41A1C", 
                   "dendritic cell" = "#377EB8", 
                   "macrophage" = "#4DAF4A", 
                   "natural killer cell" = "#984EA3", 
                   "T cell" = "#FF7F00")

umap_data <- seurat.combined@reductions$umap@cell.embeddings
plot(umap_data[,1], umap_data[,2], pch=20, col=manual_colors[seurat.combined$SingleR], asp=1, xlab="UMAP1", ylab="UMAP2")

legend("topright", legend=names(manual_colors), fill=manual_colors, border="white", bty="n", cex=0.8)





####Differential Expression Analysis####
##Find markers between conditions

unique(seurat.combined@meta.data$orig.ident)

seurat.combined <- SetIdent(seurat.combined, value = "orig.ident")

tumor_vs_normal_markers <- FindMarkers(seurat.combined, ident.1 = "tumor", ident.2 = "normal", min.pct = 0.25, logfc.threshold = 0.25)

#Sort by avg_log2FC
sorted_markers <- tumor_vs_normal_markers[order(-tumor_vs_normal_markers$avg_log2FC), ]

head(sorted_markers, n = 20)


####To identify the top 20 upregulated and downregulated genes in the tumor compared to normal####

# Sorting the markers by avg_log2FC to get top upregulated genes in tumor
top20_upregulated <- head(tumor_vs_normal_markers[order(-tumor_vs_normal_markers$avg_log2FC),], 20)

# Sorting the markers by avg_log2FC to get top downregulated genes in tumor
top20_downregulated <- head(tumor_vs_normal_markers[order(tumor_vs_normal_markers$avg_log2FC),], 20)

print(top20_upregulated)
print(top20_downregulated)

#View the plot
VlnPlot(seurat.combined, features = rownames(top20_upregulated))
VlnPlot(seurat.combined, features = rownames(top20_downregulated))


# Saving top 20 upregulated genes
write.csv(top20_upregulated, file = "top20_upregulated_genes.csv")

# Saving top 20 downregulated genes
write.csv(top20_downregulated, file = "top20_downregulated_genes.csv")

################### Cell annotation###########

library(ggplot2)
##First extract and combine the top 20 up and dwon regualted
# Extract top 20 upregulated genes
top_upregulated <- head(sorted_markers, n=20)

# Extract top 20 downregulated genes
top_downregulated <- tail(sorted_markers, n=20)

# Combine them for visualization
combined_genes <- rbind(top_upregulated, top_downregulated)

# Combine the top upregulated and downregulated genes for saving
top_genes_df <- rbind(top_upregulated, top_downregulated)

# Write to CSV
write.csv(top_genes_df, file="top_differentially_expressed_genes.csv", row.names=FALSE)


# Plot
head(combined_genes)

ggplot(combined_genes, aes(x=reorder(GeneName, avg_log2FC), y=avg_log2FC, fill=avg_log2FC)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_gradient2(low="green", mid="white", high="red", midpoint=0) +
  labs(title="Top Differentially Expressed Genes", x="Genes", y="Average Log2 Fold Change") +
  theme_minimal()

#Heatmap
top_genes <- rownames(combined_genes)

DoHeatmap(seurat.combined, features = top_genes) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low="green", mid="white", high="red", midpoint=0)


#### Functional Enrichment Analysis ####

# Using the `clusterProfiler` package for GSEA:
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("org.Mm.eg.db")

library(org.Mm.eg.db)
library(clusterProfiler)
enrich_result <- enrichGO(gene = rownames(top20_upregulated), universe = rownames(tumor_vs_normal_markers), 
                          OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
barplot(enrich_result)


enrich_result_1 <- enrichGO(gene = rownames(top20_downregulated), universe = rownames(tumor_vs_normal_markers), 
                          OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
barplot(enrich_result_1)


####Gene Set Enrichment Analysis (GSEA) additiional analysis####

#Dotplot Visualization
dotplot(enrich_result, showCategory=15) # for upregulated genes
dotplot(enrich_result_1, showCategory=15) # for downregulated genes


#Cnetplot Visualization

cnetplot(enrich_result, foldChange = top20_upregulated$avg_log2FC)
cnetplot(enrich_result_1, foldChange = top20_downregulated$avg_log2FC)


#### Feature gene plot for top5 unregulated ####

genes.of.interest <- c("Ifitm1", "Wfdc172", "BC100530", "Hbb-bs", "Gm5483")
Seurat::FeaturePlot(object = seurat.combined, features = genes.of.interest, reduction = "umap")


#### Trajectory analysis or Pseudotime Analysis####

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("monocle")
library(monocle)

##Convert the Seurat object to a Monocle CellDataSet (CDS) object

# Convert the Seurat object to a Monocle CDS
library(SingleCellExperiment)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("slingshot")


sce <- Seurat::as.SingleCellExperiment(seurat.combined)
# Extract the Seurat clusters
sce$cluster <- seurat.combined$seurat_clusters
library(slingshot)
sds <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'UMAP')


library(ggplot2)

# Extract UMAP data and add Slingshot results
umap_data <- as.data.frame(reducedDims(sce)$UMAP)
colnames(umap_data) <- c("UMAP1", "UMAP2")



# Extract tumor cells
sce_tumor <- sce[sce$condition == "tumor",]

# Extract normal cells
sce_normal <- sce[sce$condition == "normal",]

# For tumor:
start.clus_tumor <- 1  # Choose based on your knowledge
sds_tumor <- slingshot(sce_tumor, clusterLabels = 'seurat_clusters', start.clus = start.clus_tumor, reducedDim = 'UMAP')

# For normal:
start.clus_normal <- 1  # Choose based on your knowledge
sds_normal <- slingshot(sce_normal, clusterLabels = 'seurat_clusters', start.clus = start.clus_normal, reducedDim = 'UMAP')

# Visualizing tumor trajectory
# Plot UMAP for tumor cells
plot(reducedDims(sds_tumor)$UMAP, 
     col = sds_tumor$seurat_clusters, 
     pch=16, 
     cex=0.5, 
     xlab = "UMAP1", 
     ylab = "UMAP2",
     main = "Tumor Trajectory")

# Add trajectory for tumor
curve_set <- slingCurves(sds_tumor)[[1]]
for (i in seq_along(curve_set)) {
  lines(curve_set[[i]], col='red', lwd=2)
}

# Add cluster labels
text(reducedDims(sds_tumor)$UMAP, labels = sds_tumor$seurat_clusters, cex=0.8, pos=3)


























