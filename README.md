# 🔬 Single-Cell RNA-Seq Analysis of Breast Cancer | Seurat + Harmony + Monocle

This project presents an end-to-end **single-cell transcriptomic analysis** of breast cancer using **Seurat**, **Harmony**, **SingleR**, and **Monocle** in R. It integrates tumor and normal breast tissue datasets, performs quality control, cell clustering, differential gene expression, cell-type annotation, functional enrichment, and pseudotime trajectory analysis.

---

## 🎯 Objective

To analyze and compare cellular heterogeneity between tumor and normal breast tissues, identify differentially expressed genes (DEGs), annotate cell types, and explore developmental trajectories via pseudotime.

---

## 📚 Dataset

- **Input**: Gene expression matrices (`normal.txt`, `tumor.txt`)
- **Samples**: scRNA-seq from breast cancer and adjacent normal tissue
- **Format**: Tab-delimited count tables

---

## 🛠️ Tools & Packages

- `Seurat`, `Harmony`, `SingleR`, `Monocle`, `slingshot`
- `celldex`, `clusterProfiler`, `org.Mm.eg.db`, `enrichplot`
- `ggplot2`, `pheatmap`, `viridis`, `scRNAseq`, `TabulaMurisData`

---

## 🔄 Analysis Workflow

### 🧹 1. Preprocessing & QC
- Imported normal/tumor data into Seurat
- Visualized nGenes/UMIs & % mitochondrial and rRNA content
- Filtered low-quality cells (genes <200, %MT >5%)

### 🔗 2. Integration & Normalization
- Merged normal + tumor cells
- Performed log-normalization and scaling
- Selected highly variable genes

### ⚗️ 3. Dimensionality Reduction & Batch Correction
- Performed PCA and visualized ElbowPlot
- Used **Harmony** for batch correction across conditions
- Visualized **UMAP** embeddings

### 🔎 4. Clustering & Annotation
- Identified clusters (resolution = 0.5)
- Annotated cells with **SingleR** using Tabula Muris (spleen) reference
- Visualized cell types in UMAP plots

### 🧬 5. Differential Expression
- Compared tumor vs. normal clusters
- Extracted top 20 upregulated & downregulated genes
- Visualized using bar plots, heatmaps, violin plots
- Saved results as CSV

### 🧠 6. Functional Enrichment (GO & GSEA)
- Used `clusterProfiler` and `org.Mm.eg.db` for:
  - GO Biological Process enrichment
  - Dotplots & cnetplots for top genes

### ⏱️ 7. Pseudotime Trajectory Inference
- Converted Seurat → SingleCellExperiment → Monocle object
- Used **Slingshot** to infer trajectories
- Visualized tumor vs. normal developmental paths

---

## 📊 Output Files (Examples)

| File | Description |
|------|-------------|
| `top20_upregulated_genes.csv` | Top 20 tumor genes |
| `top20_downregulated_genes.csv` | Top 20 normal genes |
| `top_differentially_expressed_genes.csv` | Combined up/down |
| `*.pdf` | UMAPs, barplots, dotplots, heatmaps |
| `breast.Rdata`, `breast.Rhistory` | Workspace/image |

---

## 🧠 Skills Demonstrated

- Single-cell data preprocessing & QC  
- Batch effect correction with Harmony  
- Cell clustering & annotation using Seurat + SingleR  
- Differential gene expression analysis  
- Functional enrichment (GO) and gene set analysis  
- Trajectory inference with Monocle + Slingshot  
- Custom ggplot2 visualizations

---
