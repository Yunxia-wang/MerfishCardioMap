# =============================================================================
# Single-cell RNA-seq Analysis Workflow with Seurat
#
# This script is organized into three main sections:
# 1. Setup and Helper Functions
# 2. Analysis of the LPM Subclusters
# 3. Integrated Analysis of Combined Samples
# 4. Analysis of the Combined Annotated Samples
#
# The `nkx2.5_project_R_functions.R` script is sourced for custom functions.
# =============================================================================


# =============================================================================
# Section 1: Setup and Helper Functions
# =============================================================================

# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

# It's generally better to not remove packages unless there's a specific conflict.
# The user's original code included this, but it's commented out for a cleaner workflow.
# remove.packages("ggplot2")
# install.packages("https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_4.3.0.tar.gz", repos = NULL, type = "source")
# BiocManager::install("ggplot2")
# install.packages("ggplot2")
# remotes::install_github("satijalab/seurat@v4.3.0")

# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork) # For combining plots
library(biomaRt) # For the user's custom R function

# Set options for plot display
options(repr.plot.width=15, repr.plot.height=9)

# Define a function to perform UMAP flipping and plotting with consistent styling
plot_flipped_umap <- function(seurat_obj,
                              reduction = "umap",
                              pt.size = 0.1,
                              custom_colors = NULL,
                              title = "",
                              xlim = c(NA, NA),
                              ylim = c(NA, NA)) {
  
  # Check if the reduction exists
  if (!reduction %in% Reductions(seurat_obj)) {
    stop(paste("Reduction", reduction, "not found in the Seurat object."))
  }

  # Extract and flip UMAP coordinates
  umap_coords <- Embeddings(seurat_obj, reduction = reduction)
  umap_coords[, 1] <- -umap_coords[, 1]

  # Create a temporary object with the flipped UMAP for plotting
  seurat_obj_flipped <- seurat_obj
  seurat_obj_flipped[[reduction]] <- CreateDimReducObject(
    embeddings = umap_coords, 
    key = paste0(toupper(reduction), "_"), 
    assay = DefaultAssay(seurat_obj)
  )

  # Generate the base plot
  p <- DimPlot(seurat_obj_flipped, reduction = reduction, pt.size = pt.size)

  # Apply custom colors if provided
  if (!is.null(custom_colors)) {
    p <- p + scale_color_manual(values = custom_colors)
  }
  
  # Apply consistent theme and axis limits
  p <- p + theme(
    axis.line = element_line(size = 1),
    axis.ticks.x = element_line(size = 1, lineend = "round"),
    axis.ticks.y = element_line(size = 1, lineend = "round"),
    axis.ticks.length = unit(0.08, "inches"),
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold")
  )
  
  # Set axis limits if provided
  if (!is.na(xlim[1]) && !is.na(xlim[2])) {
    p <- p + xlim(xlim)
  }
  if (!is.na(ylim[1]) && !is.na(ylim[2])) {
    p <- p + ylim(ylim)
  }

  # Add a title if provided
  if (nchar(title) > 0) {
    p <- p + labs(title = title)
  }

  return(p)
}

# Define a function to save a plot to both PNG and PDF formats
save_plot <- function(plot, filename, width = 6, height = 6, dpi = 600) {
  # Save as PNG
  ggsave(
    paste0("../result/", filename, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
  # Save as PDF
  ggsave(
    paste0("../result/", filename, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}

# =============================================================================
# Section 2: Analysis of the LPM Subclusters
# =============================================================================

# Load the annotated LPM subclusters Seurat object
seurat_lpm_obj <- readRDS("../data/LPM_subclusters_annoated.rds")

# Plot and save the UMAP for all LPM subclusters
p_lpm_merged <- plot_flipped_umap(seurat_lpm_obj,
                                  xlim = c(-5, 2.5),
                                  ylim = c(-5, 5))
print(p_lpm_merged)
save_plot(p_lpm_merged, "LPM_umap_merged")

# Subset the LPM data for 'somites_20' sample
seurat_20ss_lpm <- subset(seurat_lpm_obj, subset = sample == "somites_20")

# Plot and save the UMAP for the subsetted data
p_20ss_lpm <- plot_flipped_umap(seurat_20ss_lpm,
                                xlim = c(-5, 2.5),
                                ylim = c(-5, 5))
print(p_20ss_lpm)
save_plot(p_20ss_lpm, "LPM_umap_somites_20")

# Plot and save a split UMAP for the LPM data
pdf("../result/LPM_separated_clustering.pdf", height = 4, width = 8)
DimPlot(seurat_lpm_obj, reduction = "umap", split.by = "sample", pt.size = 0.1) +
  coord_cartesian(xlim = c(-5, 2.5), ylim = c(-6, 4))
dev.off()

# Plot and save a feature plot for LPM genes
DefaultAssay(seurat_lpm_obj) <- "RNA"
p_lpm_genes <- FeaturePlot(seurat_lpm_obj, features = c("myh7", "myh6"), split.by = "sample", pt.size = 0.1, combine = FALSE, min.cutoff = 0, max.cutoff = 5)
plots_lpm <- lapply(p_lpm_genes, function(x) x + xlim(c(-5, 2.5)) + ylim(c(-6, 4)))
pdf('../result/LPM_gene_expression.pdf', height = 6, width = 8)
wrap_plots(plots_lpm)
dev.off()


# =============================================================================
# Section 3: Integrated Analysis of Combined Samples
# =============================================================================

# The following functions are assumed to be in the sourced R script.
# source("../data/nkx2.5_project_R_functions.R")
# These functions were present in the original user code and are necessary for
# the integrated analysis to run.
# convertGeneID2Symbol()
# removeRP()
# generateSeuratObject()

# Load and process raw count matrices
count_matrix <- list(
  convertGeneID2Symbol("../data/Z5EL-014.csv"),
  convertGeneID2Symbol("../data/Z5EL-014_2.csv"),
  convertGeneID2Symbol("../data/20_Somites_1.csv"),
  convertGeneID2Symbol("../data/20_Somites_2.csv")
)
count_matrix <- lapply(count_matrix, removeRP, RPGene = RPGene)

# Generate Seurat objects
seu_14ss_1 <- generateSeuratObject(count_matrix[[1]], "somites_14", 1)
seu_14ss_2 <- generateSeuratObject(count_matrix[[2]], "somites_14", 2)
seu_20ss_1 <- generateSeuratObject(count_matrix[[3]], "somites_20", 1)
seu_20ss_2 <- generateSeuratObject(count_matrix[[4]], "somites_20", 2)
rm(count_matrix)

# Create a list of Seurat objects for integration
nkx25.list <- list(seu_14ss_1, seu_14ss_2, seu_20ss_1, seu_20ss_2)

# Normalize and find variable features for each object
nkx25.list <- lapply(X = nkx25.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select integration features and perform integration
features <- SelectIntegrationFeatures(object.list = nkx25.list)
anchors <- FindIntegrationAnchors(object.list = nkx25.list, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)

# Set the default assay and run the standard Seurat workflow
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.1)

# QC plot
pdf("../result/QC.pdf", height = 4, width = 6)
VlnPlot(combined, features = c("percent.mt", "nCount_RNA", "nFeature_RNA"), group.by = "sample", pt.size = 0)
dev.off()

# Annotation of cell types
combined <- subset(combined, idents = 0:5)
Idents(combined) <- combined@meta.data$seurat_clusters
cluster.ids <- c("LPM", "Neural", "CM", "EC", "Epidermis", "Tailbud")
names(cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, cluster.ids)

# Find all cluster markers
DefaultAssay(combined) <- "RNA"
cluster.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(cluster.markers, file="../result/integrated_cluster_markers.txt", quote=F, sep="\t", row.names=F, col.names=T)

# Visualization
pdf("../result/integrated_clustering_colored_by_sample.pdf", height = 4, width = 6)
DimPlot(combined, reduction = "umap", group.by = "sample")
dev.off()
pdf("../result/integrated_clustering_colored_by_cell_type.pdf", height = 4, width = 6)
DimPlot(combined, reduction = "umap", label = FALSE)
dev.off()
pdf("../result/separated_clustering.pdf", height = 4, width = 8)
DimPlot(combined, reduction = "umap", split.by = "sample")
dev.off()


# =============================================================================
# Section 4: Analysis of the Combined Annotated Samples
# =============================================================================

# Load the combined annotated Seurat object
seurat_combined_obj <- readRDS("../data/combined_samples_annotated.rds")

# Define custom color palette for cell types
custom_colors <- c(
  "Tailbud"   = "#ec6cc7",
  "EC"        = "#55bcc2",
  "CM"        = "#4DAF4A",
  "LPM"       = "#e87e71",
  "Epidermis" = "#6f9bf8",
  "Neural"    = "#c59832",
  "LPM-3"     = "#984EA3",
  "LPM-1"     = "#FF7F00",
  "LPM-2"     = "#A65628"
)

# Plot and save UMAP for the full combined object
p_full <- plot_flipped_umap(seurat_combined_obj)
print(p_full)
save_plot(p_full, "LPM_umap", width = 6, height = 6)

# Subset and plot for 'somites_14'
seurat_14ss <- subset(seurat_combined_obj, subset = sample == "somites_14")
p_14ss <- plot_flipped_umap(seurat_14ss, custom_colors = custom_colors)
print(p_14ss)
save_plot(p_14ss, "14SS_cluster_umap", width = 6, height = 6)

# Subset and plot for 'somites_20'
seurat_20ss <- subset(seurat_combined_obj, subset = sample == "somites_20")
p_20ss <- plot_flipped_umap(seurat_20ss, custom_colors = custom_colors)
print(p_20ss)
save_plot(p_20ss, "20SS_cluster_umap", width = 6, height = 6)

# Additional plots for somites_20 with different limits
p_20ss_lim <- plot_flipped_umap(seurat_20ss,
                                custom_colors = custom_colors,
                                xlim = c(-13, 13),
                                ylim = c(-10, 10))
save_plot(p_20ss_lim, "somites_20_integrated_clustering_colored_by_cell_type", width = 7, height = 6)

# Plot separated UMAP
pdf("../result/separated_clustering.pdf", height = 5, width = 10)
DimPlot(seurat_combined_obj, reduction = "umap", split.by = "sample", pt.size = 0.1) +
  coord_cartesian(xlim = c(-13, 13), ylim = c(-10, 10))
dev.off()

# Feature plots for combined samples
DefaultAssay(seurat_combined_obj) <- "RNA"
p_combined_genes <- FeaturePlot(seurat_combined_obj, features = c("myh7", "myh6"), split.by = "sample", pt.size = 0.1, combine = FALSE, min.cutoff = 0, max.cutoff = 5)
plots_combined <- lapply(p_combined_genes, function(x) x + xlim(c(-13, 13)) + ylim(c(-10, 10)))
pdf('../result/gene_expression.pdf', height = 9, width = 12)
wrap_plots(plots_combined)
dev.off()
