library("Seurat")
library("ggplot2")
library("tidyverse")
library("gridExtra")
library("dplyr")
library("patchwork")
library("httr")
library("biomaRt")
library("rentrez")

setwd("~/Desktop/Single_cell_RNA-seq/GSE161277_RAW")

matrix_files <- list.files(pattern = "Patient[0-9]+_para-cancer_?\\d*_matrix\\.mtx\\.gz")
print(matrix_files)

# Initialize list to store Seurat objects
seurat_objects_05 <- list()

# Loop through each matrix file and create/define directory pot each patient
for(matrix_file in matrix_files) {
  base_name <- gsub("_matrix\\.mtx\\.gz", "", matrix_file)
  patient_dir <- paste0(getwd(), "/", base_name)
  
  if(!dir.exists(patient_dir)) {
    dir.create(patient_dir)
  }
  
# Move the files into the directory (Easy to organize)
  file.rename(from = matrix_file, to = paste0(patient_dir, "/matrix.mtx.gz"))
  file.rename(from = paste0(base_name, "_features.tsv.gz"), to = paste0(patient_dir, "/features.tsv.gz"))
  file.rename(from = paste0(base_name, "_barcodes.tsv.gz"), to = paste0(patient_dir, "/barcodes.tsv.gz"))
  
# Load the data using Seurat
  counts <- Read10X(data.dir = patient_dir)
  
# Create Seurat object for each dataset
  seurat_objects_05[[base_name]] <- CreateSeuratObject(counts = counts)
  seurat_objects_05[[base_name]]@meta.data$orig.ident <- base_name
}

########

all_seurat_objects <- c(seurat_objects_01, seurat_objects_02, seurat_objects_03, seurat_objects_04, seurat_objects_05)
ids <- paste("dataset", 1:length(all_seurat_objects), sep = "_") # To identify

# Merging step
merged_seurat <- merge(x = all_seurat_objects[[1]], y = all_seurat_objects[-1], add.cell.ids = ids)
merged_seurat
View(merged_seurat@meta.data)

# Separation for further identification of type or patients
merged_seurat$sample <- rownames(merged_seurat@meta.data)
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'orig.ident', into = c('GSM', 'Patient', 'Type', 'Num'), sep = '_', fill = 'right', remove = FALSE)

# Remove columns by their names
merged_seurat@meta.data <- merged_seurat@meta.data %>%
  select(-GSM, -sample)
View(merged_seurat@meta.data)

# Save the merged Seurat object to disk
saveRDS(merged_seurat, file = "~/Desktop/Single_cell_RNA-seq/merged_seurat.rds")

######### MERGED FINE #######
#### NOW TO PERFORM QUALITY CONTROL (QC) #####

# % MT reads
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
View(merged_seurat@meta.data)

merged_seurat@meta.data$sample_type <- merged_seurat@meta.data$Type

# Function to plot and save violin plots for each sample type
plot_and_save <- function(ident, filtered_seurat_obj) {
  subset_seurat <- subset(filtered_seurat_obj, Type == ident)

  p1 <- VlnPlot(subset_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Type")
  p2 <- FeatureScatter(subset_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = 'lm') + ggtitle(paste("Scatter Plot for", ident)) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5))
  
  pdf(paste0(ident, "_filtered_plots.pdf"))
  print(p1)
  plot.new()
  print(p2)
  dev.off()
}

# Extract unique types from the data
unique_identifiers <- unique(filtered_seurat_obj@meta.data$Type)
lapply(unique_identifiers, plot_and_save, filtered_seurat_obj)

#### PERFORM FILTERING - THIS IS ONLY FOR CALCULATION OF OUTLIERS #########

calculate_outliers <- function(type, merged_seurat) {
  subset_seurat <- subset(merged_seurat, subset = Type == type)
  
  feature_lower_threshold <- quantile(subset_seurat$nFeature_RNA, 0.05)
  feature_upper_threshold <- quantile(subset_seurat$nFeature_RNA, 0.95)
  count_lower_threshold <- quantile(subset_seurat$nCount_RNA, 0.05)
  count_upper_threshold <- quantile(subset_seurat$nCount_RNA, 0.95)

}

# This is only for calling the variables of the type again..
lapply(unique_identifiers, calculate_outliers, merged_seurat)

###### APPLY THE FILTERING AND SEE WHAT HAPPENS WITH THE SAMPLES ######

adenoma_filter <- subset(merged_seurat, Type == "adenoma" & 
                           nFeature_RNA >= 110 & nFeature_RNA <= 6000 & 
                           nCount_RNA <= 45000 & percent.mt <= 20)

carcinoma_filter <- subset(merged_seurat, Type == "carcinoma" & 
                             nFeature_RNA >= 120 & nFeature_RNA <= 6500 & 
                             nCount_RNA <= 50000 & percent.mt <= 20)

normal_filter <- subset(merged_seurat, Type == "normal" & 
                          nFeature_RNA >= 110 & nFeature_RNA <= 5000 & 
                          nCount_RNA <= 30000 & percent.mt <= 15)

blood_filter <- subset(merged_seurat, Type == "blood" & 
                         nFeature_RNA >= 400 & nFeature_RNA <= 2000 & 
                         nCount_RNA >= 500 & nCount_RNA <= 8000 & 
                         percent.mt <= 8)

paracancer_filter <- subset(merged_seurat, Type == "para-cancer" & 
                              nFeature_RNA >= 70 & nFeature_RNA <= 4500 & 
                              nCount_RNA >= 400 & nCount_RNA <= 25000 & percent.mt <= 25)

filtered_seurat_obj <- merge(adenoma_filter, y = c(carcinoma_filter, normal_filter, blood_filter, paracancer_filter))

### PLOT THEM WITH THE SAME FUNCTION THAT I USED PREVIOUSLY FOR THE VIOLIN AND DISPERSION ####

# To read it..
saveRDS(filtered_seurat_obj, file = "~/Desktop/Single_cell_RNA-seq/filtered_seurat_obj.rds")
filtered_seurat_obj <- readRDS("~/Desktop/Single_cell_RNA-seq/filtered_seurat_obj.rds")

####### Normalization, High variable features finding, data scalation, PCA ########

Normalized_seurat <- NormalizeData(object = filtered_seurat_obj)
Normalized_seurat <- FindVariableFeatures(object = Normalized_seurat)
Normalized_seurat <- ScaleData(object = Normalized_seurat)
Normalized_seurat <- RunPCA(object = Normalized_seurat)
ElbowPlot(object = Normalized_seurat)
Normalized_seurat <- FindNeighbors(object = Normalized_seurat, dims = 1:20)
Normalized_seurat <- FindClusters(object = Normalized_seurat)
Normalized_seurat <- RunUMAP(object = Normalized_seurat, dims = 1:20)

p1 <- DimPlot(Normalized_seurat, reduction = 'umap', group.by = 'Patient')
p2 <- DimPlot(Normalized_seurat, reduction = 'umap', group.by = 'Type', cols = c('yellow', 'green', 'blue', 'pink', 'purple'))
grid.arrange(p1, p2, ncol = 2, nrow = 2)

# Save the Processed Seurat object to disk
saveRDS(Normalized_seurat, file = "~/Desktop/Single_cell_RNA-seq/Normalized_seurat.rds")
Normalized_seurat <- readRDS("~/Desktop/Single_cell_RNA-seq/Normalized_seurat.rds")
class(Normalized_seurat[["RNA"]])  
Normalized_seurat

####### I HAVE SEURAT V5 ##### I WILL DO Anchor based integration (CCA) with IntegrateLayers instead of IntegrateData

############# CCA integration #############
obj <- IntegrateLayers(object = Normalized_seurat, method = CCAIntegration,
                            orig.reduction = "pca", new.reduction = "integrated.cca",
                            verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:20)
obj <- FindClusters(obj, resolution = 0.8, cluster.name = "cca_clusters")
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")

saveRDS(obj, file = "~/Desktop/Single_cell_RNA-seq/obj.rds")
obj <- readRDS("~/Desktop/Single_cell_RNA-seq/obj.rds")
obj

p3 <- DimPlot(obj, reduction = "umap.cca", group.by = "Patient") + ggtitle("cca_integrated_patient")
p4 <- DimPlot(obj, reduction = "umap.cca", group.by = "Type") + ggtitle("cca_integrated_type")
p3 + p4

obj <- JoinLayers(obj)

# For CCA integrated data
cca_markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save cca_markers as RDS
saveRDS(cca_markers, file = "~/Desktop/Single_cell_RNA-seq/cca_markers.rds")
cca_markers <- readRDS("~/Desktop/Single_cell_RNA-seq/cca_markers.rds")

# Group and visualize Top Markers for Each cluster
top10_cca <- cca_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
print(top10_cca)

# Visualize marker expression for cluster 0
FeaturePlot(obj, features = top10_cca$gene[top10_cca$cluster == 0], reduction = "umap.cca")

# Known markers from the paper
known_markers <- c(
  "CCL5" = "CD8+ T cells",
  "CD8A" = "CD8+ T cells",
  "CD8B" = "CD8+ T cells",
  "KLRC2" = "Natural Killer (NK) cells",
  "CXCR6" = "Follicular helper T cells",
  "CLDN2" = "Tight junctions, potentially involved in carcinogenesis",
  "MMP7" = "Matrix remodeling, cancer progression",
  "IL7R" = "Central memory T cells",
  "LEFTY1" = "Transition from normal to malignant colorectal tissues",
  "PCCA" = "Transition from normal to malignant colorectal tissues",
  "LGR5" = "Intestinal stem cells",
  "APCDD1" = "Physiological and benign epithelial cells",
  "SMOC2" = "Benign epithelial cells",
  "REG1A" = "Benign epithelial cells",
  "REG4" = "Carcinogenesis, carcinoma precursor cells",
  "C1QA" = "Immune response, complement system",
  "C1QB" = "Complement system, immune response",
  "C1QC" = "Complement system, immune response",
  "GUCA2A" = "Enterocytes",
  "GUCA2B" = "Enterocytes",
  "SLC26A3" = "Enterocytes",
  "BMX" = "Adenoma precursor cells, carcinogenesis",
  "SH2D6" = "Adenoma precursor cells, carcinogenesis",
  "MUC2" = "Goblet cells"
)

# Function to query NCBI for gene annotations with error handling
get_gene_annotation_ncbi <- function(gene_name) {
  print(paste("Fetching annotation for:", gene_name))
  tryCatch({
    result <- entrez_search(db = "gene", term = paste(gene_name, "[Symbol] AND Homo sapiens[Organism]"))
    if (result$count > 0) {
      gene_id <- result$ids[1]
      summary <- entrez_summary(db = "gene", id = gene_id)
      return(summary$description)
    }
    return("Unknown")
  }, error = function(e) {
    message(paste("Error fetching annotation for", gene_name, ":", e$message))
    return("Unknown")
  })
}

# Annotate markers
annotate_markers <- function(gene_list, known_markers) {
  annotations <- sapply(gene_list, function(gene) {
    if (gene %in% names(known_markers)) {
      return(known_markers[gene])
    } else {
      return(get_gene_annotation_ncbi(gene))
    }
  })
  return(annotations)
}

### Top 10 markers previously done #####
top10_cca$Annotation <- annotate_markers(top10_cca$gene, known_markers)
print(top10_cca)
saveRDS(top10_cca, file = "~/Desktop/Single_cell_RNA-seq/top10_cca_annotated.rds")
top10_cca <- readRDS("~/Desktop/Single_cell_RNA-seq/top10_cca_annotated.rds")


##### I already have the annotation for each cluster - cell type, What does it mean?, each cluster belongs to which type of tissue? ####

# Seurat_clusters as factors and character for mapping
obj$seurat_clusters <- as.character(obj$seurat_clusters)

# Vector to map clusters to cell types based on known markers, I did it manually..
cluster_to_celltype <- c(
  "0" = "CD8+ T cells",
  "1" = "B cells",
  "2" = "Epithelial cells",
  "3" = "Central memory - T cells",
  "4" = "T cells and NK cells",
  "5" = "Epithelial cells",
  "6" = "Central memory - T cells",
  "7" = "B cells",
  "8" = "Myeloid cells (monocytes or neutrophils)",
  "9" = "B cells",
  "10" = "NK cells",
  "11" = "Epithelial cells",
  "12" = "NK cells",
  "13" = "Epithelial cells",
  "14" = "Proliferating cells (likely cycling cells or stem cells)",
  "15" = "Epithelial cells",
  "16" = "Proliferating cells (likely cycling cells or stem cells)",
  "17" = "Myeloid cells (monocytes or neutrophils)",
  "18" = "Enterocytes",
  "19" = "Epithelial cells",
  "20" = "Endothelial cells",
  "21" = "Adenoma precursor cells",
  "22" = "Goblet cells",
  "23" = "Stromal cells"
)

# Map clusters to cell types
cell_types <- cluster_to_celltype[obj$seurat_clusters]
cell_types_df <- data.frame(cell_type = cell_types, row.names = colnames(obj))

# Add cell type annotations to the Seurat object metadata
obj <- AddMetaData(obj, metadata = cell_types_df)
saveRDS(obj, file = "~/Desktop/Single_cell_RNA-seq/annotated_obj.rds")

# Visualize UMAP with annotated clusters
p1 <- DimPlot(obj, reduction = "umap.cca", group.by = "cell_type") + ggtitle("UMAP with Cluster Annotations by Cell Type")
p2 <- DimPlot(obj, reduction = "umap.cca", group.by = "seurat_clusters") + ggtitle("UMAP with Cluster Annotations by Seurat Clusters")

p1+p4

