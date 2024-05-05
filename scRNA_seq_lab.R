# install packages from bioconductor
BiocManager::install("AnnotationHub")
BiocManager::install("ensembldb")
BiocManager::install("multtest")
BiocManager::install("glmGamPoi")


# install other packages from CRAN
install.packages("Matrix")
install.packages("RCurl") 
install.packages("scales")
install.packages("cowplot") 
install.packages("BiocManager")
install.packages("Seurat")
install.packages("metap")


# load libraries
library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(cowplot)
library(metap)
library(AnnotationHub)
library(ensembldb)
library(multtest)
library(glmGamPoi)


# Use for loop to create seurat objects for both samples
for (file in c("ctrl_raw_feature_bc_matrix", "stim_raw_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}



# examine metadata
head(ctrl_raw_feature_bc_matrix@meta.data)
head(stim_raw_feature_bc_matrix@meta.data)


# merge Seurat objects
merged_seurat <- merge(x = ctrl_raw_feature_bc_matrix,
                       y = stim_raw_feature_bc_matrix,
                       add.cell.id = c("ctrl", "stim"))



# examine metadata of merged object
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)



# open the metadata in a new tab
View(merged_seurat@meta.data)



# novelty score (# of genes per UMI for each cell)
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)



# add the percent of transcripts mapping to mitochondrial genes
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100



# make a dataframe of metadata
metadata <- merged_seurat@meta.data



# Add cell IDs to metadata
metadata$cells <- rownames(metadata)



# Add sample type column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells, "^stim_"))] <- "stim"



# rename the columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)



# add metadata back to seurat object and save as an RData object
merged_seurat@meta.data <- metadata
save(merged_seurat, file = "data/merged_filtered_seurat.RData")



# create bar chart of counts/sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


# Plot UMIs/transcripts per cell as a density plot
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)


# Use density plot to show distribution of genes detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)


# Use a density plot to show overall complexity of gene expression (novelty score)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


# Density plot of percent of mitochondrial genes
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)



# Display correlation between genes detected and UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)



# Filter low quality cells
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))




# Extract counts
filtered_seurat <- JoinLayers(filtered_seurat)
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output logical matrix showing which genes have 0 count
nonzero <- counts > 0



# Keep only genes expressed in 10 or more cells
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data

# recheck counts/sample using the filtered data
metadata_clean %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")



# Plot UMIs/transcripts per cell as a density plot using cleaned metadata
metadata_clean %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)



# Use density plot to show distribution of genes detected per cell
metadata_clean %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)


# Use a density plot to show overall complexity of gene expression (novelty score) using filtered data
metadata_clean %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


# Density plot of percent of mitochondrial genes
metadata_clean %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)


# Display correlation between genes detected and UMIs in the filtered dataset
metadata_clean %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~sample)


# Create .RData object
save(filtered_seurat, file="data/seurat_filtered.RData")


# Normalize counts
seurat_phase <- NormalizeData(filtered_seurat)



# Load cell cycle markers
load("data/cycle.rda")



# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells        
View(seurat_phase@meta.data)



# Identify most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale counts
seurat_phase <- ScaleData(seurat_phase)


# Identify 15 most highly variable genes
ranked_variable_genes <- VariableFeatures(seurat_phase)
top_genes <- ranked_variable_genes[1:15]

# Plot the avg expression and variance of these genes
p <- VariableFeaturePlot(seurat_phase)
LabelPoints(plot = p, points = top_genes, repel = TRUE)


# Run PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")




# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                   breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                   labels=c("Low","Medium","Medium high", "High"))
			


# Plot PCA grouped by amount of mitochondrial genes identified
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr")



############## Part 3 ############

## Integration

  
# Read in split_seurat
split_seurat <- readRDS("split_seurat.rds")
 

  
# Select highly variable features
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)
 

  
#prepare object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = integ_features)
 
  
# Find integration anchors
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
 

  
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
 

  
# Run and plot PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

PCAPlot(seurat_integrated,
        split.by = "sample")
 



  
# Set specific seed to ensure consistent umap
set.seed(123456)

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
			     reduction = "pca")

# Plot UMAP                             
DimPlot(seurat_integrated)    
 



  
# Plot UMAP split by sample
DimPlot(seurat_integrated,
        split.by = "sample")  
 



  
# Save the integrated dataset
saveRDS(seurat_integrated, "results/integrated_seurat.rds")
 

####### Clustering ######

  
# Create heatmap of integrated data

DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)
 



  
# Print out most variable genes driving each PC
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)
 
  
# Plot the elbow plot
ElbowPlot(object = seurat_integrated, 
          ndims = 40)
 


  
# Make k-nearest neighbot graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                dims = 1:40)
 

  
# Find clusters at different resolutions
seurat_integrated <- FindClusters(object = seurat_integrated,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
 
  
# Explore resolutions
seurat_integrated@meta.data %>% 
        View()
 

  
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

 

  
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
 

  
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
 


  
# Replot UMAP with resolution 0.8
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
 

###### Clustering quality control ######

  
# Determine no. of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "sample")) %>%
        dplyr::count(ident, sample)

# Barplot of number of cells per cluster by sample
ggplot(n_cells, aes(x=ident, y=n, fill=sample)) +
    geom_bar(position=position_dodge(), stat="identity") +
    geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))
 

  
# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()
 


  
# Barplot of proportion of cells in each cluster by sample
ggplot(seurat_integrated@meta.data) +
    geom_bar(aes(x=integrated_snn_res.0.8, fill=sample), position=position_fill()) 
 

  
# Use UMAP to see if clusters segregate by cell cycle phase
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()
 


  
# Determine metrics to plot present 
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
 


  
# Boxplot of nGene per cluster
ggplot(seurat_integrated@meta.data) +
    geom_boxplot(aes(x=integrated_snn_res.0.8, y=nGene, fill=integrated_snn_res.0.8)) +
    NoLegend()
 


  
# get data for principal components

columns <- c(paste0("PC_", 1:16), "ident", "umap_1", "umap_2")

# Extracting data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)
 

  
# Examine top 16 PCS

umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "umap_1", "umap_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(umap_1), y=mean(umap_2))
  
# Plot UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
        ggplot(pc_data, 
               aes(umap_1, umap_2)) +
                geom_point(aes_string(color=pc), 
                           alpha = 0.7) +
                scale_color_gradient(low = "grey90", 
                                     high = "blue")  +
                geom_text(data=umap_label, 
                          aes(label=ident, x, y)) +
                ggtitle(pc)
}) %>% 
        plot_grid(plotlist = .)
 

  
# Examine results of first 5 PCs
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)
 

  
# View the UMAP clustering again 

DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE) + NoLegend()
 

  
# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
 

  
# Plot umap with for CD14+ monocyte markers
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD14", "LYZ"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
 


  
# UMAP for FCGR3A+ monocyte markers
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("FCGR3A", "MS4A7"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
 

  
# Use UMAP to examine macrophage expression

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("MARCO", "ITGAM", "ADGRE1"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
 




  
# Use UMAP to examine conventional dendritic cell markers

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("FCER1A", "CST3"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
 


  
# Examine expression of plasmacytoid dendritic cell markers

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
 


  
# Create RNA expression dot plot
markers <- list()
markers[["CD14+ monocytes"]] <- c("CD14", "LYZ")
markers[["FCGR3A+ monocyte"]] <- c("FCGR3A", "MS4A7")
markers[["Macrophages"]] <- c("MARCO", "ITGAM", "ADGRE1")
markers[["Conventional dendritic"]] <- c("FCER1A", "CST3")
markers[["Plasmacytoid dendritic"]] <- c("IL3RA", "GZMB", "SERPINF1", "ITM2C")


DotPlot(seurat_integrated, markers, assay="RNA")
 


  
# Create RNA expression dot plot 
markers <- list()
markers[["B Cells"]] <- c("CD79A", "MS4A1")
markers[["T cells"]] <- c("CD3D")
markers[["CD4+ T cells"]] <- c("IL7R", "CCR7")
markers[["CD8+ T cells"]] <- c("CD8A")
markers[["NK cells"]] <- c("GNLY", "NKG7")
markers[["Megakaryocytes"]] <- c("PPBP")
markers[["Erythrocytes"]] <- c("HBB", "HBA2")
#markers[["Unknown"]] <- c("IL3RA", "GZMB", "SERPINF1", "ITM2C")


DotPlot(seurat_integrated, markers, assay="RNA")
 




####### Marker Identification #######

  
# Test makrer finding command on one clusters

cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                              ident.1 = 0,
                     	      grouping.var = "sample",
                              only.pos = TRUE,
		              logfc.threshold = 0.25)
 

  
# Add gene annotations

annotations <- read.csv("data/annotation.csv")
 

  
# Combine markers with gene descriptions

cluster0_ann_markers <- cluster0_conserved_markers %>% 
                rownames_to_column(var="gene") %>% 
                left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

View(cluster0_ann_markers)
 

  
# Identify conserved markers for cluster 10
cluster10_conserved_markers <- FindConservedMarkers(seurat_integrated,
                              ident.1 = 10,
                     	      grouping.var = "sample",
                              only.pos = TRUE,
		              logfc.threshold = 0.25)
 

  
# Combine markers with gene descriptions

cluster10_ann_markers <- cluster10_conserved_markers %>% 
                rownames_to_column(var="gene") %>% 
                left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

View(cluster10_ann_markers)
 


  
# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
               by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
  }

 

  
# Run function across select clusters
conserved_markers <- map_dfr(c(4,0,6,2), get_conserved)
 

  
# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (ctrl_avg_log2FC + stim_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)

# Visualize top 10 markers per cluster
View(top10)
 

  
# Plot interesting marker gene expression for cluster 4
FeaturePlot(object = seurat_integrated, 
                        features = c("HSPH1", "HSPE1", "DNAJB1"),
                         order = TRUE,
                         min.cutoff = 'q10', 
                         label = TRUE,
			 repel = TRUE)
 

  
# Vln plot - cluster 4
VlnPlot(object = seurat_integrated, 
        features = c("HSPH1", "HSPE1", "DNAJB1"))
 

  

# Determine differentiating markers for CD4+ T cell
cd4_tcells <- FindMarkers(seurat_integrated,
                          ident.1 = 2,
                          ident.2 = c(0,4,6))                  

# Add gene symbols to the DE table
cd4_tcells <- cd4_tcells %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
             by = c("gene" = "gene_name"))

# Reorder columns and sort by padj      
cd4_tcells <- cd4_tcells[, c(1, 3:5,2,6:7)]

cd4_tcells <- cd4_tcells %>%
  dplyr::arrange(p_val_adj) 

# View data
View(cd4_tcells)
 

  
# Rename all identities
seurat_integrated <- RenameIdents(object = seurat_integrated, 
                               "0" = "Naive or memory CD4+ T cells",
                               "1" = "CD14+ monocytes",
                               "2" = "Activated T cells",
                               "3" = "CD14+ monocytes",
                               "4" = "Stressed cells / Unknown",
                               "5" = "CD8+ T cells",
                               "6" = "Naive or memory CD4+ T cells",
                               "7" = "B cells",
                               "8" = "NK cells",
                               "9" = "CD8+ T cells",
                               "10" = "FCGR3A+ monocytes",
                               "11" = "B cells",
                               "12" = "NK cells",
                               "13" = "B cells",
                               "14" = "Conventional dendritic cells",
                               "15" = "Megakaryocytes",
			       "16" = "Plasmacytoid dendritic cells")


# Plot the UMAP
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)
 

  
# Remove the stressed or dying cells
seurat_subset_labeled <- subset(seurat_integrated,
                               idents = "Stressed cells / Unknown", invert = TRUE)

# Re-visualize the clusters
DimPlot(object = seurat_subset_labeled, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
	repel = TRUE)
 

  
# Save final R object
write_rds(seurat_integrated,
          file = "results/seurat_labelled.rds")

# Create and save a text file with sessionInfo
sink("sessionInfo_scrnaseq_Feb2023.txt")
sessionInfo()
sink()
 

  
# Add celltype annotation as a column in meta.data 
seurat_subset_labeled$celltype <- Idents(seurat_subset_labeled)

# Compute number of cells per celltype
n_cells <- FetchData(seurat_subset_labeled, 
                     vars = c("celltype", "sample")) %>%
        dplyr::count(celltype, sample)

# Barplot of number of cells per celltype by sample
ggplot(n_cells, aes(x=celltype, y=n, fill=sample)) +
    geom_bar(position=position_dodge(), stat="identity") +
    geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))
 

  
# Make sure labels are correct on seurat object
Idents(seurat_subset_labeled) <- seurat_subset_labeled$sample

# get subset of just bcells
seurat_b_cells <- subset(seurat_subset_labeled, subset = (celltype == "B cells"))

# Run a wilcox test to compare ctrl vs stim
b_markers <- FindMarkers(seurat_b_cells,
                                    ident.1 = "ctrl",
                                    ident.2 = "stim",
                                    grouping.var = "sample",
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
 
  
# Install enhanced volcano package
#BiocManager::install("EnhancedVolcano")
 


  
#Create volcano plot to show differentially expressed genes

library(EnhancedVolcano)
EnhancedVolcano(b_markers,
    row.names(b_markers),
    x="avg_log2FC",
    y="p_val_adj"
)
 
