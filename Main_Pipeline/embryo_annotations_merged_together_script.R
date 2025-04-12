# New L1 with P data and the SAA data all together


# Loading datasets - L1 data I generated first and SAA data
setwd("./Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/")

L1.n <- readRDS("./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/121923_L1_full_neuron_with_saa_doublets.rds")
L1.all <- readRDS("./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/062223_L1_all_cells_none_removed.cds")
new.saa <- readRDS("~/Dropbox (VU Basic Sciences)/L1 Cengen manuscript/Claire data and presentations/no_intron_analysis/data/042324_SAA_all_cells_reannotated_cds.rds.RDS")

# Loading in P lineage data
P.all <- readRDS("~/Dropbox (VU Basic Sciences)/L1 Cengen manuscript/Claire data and presentations/no_intron_analysis/data/041824_P-cell_full_dataset_cds.rds.RDS")

library(monocle3)
library(dplyr)
library(ggplot2)
library(wbData)

library(SingleCellExperiment) 
?monocle3::new_cell_data_set

#------A: Converting and combining all 14 datasets:
# Define the file paths for all datasets

file_paths <- c(
  "~/Lab_Data/SoupX/Murray_b01/Murray_b01_SoupX_corrected.rds",
  "~/Lab_Data/SoupX/Murray_b02/Murray_b02_SoupX_corrected.rds",
  "~/Lab_Data/SoupX/Murray_r17/Murray_r17_SoupX_corrected.rds",
  "~/Lab_Data/SoupX/Waterson_300_lane1/Waterson_300_lane1_sce_SoupX_corrected.rds",
  "~/Lab_Data/SoupX/Waterson_300_lane2/Waterson_300_lane2_sce_SoupX_corrected.rds",
  "~/Lab_Data/SoupX/Waterson_300_lane3/Waterson_300_lane3_sce_SoupX_corrected.rds",
  "~/Lab_Data/SoupX/Waterson_400_lane1/Waterson_400_lane1_sce_SoupX_corrected.rds",
  "~/Lab_Data/SoupX/Waterson_400_lane2/Waterson_400_lane2_sce_SoupX_corrected.rds",
  "~/Lab_Data/SoupX/Waterson_400_lane3/Waterson_400_lane3_sce_SoupX_corrected.rds",
  "~/Lab_Data/SoupX/Waterson_500_lane1_batch1/Waterson_500_lane1_batch1_sce_SoupX_corrected.rds",
  "~/Lab_Data/SoupX/Waterson_500_lane2_batch1/Waterson_500_lane2_batch1_sce_SoupX_corrected.rds",
  "~/Lab_Data/SoupX/Waterson_500_lane3_batch1/Waterson_500_lane3_batch1_sce_SoupX_corrected.rds",
  "~/Lab_Data/SoupX/Waterson_500_lane1_batch2/Waterson_500_lane1_batch2_sce_SoupX_corrected.rds",
  "~/Lab_Data/SoupX/Waterson_500_lane2_batch2/Waterson_500_lane2_batch2_sce_SoupX_corrected.rds"
)

sce <- readRDS("~/Lab_Data/SoupX/Murray_b01/Murray_b01_SoupX_corrected.rds")
assayNames(sce)

cds_list <- lapply(file_paths, function(path) {
  tryCatch({
    sce <- readRDS(path)
    
    #Changing col barcodes to given barcodes
    colData(sce)$matching_barcodes
    colnames(sce)
    colnames(sce) <- colData(sce)$matching_barcodes
    
    # Getting raw counts matrix(?)
    expression_matrix <- assay(sce, "counts")  #logcounts" for log-transformed data
    
    # Convert to CellDataSet, append
    cds <- new_cell_data_set(
      expression_data = expression_matrix,
      cell_metadata = colData(sce),
      gene_metadata = rowData(sce)
    )
    cds
  }, error = function(e) {
    message(sprintf("Something broke"))
    NULL  #Failsafe
  })
})

#colnames(sce_corrected) <- colData(sce_corrected)$matching_barcodes
# Filter out NULLs
cds_list <- cds_list[!sapply(cds_list, is.null)]
cds_list

#Single set from cds_list
combined_cds <- combine_cds(cds_list, keep_all_genes = TRUE, cell_names_unique = TRUE)
combined_cds
saveRDS(combined_cds, "~/Lab_Data/SoupX/combined_cds.rds")

print(combined_cds)

#-------A: End converting and combining

gids <- wb_load_gene_ids("WS273")

source("./Monocle3Functions.txt")

b01 = readRDS("~/Lab_Data/SoupX/Murray_b01/Murray_b01_SoupX_corrected.rds")
b02 = readRDS("~/Lab_Data/SoupX/Murray_b02/Murray_b02_SoupX_corrected.rds")
r17 = readRDS("~/Lab_Data/SoupX/Murray_r17/Murray_r17_SoupX_corrected.rds")
W300_lane1 = readRDS("~/Lab_Data/SoupX/Waterson_300_lane1/Waterson_300_lane1_sce_SoupX_corrected.rds")
W300_lane2 = readRDS("~/Lab_Data/SoupX/Waterson_300_lane2/Waterson_300_lane2_sce_SoupX_corrected.rds")
W300_lane3 = readRDS("~/Lab_Data/SoupX/Waterson_300_lane3/Waterson_300_lane3_sce_SoupX_corrected.rds")
W400_lane1 = readRDS("~/Lab_Data/SoupX/Waterson_400_lane1/Waterson_400_lane1_sce_SoupX_corrected.rds")
W400_lane2 = readRDS("~/Lab_Data/SoupX/Waterson_400_lane2/Waterson_400_lane2_sce_SoupX_corrected.rds")
W400_lane3 = readRDS("~/Lab_Data/SoupX/Waterson_400_lane3/Waterson_400_lane3_sce_SoupX_corrected.rds")
W500_lane1_b1 = readRDS("~/Lab_Data/SoupX/Waterson_500_lane1_batch1/Waterson_500_lane1_batch1_sce_SoupX_corrected.rds")
W500_lane2_b1 = readRDS("~/Lab_Data/SoupX/Waterson_500_lane2_batch1/Waterson_500_lane2_batch1_sce_SoupX_corrected.rds")
W500_lane3_b1 = readRDS("~/Lab_Data/SoupX/Waterson_500_lane3_batch1/Waterson_500_lane3_batch1_sce_SoupX_corrected.rds")
W500_lane1_b2 = readRDS("~/Lab_Data/SoupX/Waterson_500_lane1_batch2/Waterson_500_lane1_batch2_sce_SoupX_corrected.rds")
W500_lane2_b2 = readRDS("~/Lab_Data/SoupX/Waterson_500_lane2_batch2/Waterson_500_lane2_batch2_sce_SoupX_corrected.rds")
typeof(r17)
combined_cds <- combine_cds(list(b01, b02, r17, W300_lane1, W300_lane2, W300_lane3, W400_lane1, W400_lane2, W400_lane3, W500_lane1_b1, W500_lane2_b1, W500_lane3_b1, W500_lane1_b2, W500_lane2_b2), keep_all_genes = T, cell_names_unique = T)

table(colData(combined_cds)$Experiment, exclude = NULL)

# ---- Cont. A.A
combined_cds <- detect_genes(combined_cds)
combined_cds <- combined_cds[rowData(combined_cds)$num_cells_expressed > 5,]
combined_cds
# features in # cells
# 20191  143213 

colnames(colData(combined_cds))
table(combined_cds$in_processed)

combined_cds <- preprocess_cds(combined_cds, num_dim = 50) #A:A This 50 value should be where the elbow plot flattens. Should be in the range of 50-75
plot_pc_variance_explained(combined_cds)

# combined_cds <- align_cds(combined_cds, aligment_group = "batch", alignment_k = 5) #Correcting for batch differences (waterson vs murray)
# K value must be less than 14 (number of batches). Can adjust this (iteration, guess and check) to see if we get better results.
# If we see a cluster of many muscle cells, all in the same spot which are being treated differently, we can adjust this.

combined_cds <- reduce_dimension(combined_cds, #Takes a long time
                              reduction_method = "UMAP",
                              preprocess_method = "PCA",
                              umap.min_dist = 0.3,
                              umap.n_neighbors = 75)

combined_cds <- cluster_cells(combined_cds, reduction_method = "UMAP", res = 3e-5) #Res changes how many clusters there are. Default: 3e-5

plot_cells(combined_cds, 
           color_cells_by = "cell.subtype", #Change to whatever column we use for cell_type_data(cell.type or cell.subtype)
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(combined_cds, 
           color_cells_by = "cluster", # Plotting graph by cluster, with 31 clusters by default
           group_label_size = 4,
           cell_size = 0.5)

# sbt-1 and egl-21 are good markers of neuron cells that have finished dividing
# cdk-1 and mcm-4 are markers for actively diving cells (ANY cells)
plot_cells(combined_cds, genes = c("sbt-1","egl-21","cdk-1","mcm-4"), cell_size = 0.5)

colData(combined_cds)$cluster <- monocle3::clusters(combined_cds)

#plot_cells(combined_cds,  #Unused
#           color_cells_by = "hph",
#           label_cell_groups = F,
#           cell_size = 0.5)

colData(combined_cds)$UMAP_1 <- reducedDims(combined_cds)[["UMAP"]][,1]
colData(combined_cds)$UMAP_2 <- reducedDims(combined_cds)[["UMAP"]][,2]
plot.expr.UMAP(combined_cds, "sbt-1", size = 0.5)

#saveRDS(combined_cds, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/042324_L1_all_live_cells_MeOH_new_saa_x2_P_lineage_with_doublets_UPR_cds_no_intronic_reads.rds")

neuron <- combined_cds[,colData(combined_cds)$cluster %in% c(1, 2, 3, 19, 29, 26, 24)] #Taking these clusters from the cluster plot
# Can repeat this inside of neuron in order to get more information (i.e just subset certain values in the top right)

neuron <- detect_genes(neuron)
neuron <- neuron[rowData(neuron)$num_cells_expressed > 5,]
neuron
# 16504 features in 38250 cells

neuron <- preprocess_cds(neuron, num_dim = 50) 
plot_pc_variance_explained(neuron)

# neuron <- align_cds(neuron, aligment_group = "Experiment", alignment_k = 5)

neuron <- reduce_dimension(neuron, #Takes awhile
                              reduction_method = "UMAP",
                              preprocess_method = "PCA",
                              umap.min_dist = 0.3,
                              umap.n_neighbors = 75)
                              #For 3D, need to set a new parameter(max_componets?) Would need a new plotting function

neuron <- cluster_cells(neuron, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

#Plot all annot\ated cell types
sort(unique(colData(neuron)$cell.subtype))

plot_cells(neuron, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

colnames(colData(neuron))
table(neuron$in_processed)

#plot_cells(neuron, 
#           color_cells_by = "hph",
#           label_cell_groups = F,
#           cell_size = 0.5)

# Make sure to run this line vvvv
colData(neuron)$cluster <- monocle3::clusters(neuron)

colData(neuron)$PCA.partition <- monocle3::partitions(neuron, reduction_method = "PCA")

plot_cells(neuron,
           color_cells_by = "PCA.partition",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

plot_cells(neuron,
           color_cells_by = "PCA.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(neuron)$UMAP_1 <- reducedDims(neuron)[["UMAP"]][,1]
colData(neuron)$UMAP_2 <- reducedDims(neuron)[["UMAP"]][,2]

markers <- top_markers(neuron, group_cells_by = "cluster") #Takes ~12 min
markers %>% filter(cell_group == 119) %>% arrange(desc(marker_score)) %>% head(20) #126 total in neuron

rowData(neuron)$gene_short_name <- i2s(rowData(neuron)$id, gids)

plot.expr.UMAP(neuron, "sbt-1", size = 0.5)
plot.expr.UMAP(neuron, "cdk-1", size = 0.5) #For progenitors
plot.expr.UMAP(neuron, "mcm-4", size = 0.5) #For progenitors
plot.expr.UMAP(neuron, "nlp-56", size = 0.5)
plot.expr.UMAP(neuron, "gcy-33", size = 0.5)
plot.expr.UMAP(neuron, "ceh-36", size = 0.5)

#A:A Ended here, begin process of cellsubtyping and 
#1st Plot: Cluster 41 = URX
#neuron <- combined_cds[,colData(combined_cds)$cluster %in% c(1, 2, 3, 19, 29, 26, 24)] #Taking these clusters from the cluster plot

# I NEED to check if adding 36 fixes this.
neuron_2 <- neuron[, !(colData(neuron)$cluster %in% c(13, 20, 41, 51, 
                                                      57, 62, 63, 69, 
                                                      70, 73, 76, 85, 
                                                      87, 94, 101, 119, 
                                                      109, 117))]#Taking these clusters from the cluster plot
# Can repeat this inside of neuron in order to get more information (i.e just subset certain values in the top right)

neuron_2 <- detect_genes(neuron_2)
neuron_2 <- neuron_2[rowData(neuron_2)$num_cells_expressed > 5,]
neuron_2
# 15966 features in 33411  cells

neuron_2 <- preprocess_cds(neuron_2, num_dim = 50) 
plot_pc_variance_explained(neuron_2)

# neuron <- align_cds(neuron, aligment_group = "Experiment", alignment_k = 5)

neuron_2 <- reduce_dimension(neuron_2, #Takes awhile
                           reduction_method = "UMAP",
                           preprocess_method = "PCA",
                           umap.min_dist = 0.3,
                           umap.n_neighbors = 75)
#For 3D, need to set a new parameter(max_componets?) Would need a new plotting function

neuron_2 <- cluster_cells(neuron_2, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_2, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_2, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)
#This code needed to go here... vvv
colData(neuron_2)$cluster <- monocle3::clusters(neuron_2)

markers_2 <- top_markers(neuron_2, group_cells_by = "cluster") 
markers_2 %>% filter(cell_group == 99) %>% arrange(desc(marker_score)) %>% head(20) 

neuron_3 <- neuron_2[, !(colData(neuron_2)$cluster %in% c(32, 68, 70, 71, 77, 
                                                          81, 85, 90, 91, 93, 
                                                          97, 98, 100, 104, 105, 
                                                          106, 108, 109, 110))]#Taking these clusters from the cluster plot

neuron_3 <- detect_genes(neuron_3)
neuron_3 <- neuron_3[rowData(neuron_3)$num_cells_expressed > 5,]
neuron_3
# 15527 features in 30166  cells

neuron_3 <- preprocess_cds(neuron_3, num_dim = 50) 
plot_pc_variance_explained(neuron_3)

neuron_3 <- reduce_dimension(neuron_3, #Takes awhile
                             reduction_method = "UMAP",
                             preprocess_method = "PCA",
                             umap.min_dist = 0.3,
                             umap.n_neighbors = 75)

neuron_3 <- cluster_cells(neuron_3, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_3, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_3, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_3)$cluster <- monocle3::clusters(neuron_3)

markers_3 <- top_markers(neuron_3, group_cells_by = "cluster") 
markers_3 %>% filter(cell_group == 1) %>% arrange(desc(marker_score)) %>% head(20) %>% select(gene_short_name)

neuron_4 <- neuron_3[, !(colData(neuron_3)$cluster %in% c(5, 43, 80, 84, 91, 
                                                          93, 96, 98, 99, 102))]#Taking these clusters from the cluster plot
neuron_4 <- detect_genes(neuron_4)
neuron_4 <- neuron_4[rowData(neuron_4)$num_cells_expressed > 5,]
neuron_4
# 15105 in 28255 cells

neuron_4 <- preprocess_cds(neuron_4, num_dim = 50) 
plot_pc_variance_explained(neuron_4)

neuron_4 <- reduce_dimension(neuron_4, #Takes awhile
                             reduction_method = "UMAP",
                             preprocess_method = "PCA",
                             umap.min_dist = 0.3,
                             umap.n_neighbors = 75)

neuron_4 <- cluster_cells(neuron_4, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_4, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_4, 
             color_cells_by = "cluster",
             group_label_size = 4,
             cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_4)$cluster <- monocle3::clusters(neuron_4)

markers_4 <- top_markers(neuron_4, group_cells_by = "cluster") 
markers_4 %>% filter(cell_group == 88) %>% arrange(desc(marker_score)) %>% head(20) %>% select(gene_short_name)


neuron_5 <- neuron_4[, !(colData(neuron_4)$cluster %in% c(37, 46, 55, 67, 79, 
                                                          84, 88, 89, 90))]#Taking these clusters from the cluster plot

neuron_5 <- detect_genes(neuron_5)
neuron_5 <- neuron_5[rowData(neuron_5)$num_cells_expressed > 5,]
neuron_5
# 14893 features in 26592 cells

neuron_5 <- preprocess_cds(neuron_5, num_dim = 50) 
plot_pc_variance_explained(neuron_5)

neuron_5 <- reduce_dimension(neuron_5, #Takes awhile
                             reduction_method = "UMAP",
                             preprocess_method = "PCA",
                             umap.min_dist = 0.25,# moved to .25 from .3 (fascilitates more separation)
                             umap.n_neighbors = 65) # moved to 65 from 75 (fascilitates more separation)

neuron_5 <- cluster_cells(neuron_5, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_5, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_5, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_5)$cluster <- monocle3::clusters(neuron_5)

markers_5 <- top_markers(neuron_5, group_cells_by = "cluster") 
markers_5 %>% filter(cell_group == 9) %>% arrange(desc(marker_score)) %>% head(20) %>% select(gene_short_name)


neuron_6 <- neuron_5[, !(colData(neuron_5)$cluster %in% c(48, 50, 54, 71, 82, 
                                                          87, 89, 94, 96, 99))]#Taking these clusters from the cluster plot
neuron_6 <- detect_genes(neuron_6)
neuron_6 <- neuron_6[rowData(neuron_6)$num_cells_expressed > 5,]
neuron_6
# 14700 features in 25149 cells

neuron_6 <- preprocess_cds(neuron_6, num_dim = 50) 
plot_pc_variance_explained(neuron_6)

neuron_6 <- reduce_dimension(neuron_6, #Takes awhile
                             reduction_method = "UMAP",
                             preprocess_method = "PCA",
                             umap.min_dist = 0.25, # moved to .25 from .3 (fascilitates more separation)
                             umap.n_neighbors = 65)  # moved to 65 from 75 (fascilitates more separation)

neuron_6 <- cluster_cells(neuron_6, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_6, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_6, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_6)$cluster <- monocle3::clusters(neuron_6)

markers_6 <- top_markers(neuron_6, group_cells_by = "cluster")
markers_6 %>% filter(cell_group == 95) %>% arrange(desc(marker_score)) %>% head(20) %>% select(gene_short_name)

colData(neuron_6)$UMAP_1 <- reducedDims(neuron_6)[["UMAP"]][,1]
colData(neuron_6)$UMAP_2 <- reducedDims(neuron_6)[["UMAP"]][,2]

plot.expr.UMAP(neuron_6, "unc-54", size = 0.5)


neuron_7 <- neuron_6[, !(colData(neuron_6)$cluster %in% c(5, 22, 36, 55, 78, 
                                                          85, 86))]#Taking these clusters from the cluster plot
neuron_7 <- detect_genes(neuron_7)
neuron_7 <- neuron_7[rowData(neuron_7)$num_cells_expressed > 5,]
neuron_7
# 14440 features in 23401 cells

neuron_7 <- preprocess_cds(neuron_7, num_dim = 50) 
plot_pc_variance_explained(neuron_7)

neuron_7 <- reduce_dimension(neuron_7, #Takes awhile
                             reduction_method = "UMAP",
                             preprocess_method = "PCA",
                             umap.min_dist = 0.25, # moved to .25 from .3 (fascilitates more separation)
                             umap.n_neighbors = 65)  # moved to 65 from 75 (fascilitates more separation)

neuron_7 <- cluster_cells(neuron_7, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_7, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_7, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_7)$cluster <- monocle3::clusters(neuron_7)

markers_7 <- top_markers(neuron_7, group_cells_by = "cluster") 
markers_7 %>% filter(cell_group == 7) %>% arrange(desc(marker_score)) %>% head(20) %>% select(gene_short_name)

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_7)$UMAP_1 <- reducedDims(neuron_7)[["UMAP"]][,1]
colData(neuron_7)$UMAP_2 <- reducedDims(neuron_7)[["UMAP"]][,2]
plot.expr.UMAP(neuron_7, "egl-21", size = 0.5)


neuron_8 <- neuron_7[, !(colData(neuron_7)$cluster %in% c(11, 57, 61, 71, 76, 
                                                          82, 83, 84))]#Taking these clusters from the cluster plot
neuron_8 <- detect_genes(neuron_8)
neuron_8 <- neuron_8[rowData(neuron_8)$num_cells_expressed > 5,]
neuron_8
# 14282 features in 22116 cells

neuron_8 <- preprocess_cds(neuron_8, num_dim = 50) 
plot_pc_variance_explained(neuron_8)

neuron_8 <- reduce_dimension(neuron_8, #Takes awhile
                             reduction_method = "UMAP",
                             preprocess_method = "PCA",
                             umap.min_dist = 0.25, 
                             umap.n_neighbors = 65) 

neuron_8 <- cluster_cells(neuron_8, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_8, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_8, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_8)$cluster <- monocle3::clusters(neuron_8)

markers_8 <- top_markers(neuron_8, group_cells_by = "cluster") 
markers_8 %>% filter(cell_group == 48) %>% arrange(desc(marker_score)) %>% head(25) %>% select(gene_short_name)

# Extract metadata from neuron_8
metadata <- colData(neuron_8) %>% as.data.frame()
# Count the number of cells for each annotation in cluster 47
metadata %>%
  filter(cluster == 48) %>%
  count(cell.subtype, sort = TRUE)  # Sorts by highest count


# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_8)$UMAP_1 <- reducedDims(neuron_8)[["UMAP"]][,1]
colData(neuron_8)$UMAP_2 <- reducedDims(neuron_8)[["UMAP"]][,2]
plot.expr.UMAP(neuron_8, "egl-21", size = 0.5)
#actively dividing: cdk-1 and mcm-4
#finished neuronal: sbt-1 and egl-21


neuron_9 <- neuron_8[, !(colData(neuron_8)$cluster %in% c(38,57,64,66,68,71,73,
                                                          76,81,84))]#Taking these clusters from the cluster plot
neuron_9 <- detect_genes(neuron_9)
neuron_9 <- neuron_9[rowData(neuron_9)$num_cells_expressed > 5,]
neuron_9
# 14103 features in 20742 cells

neuron_9 <- preprocess_cds(neuron_9, num_dim = 50) 
plot_pc_variance_explained(neuron_9)

neuron_9 <- reduce_dimension(neuron_9, #Takes awhile
                             reduction_method = "UMAP",
                             preprocess_method = "PCA",
                             umap.min_dist = 0.25, 
                             umap.n_neighbors = 65)

neuron_9 <- cluster_cells(neuron_9, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_9, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_9, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)


#This code needed to go here... vvv
colData(neuron_9)$cluster <- monocle3::clusters(neuron_9)

markers_9 <- top_markers(neuron_9, group_cells_by = "cluster") 
markers_9 %>% filter(cell_group == 76) %>% arrange(desc(marker_score)) %>% head(25) %>% select(gene_short_name)

# Extract metadata from neuron_9
metadata <- colData(neuron_9) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 76) %>%
  count(cell.subtype, sort = TRUE)  # Sorts by highest count


# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_9)$UMAP_1 <- reducedDims(neuron_9)[["UMAP"]][,1]
colData(neuron_9)$UMAP_2 <- reducedDims(neuron_9)[["UMAP"]][,2]
plot.expr.UMAP(neuron_9, "F44E5.4", size = 0.5)
#actively dividing: cdk-1 and mcm-4
#finished neuronal: sbt-1 and egl-21


neuron_10 <- neuron_9[, !(colData(neuron_9)$cluster %in% c(20,44,50,61,65,71,
                                                           73,75,78,79))]#Taking these clusters from the cluster plot
neuron_10 <- detect_genes(neuron_10)
neuron_10 <- neuron_10[rowData(neuron_10)$num_cells_expressed > 5,]
neuron_10
# 13906 features in 19259 cells

neuron_10 <- preprocess_cds(neuron_10, num_dim = 50) 
plot_pc_variance_explained(neuron_10)

neuron_10 <- reduce_dimension(neuron_10, #Takes awhile
                             reduction_method = "UMAP",
                             preprocess_method = "PCA",
                             umap.min_dist = 0.225, # moved to .225 from .25 (fascilitates more separation)
                             umap.n_neighbors = 60)  # moved to 60 from 65 (fascilitates more separation)

neuron_10 <- cluster_cells(neuron_10, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_10, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_10, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_10)$cluster <- monocle3::clusters(neuron_10)

markers_10 <- top_markers(neuron_10, group_cells_by = "cluster") 
markers_10 %>% filter(cell_group == 48) %>% arrange(desc(marker_score)) %>% head(20) %>% select(gene_short_name)

metadata <- colData(neuron_10) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 48) %>%
  count(cell.subtype, sort = TRUE)  # Sorts by highest count

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_10)$UMAP_1 <- reducedDims(neuron_10)[["UMAP"]][,1]
colData(neuron_10)$UMAP_2 <- reducedDims(neuron_10)[["UMAP"]][,2]
plot.expr.UMAP(neuron_10, "egl-21", size = 0.5)



neuron_11 <- neuron_10[, !(colData(neuron_10)$cluster %in% c(5,20,43,52,62,65,
                                                             68,71,72,74))]#Taking these clusters from the cluster plot
neuron_11 <- detect_genes(neuron_11)
neuron_11 <- neuron_11[rowData(neuron_11)$num_cells_expressed > 5,]
neuron_11
# 13609 features in 17602 cells

neuron_11 <- preprocess_cds(neuron_11, num_dim = 50) 
plot_pc_variance_explained(neuron_11)

neuron_11 <- reduce_dimension(neuron_11, #Takes awhile
                              reduction_method = "UMAP",
                              preprocess_method = "PCA",
                              umap.min_dist = 0.225, 
                              umap.n_neighbors = 60)

neuron_11 <- cluster_cells(neuron_11, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_11, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_11, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_11)$cluster <- monocle3::clusters(neuron_11)

markers_11 <- top_markers(neuron_11, group_cells_by = "cluster") 
markers_11 %>% filter(cell_group == 39) %>% arrange(desc(marker_score)) %>% head(20) %>% select(gene_short_name)

metadata <- colData(neuron_11) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 39) %>%
  count(cell.subtype, sort = TRUE)  # Sorts by highest count

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_11)$UMAP_1 <- reducedDims(neuron_11)[["UMAP"]][,1]
colData(neuron_11)$UMAP_2 <- reducedDims(neuron_11)[["UMAP"]][,2]
plot.expr.UMAP(neuron_11, "egl-21", size = 0.5)



neuron_12 <- neuron_11[, !(colData(neuron_11)$cluster %in% c(12,18,39,40,46,47,
                                                             53,59,61,62,65,66,
                                                             68))]#Taking these clusters from the cluster plot
neuron_12 <- detect_genes(neuron_12)
neuron_12 <- neuron_12[rowData(neuron_12)$num_cells_expressed > 5,]
neuron_12
# 13299 features in 15398 cells

neuron_12 <- preprocess_cds(neuron_12, num_dim = 50) 
plot_pc_variance_explained(neuron_12)

neuron_12 <- reduce_dimension(neuron_12, #Takes awhile
                              reduction_method = "UMAP",
                              preprocess_method = "PCA",
                              umap.min_dist = 0.2, # moved to .2 from .225 (fascilitates more separation)
                              umap.n_neighbors = 25)  # moved to 50 from 60 (fascilitates more separation)

neuron_12 <- cluster_cells(neuron_12, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_12, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_12, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_12)$cluster <- monocle3::clusters(neuron_12)

markers_12 <- top_markers(neuron_12, group_cells_by = "cluster") 
markers_12 %>% filter(cell_group == 50) %>% arrange(desc(marker_score)) %>% head(20) %>% select(gene_short_name)

metadata <- colData(neuron_12) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 50) %>%
  count(cell.subtype, sort = TRUE)  # Sorts by highest count

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_12)$UMAP_1 <- reducedDims(neuron_12)[["UMAP"]][,1]
colData(neuron_12)$UMAP_2 <- reducedDims(neuron_12)[["UMAP"]][,2]
plot.expr.UMAP(neuron_12, "linc-104", size = 0.5)

colData(neuron_12)$cluster <- monocle3::clusters(neuron_12)


neuron_13 <- neuron_12[, !(colData(neuron_12)$cluster %in% c(31,38,42,44,45,49,
                                                             53,55,56,57,58,60,
                                                             62,64))]#Taking these clusters from the cluster plot
neuron_13 <- detect_genes(neuron_13)
neuron_13 <- neuron_13[rowData(neuron_13)$num_cells_expressed > 5,]
neuron_13
# 12846 features in 13574 cells

neuron_13 <- preprocess_cds(neuron_13, num_dim = 50) 
plot_pc_variance_explained(neuron_13)

neuron_13 <- reduce_dimension(neuron_13, #Takes awhile
                              reduction_method = "UMAP",
                              preprocess_method = "PCA",
                              umap.min_dist = 0.2, 
                              umap.n_neighbors = 50)  

neuron_13 <- cluster_cells(neuron_13, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_13, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_13, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_13)$cluster <- monocle3::clusters(neuron_13)

markers_13 <- top_markers(neuron_13, group_cells_by = "cluster") 
markers_13 %>% filter(cell_group == 15) %>% arrange(desc(marker_score)) %>% head(20)

metadata <- colData(neuron_13) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 15) %>%
  count(cell.subtype, sort = TRUE)  # Sorts by highest countcount(cell.subtype, sort = TRUE)  # Sorts by highest count

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_13)$UMAP_1 <- reducedDims(neuron_13)[["UMAP"]][,1]
colData(neuron_13)$UMAP_2 <- reducedDims(neuron_13)[["UMAP"]][,2]
plot.expr.UMAP(neuron_13, "egl-21", size = 0.5)


neuron_14 <- neuron_13[, !(colData(neuron_13)$cluster %in% c(15,22,35,41,48,53,56,57,58))]#Taking these clusters from the cluster plot
neuron_14 <- detect_genes(neuron_14)
neuron_14 <- neuron_14[rowData(neuron_14)$num_cells_expressed > 5,]
neuron_14
# 12696 features in 12427 cells

neuron_14 <- preprocess_cds(neuron_14, num_dim = 50) 
plot_pc_variance_explained(neuron_14)

neuron_14 <- reduce_dimension(neuron_14, #Takes awhile
                              reduction_method = "UMAP",
                              preprocess_method = "PCA",
                              umap.min_dist = 0.2, 
                              umap.n_neighbors = 50)  

neuron_14 <- cluster_cells(neuron_14, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_14, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_14, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_14)$cluster <- monocle3::clusters(neuron_14)

markers_13 <- top_markers(neuron_14, group_cells_by = "cluster") 
markers_13 %>% filter(cell_group == 13) %>% arrange(desc(marker_score)) %>% head(20)

metadata <- colData(neuron_14) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 13) %>%
  count(cell.subtype, sort = TRUE)  # Sorts by highest countcount(cell.subtype, sort = TRUE)  # Sorts by highest count

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_14)$UMAP_1 <- reducedDims(neuron_14)[["UMAP"]][,1]
colData(neuron_14)$UMAP_2 <- reducedDims(neuron_14)[["UMAP"]][,2]
plot.expr.UMAP(neuron_14, "egl-21", size = 0.5)



neuron_15 <- neuron_14[, !(colData(neuron_14)$cluster %in% c(4,13,18,21,29,37,40,
                                                             51,53))]#Taking these clusters from the cluster plot
neuron_15 <- detect_genes(neuron_15)
neuron_15 <- neuron_15[rowData(neuron_15)$num_cells_expressed > 5,]
neuron_15
# 12258 features in 10490 cells

neuron_15 <- preprocess_cds(neuron_15, num_dim = 50) 
plot_pc_variance_explained(neuron_15)

neuron_15 <- reduce_dimension(neuron_15, #Takes awhile
                              reduction_method = "UMAP",
                              preprocess_method = "PCA",
                              umap.min_dist = 0.2, 
                              umap.n_neighbors = 50)  

neuron_15 <- cluster_cells(neuron_15, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_15, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_15, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_15)$cluster <- monocle3::clusters(neuron_15)

markers_15 <- top_markers(neuron_15, group_cells_by = "cluster") 
markers_15 %>% filter(cell_group == 22) %>% arrange(desc(marker_score)) %>% head(20)

metadata <- colData(neuron_15) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 22) %>%
  count(cell.subtype, sort = TRUE)  # Sorts by highest countcount(cell.subtype, sort = TRUE)  # Sorts by highest count

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_15)$UMAP_1 <- reducedDims(neuron_15)[["UMAP"]][,1]
colData(neuron_15)$UMAP_2 <- reducedDims(neuron_15)[["UMAP"]][,2]
plot.expr.UMAP(neuron_15, "sbt-1", size = 0.5)



neuron_16 <- neuron_15[, !(colData(neuron_15)$cluster %in% c(21,22,28,29,31,34,
                                                             35,37,41,42,45,47,
                                                             48,49))]#Taking these clusters from the cluster plot
neuron_16 <- detect_genes(neuron_16)
neuron_16 <- neuron_16[rowData(neuron_16)$num_cells_expressed > 5,]
neuron_16
# 11895 features in 8928 cells

neuron_16 <- preprocess_cds(neuron_16, num_dim = 50) 
plot_pc_variance_explained(neuron_16)

neuron_16 <- reduce_dimension(neuron_16, #Takes awhile
                              reduction_method = "UMAP",
                              preprocess_method = "PCA",
                              umap.min_dist = 0.15, # moved to .15 from .2 (fascilitates more separation)
                              umap.n_neighbors = 40) # moved to 40 from 50 (fascilitates more separation)

neuron_16 <- cluster_cells(neuron_16, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_16, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_16, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_16)$cluster <- monocle3::clusters(neuron_16)

markers_16 <- top_markers(neuron_16, group_cells_by = "cluster") 
markers_16 %>% filter(cell_group == 1) %>% arrange(desc(marker_score)) %>% head(20)

metadata <- colData(neuron_16) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 1) %>%
  count(cell.subtype, sort = TRUE)  # Sorts by highest countcount(cell.subtype, sort = TRUE)  # Sorts by highest count

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_16)$UMAP_1 <- reducedDims(neuron_16)[["UMAP"]][,1]
colData(neuron_16)$UMAP_2 <- reducedDims(neuron_16)[["UMAP"]][,2]
plot.expr.UMAP(neuron_16, "sbt-1", size = 0.5)



neuron_17 <- neuron_16[, !(colData(neuron_16)$cluster %in% c(1,4,24,30,34,36,
                                                             37,38,39))]#Taking these clusters from the cluster plot
neuron_17 <- detect_genes(neuron_17)
neuron_17 <- neuron_17[rowData(neuron_17)$num_cells_expressed > 5,]
neuron_17
# 11264 features in 7198 cells

neuron_17 <- preprocess_cds(neuron_17, num_dim = 50) 
plot_pc_variance_explained(neuron_17)

neuron_17 <- reduce_dimension(neuron_17, #Takes awhile
                              reduction_method = "UMAP",
                              preprocess_method = "PCA",
                              umap.min_dist = 0.15, # moved to .15 from .2 (fascilitates more separation)
                              umap.n_neighbors = 40) # moved to 40 from 50 (fascilitates more separation)

neuron_17 <- cluster_cells(neuron_17, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_17, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_17, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_17)$cluster <- monocle3::clusters(neuron_17)

markers_17 <- top_markers(neuron_17, group_cells_by = "cluster") 
markers_17 %>% filter(cell_group == 3) %>% arrange(desc(marker_score)) %>% head(20)

metadata <- colData(neuron_17) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 3) %>%
  count(cell.subtype, sort = TRUE)  # Sorts by highest countcount(cell.subtype, sort = TRUE)  # Sorts by highest count

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_17)$UMAP_1 <- reducedDims(neuron_17)[["UMAP"]][,1]
colData(neuron_17)$UMAP_2 <- reducedDims(neuron_17)[["UMAP"]][,2]
plot.expr.UMAP(neuron_17, "F44E5.4", size = 0.5)



neuron_18 <- neuron_17[, !(colData(neuron_17)$cluster %in% c(3,6,17,21,23,24,
                                                             26,27,28,29,30))]#Taking these clusters from the cluster plot
neuron_18 <- detect_genes(neuron_18)
neuron_18 <- neuron_18[rowData(neuron_18)$num_cells_expressed > 5,]
neuron_18
# 10776 features in 5562 cells

neuron_18 <- preprocess_cds(neuron_18, num_dim = 50) 
plot_pc_variance_explained(neuron_18)

neuron_18 <- reduce_dimension(neuron_18, #Takes awhile
                              reduction_method = "UMAP",
                              preprocess_method = "PCA",
                              umap.min_dist = 0.1, # moved to .1 from .15 (fascilitates more separation)
                              umap.n_neighbors = 30) # moved to 30 from 40 (fascilitates more separation)

neuron_18 <- cluster_cells(neuron_18, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_18, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_18, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_18)$cluster <- monocle3::clusters(neuron_18)

markers_18 <- top_markers(neuron_18, group_cells_by = "cluster") 
markers_18 %>% filter(cell_group == 7) %>% arrange(desc(marker_score)) %>% head(20)

metadata <- colData(neuron_18) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 7) %>%
  count(cell.subtype, sort = TRUE)  # Sorts by highest countcount(cell.subtype, sort = TRUE)  # Sorts by highest count

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_18)$UMAP_1 <- reducedDims(neuron_18)[["UMAP"]][,1]
colData(neuron_18)$UMAP_2 <- reducedDims(neuron_18)[["UMAP"]][,2]
plot.expr.UMAP(neuron_18, "sbt-1", size = 0.5)



neuron_19 <- neuron_18[, !(colData(neuron_18)$cluster %in% c(1,7,11,13,16,18,
                                                             20,23,25,27,28,29,
                                                             30,31))]#Taking these clusters from the cluster plot
neuron_19 <- detect_genes(neuron_19)
neuron_19 <- neuron_19[rowData(neuron_19)$num_cells_expressed > 5,]
neuron_19
# 10008 features in 3651 cells

neuron_19 <- preprocess_cds(neuron_19, num_dim = 50) 
plot_pc_variance_explained(neuron_19)

neuron_19 <- reduce_dimension(neuron_19, #Takes awhile
                              reduction_method = "UMAP",
                              preprocess_method = "PCA",
                              umap.min_dist = 0.1, # moved to .1 from .15 (fascilitates more separation)
                              umap.n_neighbors = 30) # moved to 30 from 40 (fascilitates more separation)

neuron_19 <- cluster_cells(neuron_19, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_19, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_19, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_19)$cluster <- monocle3::clusters(neuron_19)

markers_19 <- top_markers(neuron_19, group_cells_by = "cluster") 
markers_19 %>% filter(cell_group == 1) %>% arrange(desc(marker_score)) %>% head(20)

metadata <- colData(neuron_19) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 1) %>%
  count(cell.subtype, sort = TRUE)  # Sorts by highest countcount(cell.subtype, sort = TRUE)  # Sorts by highest count

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_19)$UMAP_1 <- reducedDims(neuron_19)[["UMAP"]][,1]
colData(neuron_19)$UMAP_2 <- reducedDims(neuron_19)[["UMAP"]][,2]
plot.expr.UMAP(neuron_19, "sbt-1", size = 0.5)



neuron_20 <- neuron_19[, !(colData(neuron_19)$cluster %in% c(1,2,7,8,9,10,11,
                                                             12,13,14,15,16,17,
                                                             18))]#Taking these clusters from the cluster plot
neuron_20 <- detect_genes(neuron_20)
neuron_20 <- neuron_20[rowData(neuron_20)$num_cells_expressed > 5,]
neuron_20
# 7614 features in 1199 cells

neuron_20 <- preprocess_cds(neuron_20, num_dim = 50) 
plot_pc_variance_explained(neuron_20)

neuron_20 <- reduce_dimension(neuron_20, #Takes awhile
                              reduction_method = "UMAP",
                              preprocess_method = "PCA",
                              umap.min_dist = 0.1, # moved to .1 from .15 (fascilitates more separation)
                              umap.n_neighbors = 30) # moved to 30 from 40 (fascilitates more separation)

neuron_20 <- cluster_cells(neuron_20, reduction_method = "UMAP", res = 3e-3) #changed from PCA and 3e-4

plot_cells(neuron_20, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_20, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

#This code needed to go here... vvv
colData(neuron_20)$cluster <- monocle3::clusters(neuron_20)

markers_20 <- top_markers(neuron_20, group_cells_by = "cluster") 
markers_20 %>% filter(cell_group == 1) %>% arrange(desc(marker_score)) %>% head(20)

metadata <- colData(neuron_20) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 1) %>%
  count(cell.subtype, sort = TRUE)  # Sorts by highest countcount(cell.subtype, sort = TRUE)  # Sorts by highest count

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_20)$UMAP_1 <- reducedDims(neuron_20)[["UMAP"]][,1]
colData(neuron_20)$UMAP_2 <- reducedDims(neuron_20)[["UMAP"]][,2]
plot.expr.UMAP(neuron_20, "sbt-1", size = 0.5)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-- 20->19
colData(neuron_20)$Cell.type <- as.character(colData(neuron_20)$cluster)
colData(neuron_20)$Cell.type <- dplyr::recode(colData(neuron_20)$Cell.type,
                                              "1" = "Unannotated",
                                              "2" = "Unannotated",
                                              "3" = "CEP",
                                              "4" = "Unannotated",
                                              "5" = "Unannotated"
                                              
)
table(colData(neuron_20)$Cell.type, exclude = NULL)
# passing annotations up a level to neuron_19
colData(neuron_19)$Cell.type <- as.character(colData(neuron_19)$cluster)
colData(neuron_19)[colnames(neuron_20),]$Cell.type <- colData(neuron_20)$Cell.type
table(colData(neuron_19)$Cell.type)

#-- 19->18
#colData(neuron_19)$Cell.type <- as.character(colData(neuron_19)$cluster)
#colData(neuron_19)$Cell.type <- dplyr::recode(colData(neuron_19)$Cell.type,
idx <- colData(neuron_19)$Cell.type == as.character(colData(neuron_19)$cluster)
colData(neuron_19)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_19)$cluster[idx]),
                                                   "1"  = "Unannotated",
                                                   "2"  = "Unannotated",
                                                   "7"  = "Unannotated",
                                                   "8"  = "Unannotated",
                                                   "9"  = "AV*",
                                                   "10" = "AS*_Cell",
                                                   "11" = "Unannotated",
                                                   "12" = "Unannotated",
                                                   "13" = "Unannotated",
                                                   "14" = "Unannotated",
                                                   "15" = "ADL",
                                                   "16" = "Unannotated",
                                                   "17" = "AWC",
                                                   "18" = "RIA")
table(colData(neuron_19)$Cell.type, exclude = NULL)
# passing annotations up a level to neuron_18
colData(neuron_18)$Cell.type <- as.character(colData(neuron_18)$cluster)
colData(neuron_18)[colnames(neuron_19),]$Cell.type <- colData(neuron_19)$Cell.type
table(colData(neuron_18)$Cell.type)

#-- 18->17
idx <- colData(neuron_18)$Cell.type == as.character(colData(neuron_18)$cluster)
colData(neuron_18)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_18)$cluster[idx]),
                                                   "1"  = "Unannotated",
                                                   "7"  = "Unannotated",
                                                   "11" = "Unannotated",
                                                   "13" = "Unannotated_Neuroblast",
                                                   "16" = "ASK",
                                                   "18" = "AIM",
                                                   "20" = "ASK_Parental",
                                                   "23" = "Unannotated",
                                                   "25" = "URX",
                                                   "27" = "Unannotated",
                                                   "28" = "Unannotated",
                                                   "29" = "PVR",
                                                   "30" = "Unannotated",
                                                   "31" = "Unannotated")
table(colData(neuron_18)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_17
colData(neuron_17)$Cell.type <- as.character(colData(neuron_17)$cluster)
colData(neuron_17)[colnames(neuron_18),]$Cell.type <- colData(neuron_18)$Cell.type
table(colData(neuron_17)$Cell.type)


#-- 17->16
idx <- colData(neuron_17)$Cell.type == as.character(colData(neuron_17)$cluster)
colData(neuron_17)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_17)$cluster[idx]),
                                                   "3"  = "Unannotated",
                                                   "6"  = "DA_DB_DD",
                                                   "17" = "SIA",
                                                   "21" = "Unannotated",
                                                   "23" = "AVE",
                                                   "24" = "Unannotated_Neuroblast",
                                                   "26" = "Unannotated",
                                                   "27" = "Unannotated",
                                                   "28" = "RMD",
                                                   "29" = "URX",
                                                   "30" = "AVL")
table(colData(neuron_17)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_16
colData(neuron_16)$Cell.type <- as.character(colData(neuron_16)$cluster)
colData(neuron_16)[colnames(neuron_17),]$Cell.type <- colData(neuron_17)$Cell.type
table(colData(neuron_16)$Cell.type)

#-- 16->15
idx <- colData(neuron_16)$Cell.type == as.character(colData(neuron_16)$cluster)
colData(neuron_16)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_16)$cluster[idx]),
                                                   "1"  = "Unannotated",
                                                   "4"  = "I_M_Parental",
                                                   "24" = "Unannotated",
                                                   "30" = "Unannotated",
                                                   "34" = "Unannotated",
                                                   "36" = "Unannotated",
                                                   "37" = "Unannotated",
                                                   "38" = "Unannotated",
                                                   "39" = "OLQ")
table(colData(neuron_16)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_15
colData(neuron_15)$Cell.type <- as.character(colData(neuron_15)$cluster)
colData(neuron_15)[colnames(neuron_16),]$Cell.type <- colData(neuron_16)$Cell.type
table(colData(neuron_15)$Cell.type)


#-- 15->14
idx <- colData(neuron_15)$Cell.type == as.character(colData(neuron_15)$cluster)
colData(neuron_15)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_15)$cluster[idx]),
                                                   "21" = "Unannotated",
                                                   "22" = "ASI_Parental",
                                                   "28" = "Unannotated",
                                                   "29" = "OL_Parental",
                                                   "31" = "Unannotated",
                                                   "34" = "AVJ",
                                                   "35" = "ADF",
                                                   "37" = "AS*_Cell",
                                                   "41" = "URA",
                                                   "42" = "Unannotated",
                                                   "45" = "Unannotated",
                                                   "47" = "Unannotated",
                                                   "48" = "Unannotated",
                                                   "49" = "M5")
table(colData(neuron_15)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_14
colData(neuron_14)$Cell.type <- as.character(colData(neuron_14)$cluster)
colData(neuron_14)[colnames(neuron_15),]$Cell.type <- colData(neuron_15)$Cell.type
table(colData(neuron_14)$Cell.type)


#-- 14->13
idx <- colData(neuron_14)$Cell.type == as.character(colData(neuron_14)$cluster)
colData(neuron_14)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_14)$cluster[idx]),
                                                   "4"  = "A*_Cell",
                                                   "13" = "ASG_or_parental",
                                                   "18" = "Il1_Il2_Parental",
                                                   "21" = "PHB_PHA",
                                                   "29" = "CPE_PDE",
                                                   "37" = "I_M",
                                                   "40" = "ADF_AWB",
                                                   "51" = "SAA_SMD",
                                                   "53" = "I5")
table(colData(neuron_14)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_13
colData(neuron_13)$Cell.type <- as.character(colData(neuron_13)$cluster)
colData(neuron_13)[colnames(neuron_14),]$Cell.type <- colData(neuron_14)$Cell.type
table(colData(neuron_13)$Cell.type)


#-- 13->12
idx <- colData(neuron_13)$Cell.type == as.character(colData(neuron_13)$cluster)
colData(neuron_13)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_13)$cluster[idx]),
                                                   "15" = "ASG_AWA",
                                                   "22" = "DA",
                                                   "35" = "IL2",
                                                   "41" = "ALM_PLM",
                                                   "48" = "RMD",
                                                   "53" = "Unannotated",
                                                   "56" = "Unannotated",
                                                   "57" = "Unannotated",
                                                   "58" = "Unannotated")
table(colData(neuron_13)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_12
colData(neuron_12)$Cell.type <- as.character(colData(neuron_12)$cluster)
colData(neuron_12)[colnames(neuron_13),]$Cell.type <- colData(neuron_13)$Cell.type
table(colData(neuron_12)$Cell.type)


#-- 12->11
idx <- colData(neuron_12)$Cell.type == as.character(colData(neuron_12)$cluster)
colData(neuron_12)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_12)$cluster[idx]),
                                                   "31" = "DA_DB",
                                                   "38" = "DA_DB",
                                                   "42" = "IL2",
                                                   "44" = "DA_DB",
                                                   "45" = "RIA",
                                                   "49" = "AWB_AWC",
                                                   "53" = "Unannotated",
                                                   "55" = "BDU",
                                                   "56" = "Unannotated",
                                                   "57" = "URA_URB",
                                                   "58" = "PVM",
                                                   "60" = "Unannotated",
                                                   "62" = "RIS",
                                                   "64" = "Unannotated")
table(colData(neuron_12)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_11
colData(neuron_11)$Cell.type <- as.character(colData(neuron_11)$cluster)
colData(neuron_11)[colnames(neuron_12),]$Cell.type <- colData(neuron_12)$Cell.type
table(colData(neuron_11)$Cell.type)


#-- 12->11 (second recoding block)
idx <- colData(neuron_12)$Cell.type == as.character(colData(neuron_12)$cluster)
colData(neuron_12)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_12)$cluster[idx]),
                                                   "2" = "Unannotated",
                                                   "5" = "Unannotated Parental",
                                                   "6" = "Unannotated",
                                                   "8" = "Unannotated",
                                                   "9" = "Unannotated Parental",
                                                   "10" = "RMD",
                                                   "13" = "SIA",
                                                   "15" = "Unannotated Parental",
                                                   "16" = "Unannotated",
                                                   "17" = "DD",
                                                   "21" = "DD",
                                                   "22" = "ALA",
                                                   "23" = "I5",
                                                   "25" = "RIC/RIM Parental",
                                                   "26" = "PVQ",
                                                   "28" = "AVL",
                                                   "29" = "SAA",
                                                   "31" = "AIA")
table(colData(neuron_12)$Cell.type, exclude = NULL)

colData(neuron_11)$Cell.type <- as.character(colData(neuron_11)$cluster)
colData(neuron_11)[colnames(neuron_12),]$Cell.type <- colData(neuron_12)$Cell.type
table(colData(neuron_11)$Cell.type)


#-- 11->10
idx <- colData(neuron_11)$Cell.type == as.character(colData(neuron_11)$cluster)
colData(neuron_11)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_11)$cluster[idx]),
                                                   "12" = "Unannotated Neuroblast",
                                                   "18" = "Unannotated",
                                                   "39" = "Unannotated Parental",
                                                   "40" = "Unannotated_Parental",
                                                   "46" = "IL1_Parental",
                                                   "47" = "DA_DB",
                                                   "53" = "URA_URB",
                                                   "59" = "IL2",
                                                   "61" = "M4 + M4_Parental",
                                                   "62" = "RIM",
                                                   "65" = "RIH",
                                                   "66" = "Unannotated Parental",
                                                   "68" = "M2")
table(colData(neuron_11)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_10
colData(neuron_10)$Cell.type <- as.character(colData(neuron_10)$cluster)
colData(neuron_10)[colnames(neuron_11),]$Cell.type <- colData(neuron_11)$Cell.type
table(colData(neuron_10)$Cell.type)


#-- 10->9
idx <- colData(neuron_10)$Cell.type == as.character(colData(neuron_10)$cluster)
colData(neuron_10)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_10)$cluster[idx]),
                                                   "5"  = "SMD_RMD",
                                                   "20" = "SMB_SMD",
                                                   "43" = "Unannotated",
                                                   "52" = "URA_URB",
                                                   "56" = "Unannotated",
                                                   "62" = "ASG",
                                                   "65" = "URA_URB",
                                                   "68" = "Unannotated",
                                                   "71" = "RIF",
                                                   "72" = "IL1_IL2",
                                                   "74" = "RME")
table(colData(neuron_10)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_9
colData(neuron_9)$Cell.type <- as.character(colData(neuron_9)$cluster)
colData(neuron_9)[colnames(neuron_10),]$Cell.type <- colData(neuron_10)$Cell.type
table(colData(neuron_9)$Cell.type)


#-- 9->8
idx <- colData(neuron_9)$Cell.type == as.character(colData(neuron_9)$cluster)
colData(neuron_9)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_9)$cluster[idx]),
                                                  "20" = "DA_DB",
                                                  "44" = "OL_Parental (OLQ, OLL)",
                                                  "50" = "RMD",
                                                  "61" = "RID",
                                                  "65" = "AVA",
                                                  "71" = "Unannotated",
                                                  "73" = "BAG",
                                                  "75" = "AWC",
                                                  "78" = "IL2",
                                                  "79" = "I6")
table(colData(neuron_9)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_8
colData(neuron_8)$Cell.type <- as.character(colData(neuron_8)$cluster)
colData(neuron_8)[colnames(neuron_9),]$Cell.type <- colData(neuron_9)$Cell.type
table(colData(neuron_8)$Cell.type)


#-- 8->7
idx <- colData(neuron_8)$Cell.type == as.character(colData(neuron_8)$cluster)
colData(neuron_8)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_8)$cluster[idx]),
                                                  "38" = "DA_DB_DD",
                                                  "57" = "IL1_IL2",
                                                  "64" = "RME",
                                                  "66" = "ASE+ASE_Parental",
                                                  "68" = "ASI",
                                                  "71" = "Neuroblast_AIZ_FLP",
                                                  "73" = "AVD",
                                                  "76" = "RIP",
                                                  "81" = "AVB",
                                                  "84" = "Unannotated_Parental")
table(colData(neuron_8)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_7
colData(neuron_7)$Cell.type <- as.character(colData(neuron_7)$cluster)
colData(neuron_7)[colnames(neuron_8),]$Cell.type <- colData(neuron_8)$Cell.type
table(colData(neuron_7)$Cell.type)


#-- 7->6
idx <- colData(neuron_7)$Cell.type == as.character(colData(neuron_7)$cluster)
colData(neuron_7)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_7)$cluster[idx]),
                                                  "11" = "OLQ",
                                                  "57" = "RIC",
                                                  "61" = "URY",
                                                  "71" = "AVE",
                                                  "76" = "AIZ",
                                                  "82" = "RIM",
                                                  "83" = "RIG",
                                                  "84" = "AVG")
table(colData(neuron_7)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_6
colData(neuron_6)$Cell.type <- as.character(colData(neuron_6)$cluster)
colData(neuron_6)[colnames(neuron_7),]$Cell.type <- colData(neuron_7)$Cell.type
table(colData(neuron_6)$Cell.type)


#-- 6->5
idx <- colData(neuron_6)$Cell.type == as.character(colData(neuron_6)$cluster)
colData(neuron_6)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_6)$cluster[idx]),
                                                  "5"  = "Unannotated (non neuronal)",
                                                  "22" = "DA_DB",
                                                  "36" = "SMB (RIV)",
                                                  "55" = "ADF",
                                                  "78" = "AVA",
                                                  "82" = "RIM",
                                                  "85" = "DVC",
                                                  "86" = "SIA")
table(colData(neuron_6)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_5
colData(neuron_5)$Cell.type <- as.character(colData(neuron_5)$cluster)
colData(neuron_5)[colnames(neuron_6),]$Cell.type <- colData(neuron_6)$Cell.type
table(colData(neuron_5)$Cell.type)


#-- 5->4
idx <- colData(neuron_5)$Cell.type == as.character(colData(neuron_5)$cluster)
colData(neuron_5)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_5)$cluster[idx]),
                                                  "48" = "DA_DB",
                                                  "50" = "IL2",
                                                  "54" = "AVB",
                                                  "71" = "RIB",
                                                  "82" = "DA_DB",
                                                  "87" = "AUA",
                                                  "89" = "ALM_PLM",
                                                  "94" = "SAA",
                                                  "96" = "DVA",
                                                  "99" = "AIA")
table(colData(neuron_5)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_4
colData(neuron_4)$Cell.type <- as.character(colData(neuron_4)$cluster)
colData(neuron_4)[colnames(neuron_5),]$Cell.type <- colData(neuron_5)$Cell.type
table(colData(neuron_4)$Cell.type)
neuron_4_backup = neuron_4

table(colData(neuron_4)$cell.subtype[colData(neuron_4)$Cell.type == "BAG"])
subset(colData(neuron_4), Cell.type == "CEP")



#-- 4->3
idx <- colData(neuron_4)$Cell.type == as.character(colData(neuron_4)$cluster)
colData(neuron_4)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_4)$cluster[idx]),
                                                  "37" = "RME",
                                                  "46" = "AIB",
                                                  "55" = "Unannotated",
                                                  "67" = "Unannotated",
                                                  "79" = "OLL",
                                                  "84" = "ASJ",
                                                  "88" = "AIN",
                                                  "89" = "I5",
                                                  "90" = "AVH")
table(colData(neuron_4)$Cell.type, exclude = NULL)

#------------------------ TEMP FIXING CODE
# Check dimensions
dim(colData(neuron))      # e.g. 38250  <-- number of cells x number of metadata columns
dim(colData(neuron_4))    # e.g. 28255
# See how many cells overlap
shared_cells <- intersect(rownames(colData(neuron)), rownames(colData(neuron_4)))
length(shared_cells)      # e.g. 28255 if that's the intersection

# Pull out the colData as data frames
cd_n  <- colData(neuron)
cd_n4 <- colData(neuron_4)

# Now assign for the shared rows, same column name
cd_n[shared_cells, "Cell.type"] <- cd_n4[shared_cells, "Cell.type"]

# Assign the DataFrame back into the neuron object
colData(neuron) <- cd_n

# Check results
table(colData(neuron)$Cell.type, exclude=NULL)
# Detect any Cell.type that is purely digits
is_numeric_anno <- grepl("^[0-9]+$", colData(neuron)$Cell.type)

# Reassign those to "To_Be_Fixed"
colData(neuron)$Cell.type[is_numeric_anno] <- "To_Be_Fixed"

# Check the updated distribution
table(colData(neuron)$Cell.type, exclude=NULL)


colnames(colData(neuron))
head(as.data.frame(colData(neuron)))
neuron_4
#------------------------ END TEMP FIXING CODE
#------------------------ NEW CLUSTERING ON "To_Be_Fixed""

to_fix_idx <- colData(neuron)$Cell.type == "To_Be_Fixed"
neuron_to_fix <- neuron[, to_fix_idx]

neuron_to_fix <- detect_genes(neuron_to_fix)

neuron_to_fix <- neuron_to_fix[rowData(neuron_to_fix)$num_cells_expressed > 5, ]
neuron_to_fix

neuron_to_fix <- preprocess_cds(neuron_to_fix, num_dim = 50)
plot_pc_variance_explained(neuron_to_fix)

# 5) Reduce dimension using UMAP (adjust umap.min_dist, umap.n_neighbors, etc. as needed)
neuron_to_fix <- reduce_dimension(neuron_to_fix,
                                  reduction_method = "UMAP",
                                  preprocess_method = "PCA",
                                  umap.min_dist = 0.1,
                                  umap.n_neighbors = 30)

neuron_to_fix <- cluster_cells(neuron_to_fix, 
                               reduction_method = "UMAP", 
                               res = 3e-3)

colData(neuron_to_fix)$cluster <- monocle3::clusters(neuron_to_fix)

plot_cells(neuron_to_fix,
           color_cells_by = "cell.subtype",
           group_label_size = 4,
           cell_size = 0.5)

plot_cells(neuron_to_fix,
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

markers_fix <- top_markers(neuron_to_fix, group_cells_by = "cluster")
markers_fix %>% filter(cell_group == 62) %>% arrange(desc(marker_score)) %>% head(20) 
#------------------------------------END NEW CLUSTERING

'''
# passing annotations up a level to neuron_3
colData(neuron_3)$Cell.type <- as.character(colData(neuron_3)$cluster)
colData(neuron_3)[colnames(neuron_4),]$Cell.type <- colData(neuron_4)$Cell.type
table(colData(neuron_3)$Cell.type)

neuron_3_backup = neuron_3
#-- 3->2 
idx <- colData(neuron_3)$Cell.type == as.character(colData(neuron_3)$cluster)
colData(neuron_3)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_3)$cluster[idx]),
                                                  "5"   = "Unannotated",
                                                  "43"  = "AWA",
                                                  "80"  = "AFD",
                                                  "84"  = "Unannotated",
                                                  "91"  = "ADL",
                                                  "93"  = "ALM_PLM",
                                                  "96"  = "PVP",
                                                  "98"  = "RIM",
                                                  "99"  = "AFD",
                                                  "102" = "AIM")
table(colData(neuron_3)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_2
colData(neuron_2)$Cell.type <- as.character(colData(neuron_2)$cluster)
colData(neuron_2)[colnames(neuron_3),]$Cell.type <- colData(neuron_3)$Cell.type
table(colData(neuron_2)$Cell.type)


#-- 2->1
idx <- colData(neuron_2)$Cell.type == as.character(colData(neuron_2)$cluster)
colData(neuron_2)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron_2)$cluster[idx]),
                                                  "32"  = "RIA",
                                                  "68"  = "AVK",
                                                  "70"  = "RMG",
                                                  "71"  = "AWC",
                                                  "77"  = "ADL",
                                                  "81"  = "ADL",
                                                  "85"  = "ASG",
                                                  "90"  = "AFD",
                                                  "91"  = "PVQ",
                                                  "93"  = "DD",
                                                  "97"  = "ALA",
                                                  "98"  = "AVK",
                                                  "100" = "AFD",
                                                  "104" = "AIN",
                                                  "105" = "AIA",
                                                  "106" = "RIS",
                                                  "108" = "PVT",
                                                  "109" = "AVK",
                                                  "110" = "AVL")
table(colData(neuron_2)$Cell.type, exclude = NULL)

# passing annotations up a level to neuron_1
colData(neuron)$Cell.type <- as.character(colData(neuron)$cluster)
colData(neuron)[colnames(neuron_2),]$Cell.type <- colData(neuron_2)$Cell.type
table(colData(neuron)$Cell.type)


#-- Final mapping for neuron_1
idx <- colData(neuron)$Cell.type == as.character(colData(neuron)$cluster)
colData(neuron)$Cell.type[idx] <- dplyr::recode(as.character(colData(neuron)$cluster[idx]),
                                                "13"  = "Unannotated",
                                                "20"  = "OLQ",
                                                "41"  = "URX",
                                                "51"  = "AIY",
                                                "57"  = "AWB",
                                                "62"  = "Unannotated",
                                                "63"  = "CAN",
                                                "69"  = "CAN",
                                                "70"  = "ASH",
                                                "73"  = "ASER, ASEL",
                                                "76"  = "BAG",
                                                "85"  = "ASK",
                                                "87"  = "ASJ",
                                                "94"  = "ASI, ASJ",
                                                "101" = "ASER, ASEL",
                                                "119" = "AUA",
                                                "109" = "IRC",
                                                "117" = "Unannotated")
table(colData(neuron)$Cell.type, exclude = NULL)

colData(neuron)[as.character(colData(neuron)$cluster) == "19",]
'''

colData(neuron_to_fix)$Cell.type <- as.character(colData(neuron_to_fix)$cluster)
colData(neuron_to_fix)$Cell.type <- dplyr::recode(
  colData(neuron_to_fix)$Cell.type,
  "1"  = "Unannotated",
  "2"  = "Unannotated",
  "3"  = "Unannotated",
  "4"  = "OLQ",
  "5"  = "AFD",
  "6"  = "ASE",
  "7"  = "BAG",
  "8"  = "ADL",
  "9"  = "RMG",
  "10" = "URX",
  "11" = "AWA",
  "12" = "AWC",
  "13" = "ADL",
  "14" = "ASK",
  "15" = "CAN",
  "16" = "Unannotated",
  "17" = "ASG",
  "18" = "AIY",
  "19" = "Unannotated",
  "20" = "ASJ",
  "21" = "VD_DD",
  "22" = "AVK_DVA",
  "23" = "RIC",
  "24" = "ASE",
  "25" = "AVK",
  "26" = "ASH",
  "27" = "OLQ",
  "28" = "CAN",
  "29" = "CAN",
  "30" = "ALM_PLM",
  "31" = "ALA",
  "32" = "RIA",
  "33" = "AWB",
  "34" = "AVK",
  "35" = "AFD",
  "36" = "URX",
  "37" = "Unannotated",
  "38" = "PVT",
  "39" = "ASI_ASJ",
  "40" = "PVQ",
  "41" = "ASH",
  "42" = "AWB",
  "43" = "RIM",
  "44" = "RIA",
  "45" = "ADL",
  "46" = "AIA",
  "47" = "RIS",
  "48" = "RIA",
  "49" = "AUA",
  "50" = "AIN",
  "51" = "ASI_ASJ",
  "52" = "Unannotated",
  "53" = "AWA",
  "54" = "AIY",
  "55" = "PVQ",
  "56" = "PVP",
  "57" = "PVP",
  "58" = "AVL",
  "59" = "ASJ",
  "60" = "AWB",
  "61" = "Unannotated",
  "62"="AIY",
  "63"="PVP",
  "64"="CAN"
)

#Combining into final neuron:
cd_fix <- colData(neuron_to_fix)
cd_n   <- colData(neuron)

shared_fix <- intersect(rownames(cd_fix), rownames(cd_n))
length(shared_fix)  # should be ~10k

#Pverwrotomg
cd_n[shared_fix, "Cell.type"] <- cd_fix[shared_fix, "Cell.type"]

colData(neuron) <- cd_n

cd_n4 <- colData(neuron_4)
cd_n  <- colData(neuron)

#Find overlapping
shared_n4 <- intersect(rownames(cd_n4), rownames(cd_n))
length(shared_n4) # ~28k

# Overwrite
cd_n[shared_n4, "Cell.type"] <- cd_n4[shared_n4, "Cell.type"]

# Assign back
colData(neuron) <- cd_n

table(colData(neuron)$Cell.type, exclude=NULL)
final_neuron <- neuron
#=======================================================================================
#Annotating the old annotations:
#A*_Cell
neuron_A_Cell <- final_neuron[, colData(final_neuron)$Cell.type == "A*_Cell"]

neuron_A_Cell <- detect_genes(neuron_A_Cell)
neuron_A_Cell <- neuron_A_Cell[rowData(neuron_A_Cell)$num_cells_expressed > 5, ]
neuron_A_Cell

neuron_A_Cell <- preprocess_cds(neuron_A_Cell, num_dim = 50)
plot_pc_variance_explained(neuron_A_Cell)

neuron_A_Cell <- reduce_dimension(neuron_A_Cell,
                                   reduction_method = "UMAP",
                                   preprocess_method = "PCA",
                                   umap.min_dist = 0.05,
                                   umap.n_neighbors = 15)

neuron_A_Cell <- cluster_cells(neuron_A_Cell, reduction_method = "UMAP", res = 9e-3)
table(monocle3::clusters(neuron_A_Cell))

plot_cells(neuron_A_Cell, 
           color_cells_by = "cluster",
           label_groups_by_cluster = TRUE,
           cell_size = 1)

colData(neuron_A_Cell)$cluster <- monocle3::clusters(neuron_A_Cell)

markers <- top_markers(neuron_A_Cell, group_cells_by = "cluster")
markers %>% filter(cell_group == 3) %>% arrange(desc(marker_score)) %>% head(20) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_A_Cell)$UMAP_1 <- reducedDims(neuron_A_Cell)[["UMAP"]][,1]
colData(neuron_A_Cell)$UMAP_2 <- reducedDims(neuron_A_Cell)[["UMAP"]][,2]
plot.expr.UMAP(neuron_A_Cell, "flp-11", size = 0.5)


colData(neuron_A_Cell)$Cell.type <- as.character(colData(neuron_A_Cell)$cluster)
colData(neuron_A_Cell)$Cell.type <- dplyr::recode(
  colData(neuron_A_Cell)$Cell.type,
  "5" = "RIS",
  "9" = "AIY",

)

cd_ris_aiy  <- colData(neuron_A_Cell)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_ris_aiy), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_ris_aiy[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)


#ADF_AWB
neuron_ADF_AWB <- final_neuron[, colData(final_neuron)$Cell.type == "ADF_AWB"]

neuron_ADF_AWB <- detect_genes(neuron_ADF_AWB)
neuron_ADF_AWB <- neuron_ADF_AWB[rowData(neuron_ADF_AWB)$num_cells_expressed > 5, ]
neuron_ADF_AWB

neuron_ADF_AWB <- preprocess_cds(neuron_ADF_AWB, num_dim = 50)
plot_pc_variance_explained(neuron_ADF_AWB)

neuron_ADF_AWB <- reduce_dimension(neuron_ADF_AWB,
                                  reduction_method = "UMAP",
                                  preprocess_method = "PCA",
                                  umap.min_dist = 0.05,
                                  umap.n_neighbors = 15)

neuron_ADF_AWB <- cluster_cells(neuron_ADF_AWB, reduction_method = "UMAP", res = 3e-1)
table(monocle3::clusters(neuron_ADF_AWB))

plot_cells(neuron_ADF_AWB, 
           color_cells_by = "cluster",
           label_groups_by_cluster = TRUE,
           cell_size = 1)

colData(neuron_ADF_AWB)$cluster <- monocle3::clusters(neuron_ADF_AWB)

markers <- top_markers(neuron_ADF_AWB, group_cells_by = "cluster")
markers %>% filter(cell_group == 5) %>% arrange(desc(marker_score)) %>% head(20) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_ADF_AWB)$UMAP_1 <- reducedDims(neuron_ADF_AWB)[["UMAP"]][,1]
colData(neuron_ADF_AWB)$UMAP_2 <- reducedDims(neuron_ADF_AWB)[["UMAP"]][,2]
plot.expr.UMAP(neuron_ADF_AWB, "flp-11", size = 0.5)

#Unannotated
neuron_Unannotated <- final_neuron[, colData(final_neuron)$Cell.type == "Unannotated"]

neuron_Unannotated <- detect_genes(neuron_Unannotated)
neuron_Unannotated <- neuron_Unannotated[rowData(neuron_Unannotated)$num_cells_expressed > 5, ]
neuron_Unannotated

neuron_Unannotated <- preprocess_cds(neuron_Unannotated, num_dim = 50)
plot_pc_variance_explained(neuron_Unannotated)

neuron_Unannotated <- reduce_dimension(neuron_Unannotated,
                                  reduction_method = "UMAP",
                                  preprocess_method = "PCA",
                                  umap.min_dist = 0.1,
                                  umap.n_neighbors = 30)

neuron_Unannotated <- cluster_cells(neuron_Unannotated, reduction_method = "UMAP", res = 3e-3)
table(monocle3::clusters(neuron_Unannotated))

plot_cells(neuron_Unannotated, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(neuron_Unannotated, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = .5)

colData(neuron_Unannotated)$cluster <- monocle3::clusters(neuron_Unannotated)

markers <- top_markers(neuron_Unannotated, group_cells_by = "cluster")
markers %>% filter(cell_group == 47) %>% arrange(desc(marker_score)) %>% head(20) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_Unannotated)$UMAP_1 <- reducedDims(neuron_Unannotated)[["UMAP"]][,1]
colData(neuron_Unannotated)$UMAP_2 <- reducedDims(neuron_Unannotated)[["UMAP"]][,2]
plot.expr.UMAP(neuron_Unannotated, "flp-11", size = 0.5)

colData(neuron_Unannotated)$Cell.type <- as.character(colData(neuron_Unannotated)$cluster)
colData(neuron_Unannotated)$Cell.type <- dplyr::recode(
  colData(neuron_Unannotated)$Cell.type,
  "1" = "Unannotated",
  "2" = "Unannotated",
  "3" = "Unannotated",
  "4" = "Unannotated",
  "5" = "Unannotated",
  "6" = "Unannotated",
  "7" = "Unannotated",
  "8" = "Unannotated",
  "9" = "Unannotated",
  "10" = "Unannotated",
  "11" = "Unannotated",
  "12" = "Unannotated",
  "13" = "Unannotated",
  "14" = "Unannotated",
  "15" = "Unannotated",
  "16" = "Unannotated",
  "17" = "Unannotated",
  "18" = "Unannotated",
  "19" = "Unannotated",
  "20" = "AFD",
  "21" = "Unannotated",
  "22" = "Unannotated",
  "23" = "Unannotated",
  "24" = "Unannotated",
  "25" = "Unannotated",
  "26" = "Unannotated",
  "27" = "Unannotated",
  "28" = "Unannotated",
  "29" = "Unannotated",
  "30" = "Unannotated",
  "31" = "Unannotated",
  "32" = "Unannotated",
  "33" = "Unannotated",
  "34" = "Unannotated",
  "35" = "Unannotated",
  "36" = "Unannotated",
  "37" = "Unannotated",
  "38" = "Unannotated",
  "39" = "Unannotated",
  "40" = "Unannotated",
  "41" = "Unannotated",
  "42" = "Unannotated",
  "43" = "Unannotated",
  "44" = "Unannotated",
  "45" = "Unannotated",
  "46" = "Unannotated",
  "47" = "AVH",
  "48" = "Unannotated",
  "49" = "Unannotated",
  "50" = "Unannotated",
  "51" = "Unannotated",
  "52" = "ADA",
  "53" = "Unannotated"
)

cd_new  <- colData(neuron_Unannotated)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)


#URA_URB
neuron_URA_URB <- final_neuron[, colData(final_neuron)$Cell.type == "URA_URB"]

neuron_URA_URB <- detect_genes(neuron_URA_URB)
neuron_URA_URB <- neuron_URA_URB[rowData(neuron_URA_URB)$num_cells_expressed > 5, ]
neuron_URA_URB

neuron_URA_URB <- preprocess_cds(neuron_URA_URB, num_dim = 50)
plot_pc_variance_explained(neuron_URA_URB)

neuron_URA_URB <- reduce_dimension(neuron_URA_URB,
                                       reduction_method = "UMAP",
                                       preprocess_method = "PCA",
                                       umap.min_dist = 0.05,
                                       umap.n_neighbors = 15)

neuron_URA_URB <- cluster_cells(neuron_URA_URB, reduction_method = "UMAP", res = 3e-2)
table(monocle3::clusters(neuron_URA_URB))

plot_cells(neuron_URA_URB, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_URA_URB, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_URA_URB)$cluster <- monocle3::clusters(neuron_URA_URB)

markers <- top_markers(neuron_URA_URB, group_cells_by = "cluster")
markers %>% filter(cell_group == 6) %>% arrange(desc(marker_score)) %>% head(20) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_URA_URB)$UMAP_1 <- reducedDims(neuron_URA_URB)[["UMAP"]][,1]
colData(neuron_URA_URB)$UMAP_2 <- reducedDims(neuron_URA_URB)[["UMAP"]][,2]
plot.expr.UMAP(neuron_URA_URB, "flp-11", size = 0.5)

colData(neuron_URA_URB)$Cell.type <- as.character(colData(neuron_URA_URB)$cluster)
colData(neuron_URA_URB)$Cell.type <- dplyr::recode(
  colData(neuron_URA_URB)$Cell.type,
  "1" = "URB",
  "2" = "URA_URB",
  "3" = "Unannotated",
  "4" = "Unannotated",
  "5" = "URY?",
  "6" = "URA_URB"
)

cd_new  <- colData(neuron_URA_URB)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)


#Unannotated_Parental2
neuron_Unannotated_Parental <- final_neuron[, colData(final_neuron)$Cell.type == "Unannotated_Parental"]

neuron_Unannotated_Parental <- detect_genes(neuron_Unannotated_Parental)
neuron_Unannotated_Parental <- neuron_Unannotated_Parental[rowData(neuron_Unannotated_Parental)$num_cells_expressed > 5, ]
neuron_Unannotated_Parental

neuron_Unannotated_Parental <- preprocess_cds(neuron_Unannotated_Parental, num_dim = 50)
plot_pc_variance_explained(neuron_Unannotated_Parental)

neuron_Unannotated_Parental <- reduce_dimension(neuron_Unannotated_Parental,
                                   reduction_method = "UMAP",
                                   preprocess_method = "PCA",
                                   umap.min_dist = 0.05,
                                   umap.n_neighbors = 15)

neuron_Unannotated_Parental <- cluster_cells(neuron_Unannotated_Parental, reduction_method = "UMAP", res = 6e-2)
table(monocle3::clusters(neuron_Unannotated_Parental))

plot_cells(neuron_Unannotated_Parental, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_Unannotated_Parental, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_Unannotated_Parental)$cluster <- monocle3::clusters(neuron_Unannotated_Parental)

markers <- top_markers(neuron_Unannotated_Parental, group_cells_by = "cluster")
markers %>% filter(cell_group == 4) %>% arrange(desc(marker_score)) %>% head(20) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_Unannotated_Parental)$UMAP_1 <- reducedDims(neuron_Unannotated_Parental)[["UMAP"]][,1]
colData(neuron_Unannotated_Parental)$UMAP_2 <- reducedDims(neuron_Unannotated_Parental)[["UMAP"]][,2]
plot.expr.UMAP(neuron_Unannotated_Parental, "flp-11", size = 0.5)

colData(neuron_Unannotated_Parental)$Cell.type <- as.character(colData(neuron_Unannotated_Parental)$cluster)
colData(neuron_Unannotated_Parental)$Cell.type <- dplyr::recode(
  colData(neuron_Unannotated_Parental)$Cell.type,
  "1" = "Neuronblast_ASG_AWA",
  "2" = "Unannotated_Neuroblast",
  "3" = "ADL_Parental",
  "4" = "Unannotated_Neuroblast"
  
)

cd_new  <- colData(neuron_Unannotated_Parental)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)


#Unannotated_Parental2
neuron_Unannotated_Parental2 <- final_neuron[, colData(final_neuron)$Cell.type == "Unannotated Parental"]

neuron_Unannotated_Parental2 <- detect_genes(neuron_Unannotated_Parental2)
neuron_Unannotated_Parental2 <- neuron_Unannotated_Parental2[rowData(neuron_Unannotated_Parental2)$num_cells_expressed > 5, ]
neuron_Unannotated_Parental2

neuron_Unannotated_Parental2 <- preprocess_cds(neuron_Unannotated_Parental2, num_dim = 50)
plot_pc_variance_explained(neuron_Unannotated_Parental2)

neuron_Unannotated_Parental2 <- reduce_dimension(neuron_Unannotated_Parental2,
                                                reduction_method = "UMAP",
                                                preprocess_method = "PCA",
                                                umap.min_dist = 0.05,
                                                umap.n_neighbors = 15)

neuron_Unannotated_Parental2 <- cluster_cells(neuron_Unannotated_Parental2, reduction_method = "UMAP", res = 6e-2)
table(monocle3::clusters(neuron_Unannotated_Parental2))

plot_cells(neuron_Unannotated_Parental2, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_Unannotated_Parental2, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_Unannotated_Parental2)$cluster <- monocle3::clusters(neuron_Unannotated_Parental2)

markers <- top_markers(neuron_Unannotated_Parental2, group_cells_by = "cluster")
markers %>% filter(cell_group == 4) %>% arrange(desc(marker_score)) %>% head(20) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_Unannotated_Parental2)$UMAP_1 <- reducedDims(neuron_Unannotated_Parental2)[["UMAP"]][,1]
colData(neuron_Unannotated_Parental2)$UMAP_2 <- reducedDims(neuron_Unannotated_Parental2)[["UMAP"]][,2]
plot.expr.UMAP(neuron_Unannotated_Parental, "flp-11", size = 0.5)

colData(neuron_Unannotated_Parental2)$Cell.type <- as.character(colData(neuron_Unannotated_Parental2)$cluster)
colData(neuron_Unannotated_Parental2)$Cell.type <- dplyr::recode(
  colData(neuron_Unannotated_Parental2)$Cell.type,
  "1" = "Unannotated_Parental",
  "2" = "PLM_ALM_Grandparent",
  "3" = "Neuroblast_ALM_BDU",
  "4" = "Unannotated",
  "5" = "Neuroblast_ALM_BDU"
)

cd_new  <- colData(neuron_Unannotated_Parental2)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)

#Unannotated Neuroblast
neuron_Unannotated_Neuroblast <- final_neuron[, colData(final_neuron)$Cell.type == "Unannotated Neuroblast"]

neuron_Unannotated_Neuroblast <- detect_genes(neuron_Unannotated_Neuroblast)
neuron_Unannotated_Neuroblast <- neuron_Unannotated_Neuroblast[rowData(neuron_Unannotated_Neuroblast)$num_cells_expressed > 5, ]
neuron_Unannotated_Neuroblast

neuron_Unannotated_Neuroblast <- preprocess_cds(neuron_Unannotated_Neuroblast, num_dim = 50)
plot_pc_variance_explained(neuron_Unannotated_Neuroblast)

neuron_Unannotated_Neuroblast <- reduce_dimension(neuron_Unannotated_Neuroblast,
                                                 reduction_method = "UMAP",
                                                 preprocess_method = "PCA",
                                                 umap.min_dist = 0.05,
                                                 umap.n_neighbors = 15)

neuron_Unannotated_Neuroblast <- cluster_cells(neuron_Unannotated_Neuroblast, reduction_method = "UMAP", res = 6e-2)
table(monocle3::clusters(neuron_Unannotated_Neuroblast))

plot_cells(neuron_Unannotated_Neuroblast, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_Unannotated_Neuroblast, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_Unannotated_Neuroblast)$cluster <- monocle3::clusters(neuron_Unannotated_Neuroblast)

markers <- top_markers(neuron_Unannotated_Neuroblast, group_cells_by = "cluster")
markers %>% filter(cell_group == 7) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_Unannotated_Neuroblast)$cluster

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_Unannotated_Neuroblast)$UMAP_1 <- reducedDims(neuron_Unannotated_Neuroblast)[["UMAP"]][,1]
colData(neuron_Unannotated_Neuroblast)$UMAP_2 <- reducedDims(neuron_Unannotated_Neuroblast)[["UMAP"]][,2]
plot.expr.UMAP(neuron_Unannotated_Parental, "flp-11", size = 0.5)

colData(neuron_Unannotated_Neuroblast)$Cell.type <- as.character(colData(neuron_Unannotated_Neuroblast)$cluster)
colData(neuron_Unannotated_Neuroblast)$Cell.type <- dplyr::recode(
  colData(neuron_Unannotated_Neuroblast)$Cell.type,
  "1" = "Neuroblast_ADE_ADA",
  "2" = "Neuroblast_HSN_PHB",
  "3" = "Neuroblast_URX_CEPDx",
  "4" = "Unannotated",
  "5" = "Unannotated",
  "6" = "Neuroblast_URX_CEPDx",
  "7" = "Unannotated",
)

cd_new  <- colData(neuron_Unannotated_Neuroblast)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)

#Unannotated Neuroblast2
neuron_Unannotated_Neuroblast2 <- final_neuron[, colData(final_neuron)$Cell.type == "Unannotated_Neuroblast"]

neuron_Unannotated_Neuroblast2 <- detect_genes(neuron_Unannotated_Neuroblast2)
neuron_Unannotated_Neuroblast2 <- neuron_Unannotated_Neuroblast2[rowData(neuron_Unannotated_Neuroblast2)$num_cells_expressed > 5, ]
neuron_Unannotated_Neuroblast2

neuron_Unannotated_Neuroblast2 <- preprocess_cds(neuron_Unannotated_Neuroblast2, num_dim = 50)
plot_pc_variance_explained(neuron_Unannotated_Neuroblast2)

neuron_Unannotated_Neuroblast2 <- reduce_dimension(neuron_Unannotated_Neuroblast2,
                                                  reduction_method = "UMAP",
                                                  preprocess_method = "PCA",
                                                  umap.min_dist = 0.05,
                                                  umap.n_neighbors = 15)

neuron_Unannotated_Neuroblast2 <- cluster_cells(neuron_Unannotated_Neuroblast2, reduction_method = "UMAP", res = 6e-2)
table(monocle3::clusters(neuron_Unannotated_Neuroblast2))

plot_cells(neuron_Unannotated_Neuroblast2, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_Unannotated_Neuroblast2, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_Unannotated_Neuroblast2)$cluster <- monocle3::clusters(neuron_Unannotated_Neuroblast2)

markers <- top_markers(neuron_Unannotated_Neuroblast2, group_cells_by = "cluster")
markers %>% filter(cell_group == 5) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_Unannotated_Neuroblast2)$cluster

metadata <- colData(neuron_Unannotated_Neuroblast2) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 1) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_Unannotated_Neuroblast2)$UMAP_1 <- reducedDims(neuron_Unannotated_Neuroblast2)[["UMAP"]][,1]
colData(neuron_Unannotated_Neuroblast2)$UMAP_2 <- reducedDims(neuron_Unannotated_Neuroblast2)[["UMAP"]][,2]
plot.expr.UMAP(neuron_Unannotated_Neuroblast2, "flp-11", size = 0.5)

colData(neuron_Unannotated_Neuroblast2)$Cell.type <- as.character(colData(neuron_Unannotated_Neuroblast2)$cluster)
colData(neuron_Unannotated_Neuroblast2)$Cell.type <- dplyr::recode(
  colData(neuron_Unannotated_Neuroblast2)$Cell.type,
  "1" = "Neuroblast_ALA_RMED",
  "2" = "PVQ_Parental",
  "3" = "Unannotated_Neuroblast",
  "4" = "Neuroblast_PVC_LUA",
  "5" = "Unannotated",
  "6" = "PVQ_Parental"
)

cd_new  <- colData(neuron_Unannotated_Neuroblast2)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))

#Unannotated (non neuronal)
neuron_Unannotated_nonneuronal <- final_neuron[, colData(final_neuron)$Cell.type == "Unannotated (non neuronal)"]

neuron_Unannotated_nonneuronal <- detect_genes(neuron_Unannotated_nonneuronal)
neuron_Unannotated_nonneuronal <- neuron_Unannotated_nonneuronal[rowData(neuron_Unannotated_nonneuronal)$num_cells_expressed > 5, ]
neuron_Unannotated_nonneuronal

neuron_Unannotated_nonneuronal <- preprocess_cds(neuron_Unannotated_nonneuronal, num_dim = 50)
plot_pc_variance_explained(neuron_Unannotated_nonneuronal)

neuron_Unannotated_nonneuronal <- reduce_dimension(neuron_Unannotated_nonneuronal,
                                                   reduction_method = "UMAP",
                                                   preprocess_method = "PCA",
                                                   umap.min_dist = 0.05,
                                                   umap.n_neighbors = 15)

neuron_Unannotated_nonneuronal <- cluster_cells(neuron_Unannotated_nonneuronal, reduction_method = "UMAP", res = 6e-2)
table(monocle3::clusters(neuron_Unannotated_nonneuronal))

plot_cells(neuron_Unannotated_nonneuronal, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_Unannotated_nonneuronal, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_Unannotated_nonneuronal)$cluster <- monocle3::clusters(neuron_Unannotated_nonneuronal)

markers <- top_markers(neuron_Unannotated_nonneuronal, group_cells_by = "cluster")
markers %>% filter(cell_group == 9) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_Unannotated_nonneuronal)$cluster

metadata <- colData(neuron_Unannotated_nonneuronal) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 1) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_Unannotated_nonneuronal)$UMAP_1 <- reducedDims(neuron_Unannotated_nonneuronal)[["UMAP"]][,1]
colData(neuron_Unannotated_nonneuronal)$UMAP_2 <- reducedDims(neuron_Unannotated_nonneuronal)[["UMAP"]][,2]
plot.expr.UMAP(neuron_Unannotated_nonneuronal, "flp-11", size = 0.5)

colData(neuron_Unannotated_nonneuronal)$Cell.type <- as.character(colData(neuron_Unannotated_nonneuronal)$cluster)
colData(neuron_Unannotated_nonneuronal)$Cell.type <- dplyr::recode(
  colData(neuron_Unannotated_nonneuronal)$Cell.type,
  "1" = "Unannotated",
  "2" = "Unannotated",
  "3" = "Unannotated",
  "4" = "Unannotated",
  "5" = "Unannotated",
  "6" = "Unannotated",
  "7" = "Unannotated",
  "8" = "Unannotated",
  "9" = "AVB"
)

cd_new  <- colData(neuron_Unannotated_nonneuronal)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))

#SMD_RMD
neuron_SMD_RMD <- final_neuron[, colData(final_neuron)$Cell.type == "SMD_RMD"]

neuron_SMD_RMD <- detect_genes(neuron_SMD_RMD)
neuron_SMD_RMD <- neuron_SMD_RMD[rowData(neuron_SMD_RMD)$num_cells_expressed > 5, ]
neuron_SMD_RMD

neuron_SMD_RMD <- preprocess_cds(neuron_SMD_RMD, num_dim = 50)
plot_pc_variance_explained(neuron_SMD_RMD)

neuron_SMD_RMD <- reduce_dimension(neuron_SMD_RMD,
                                                   reduction_method = "UMAP",
                                                   preprocess_method = "PCA",
                                                   umap.min_dist = 0.05,
                                                   umap.n_neighbors = 15)

neuron_SMD_RMD <- cluster_cells(neuron_SMD_RMD, reduction_method = "UMAP", res = 6e-2)
table(monocle3::clusters(neuron_SMD_RMD))

plot_cells(neuron_SMD_RMD, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_SMD_RMD, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_SMD_RMD)$cluster <- monocle3::clusters(neuron_SMD_RMD)

markers <- top_markers(neuron_SMD_RMD, group_cells_by = "cluster")
markers %>% filter(cell_group == 1) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_SMD_RMD)$cluster

metadata <- colData(neuron_SMD_RMD) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 1) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_SMD_RMD)$UMAP_1 <- reducedDims(neuron_SMD_RMD)[["UMAP"]][,1]
colData(neuron_SMD_RMD)$UMAP_2 <- reducedDims(neuron_SMD_RMD)[["UMAP"]][,2]
plot.expr.UMAP(neuron_SMD_RMD, "flp-11", size = 0.5)

colData(neuron_SMD_RMD)$Cell.type <- as.character(colData(neuron_SMD_RMD)$cluster)
colData(neuron_SMD_RMD)$Cell.type <- dplyr::recode(
  colData(neuron_SMD_RMD)$Cell.type,
  "1" = "SMD_RMD",
  "2" = "SMD_RMD",
  "3" = "SMD_RMD",
  "4" = "SMD_RMD",
  "5" = "SMD_RMD",
  "6" = "SMD_RMD",
  "7" = "SAA"
)

cd_new  <- colData(neuron_SMD_RMD)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))

#SAA_SMD
neuron_SAA_SMD<- final_neuron[, colData(final_neuron)$Cell.type == "SAA_SMD"]

neuron_SAA_SMD <- detect_genes(neuron_SAA_SMD)
neuron_SAA_SMD <- neuron_SAA_SMD[rowData(neuron_SAA_SMD)$num_cells_expressed > 5, ]
neuron_SAA_SMD

neuron_SAA_SMD <- preprocess_cds(neuron_SAA_SMD, num_dim = 50)
plot_pc_variance_explained(neuron_SAA_SMD)

neuron_SAA_SMD <- reduce_dimension(neuron_SAA_SMD,
                                   reduction_method = "UMAP",
                                   preprocess_method = "PCA",
                                   umap.min_dist = 0.05,
                                   umap.n_neighbors = 15)

neuron_SAA_SMD <- cluster_cells(neuron_SAA_SMD, reduction_method = "UMAP", res = 90e-2)
table(monocle3::clusters(neuron_SAA_SMD))

plot_cells(neuron_SAA_SMD, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_SAA_SMD, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_SAA_SMD)$cluster <- monocle3::clusters(neuron_SAA_SMD)

markers <- top_markers(neuron_SAA_SMD, group_cells_by = "cluster")
markers %>% filter(cell_group == 3) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_SAA_SMD)$cluster

metadata <- colData(neuron_SAA_SMD) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 1) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_SAA_SMD)$UMAP_1 <- reducedDims(neuron_SAA_SMD)[["UMAP"]][,1]
colData(neuron_SAA_SMD)$UMAP_2 <- reducedDims(neuron_SAA_SMD)[["UMAP"]][,2]
plot.expr.UMAP(neuron_SAA_SMD, "flp-11", size = 0.5)

colData(neuron_SAA_SMD)$Cell.type <- as.character(colData(neuron_SAA_SMD)$cluster)
colData(neuron_SAA_SMD)$Cell.type <- dplyr::recode(
  colData(neuron_SAA_SMD)$Cell.type,
  "1" = "SMD",
  "2" = "SAA",
  "3" = "SAA"
  
)

cd_new  <- colData(neuron_SAA_SMD)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))

#DA_DB_DD
neuron_DA_DB_DD<- final_neuron[, colData(final_neuron)$Cell.type == "DA_DB_DD"]

neuron_DA_DB_DD <- detect_genes(neuron_DA_DB_DD)
neuron_DA_DB_DD <- neuron_DA_DB_DD[rowData(neuron_DA_DB_DD)$num_cells_expressed > 5, ]
neuron_DA_DB_DD

neuron_DA_DB_DD <- preprocess_cds(neuron_DA_DB_DD, num_dim = 10)
plot_pc_variance_explained(neuron_DA_DB_DD)

neuron_DA_DB_DD <- reduce_dimension(neuron_DA_DB_DD,
                                   reduction_method = "UMAP",
                                   preprocess_method = "PCA",
                                   umap.min_dist = 0.05,
                                   umap.n_neighbors = 15)

neuron_DA_DB_DD <- cluster_cells(neuron_DA_DB_DD, reduction_method = "UMAP", res = 1e-3)
table(monocle3::clusters(neuron_DA_DB_DD))

plot_cells(neuron_DA_DB_DD, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_DA_DB_DD, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_DA_DB_DD)$cluster <- monocle3::clusters(neuron_DA_DB_DD)

markers <- top_markers(neuron_DA_DB_DD, group_cells_by = "cluster")
markers %>% filter(cell_group == 1) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_DA_DB_DD)$cluster

metadata <- colData(neuron_DA_DB_DD) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 1) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_DA_DB_DD)$UMAP_1 <- reducedDims(neuron_DA_DB_DD)[["UMAP"]][,1]
colData(neuron_DA_DB_DD)$UMAP_2 <- reducedDims(neuron_DA_DB_DD)[["UMAP"]][,2]
plot.expr.UMAP(neuron_DA_DB_DD, "ttr-39", size = 0.5)

colData(neuron_DA_DB_DD)$Cell.type <- as.character(colData(neuron_DA_DB_DD)$cluster)
colData(neuron_DA_DB_DD)$Cell.type <- dplyr::recode(
  colData(neuron_DA_DB_DD)$Cell.type,
  "1" = "DA_DB", 
  "2" = "PDB" #?
)

cd_new  <- colData(neuron_DA_DB_DD)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))


#DA_DB
neuron_DA_DB <- final_neuron[, colData(final_neuron)$Cell.type == "DA_DB"]

neuron_DA_DB <- detect_genes(neuron_DA_DB)
neuron_DA_DB <- neuron_DA_DB[rowData(neuron_DA_DB)$num_cells_expressed > 5, ]
neuron_DA_DB

neuron_DA_DB <- preprocess_cds(neuron_DA_DB, num_dim = 20)
plot_pc_variance_explained(neuron_DA_DB)

neuron_DA_DB <- reduce_dimension(neuron_DA_DB,
                                    reduction_method = "UMAP",
                                    preprocess_method = "PCA",
                                    umap.min_dist = 0.05,
                                    umap.n_neighbors = 15)

neuron_DA_DB <- cluster_cells(neuron_DA_DB, reduction_method = "UMAP", res = 2e-3)
table(monocle3::clusters(neuron_DA_DB))

plot_cells(neuron_DA_DB, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_DA_DB, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_DA_DB)$cluster <- monocle3::clusters(neuron_DA_DB)

markers <- top_markers(neuron_DA_DB, group_cells_by = "cluster")
markers %>% filter(cell_group == 1) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_DA_DB)$cluster

metadata <- colData(neuron_DA_DB) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 1) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_DA_DB)$UMAP_1 <- reducedDims(neuron_DA_DB)[["UMAP"]][,1]
colData(neuron_DA_DB)$UMAP_2 <- reducedDims(neuron_DA_DB)[["UMAP"]][,2]
plot.expr.UMAP(neuron_DA_DB, "F44E5.5", size = 0.5)
#vab-7 is DB
#unc-4 is DA

colData(neuron_DA_DB)$Cell.type <- as.character(colData(neuron_DA_DB)$cluster)
colData(neuron_DA_DB)$Cell.type <- dplyr::recode(
  colData(neuron_DA_DB)$Cell.type,
  "1" = "DA_DB",
  "2" = "DA_DB",
  "3" = "DA_DB",
  "4" = "DA_DB",
  "5" = "DA_DB",
  "6" = "DA_DB",
  "7" = "DA_DB",
  "8" = "DB",
  "9" = "DA_DB",
  "10" = "DB",
  "11" = "DA",
  "12" = "DA",
  "13" = "DB"
)

cd_new  <- colData(neuron_DA_DB)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))

#AVK_DVA
neuron_AVK_DVA <- final_neuron[, colData(final_neuron)$Cell.type == "AVK_DVA"]

neuron_AVK_DVA <- detect_genes(neuron_AVK_DVA)
neuron_AVK_DVA <- neuron_AVK_DVA[rowData(neuron_AVK_DVA)$num_cells_expressed > 5, ]
neuron_AVK_DVA

neuron_AVK_DVA <- preprocess_cds(neuron_AVK_DVA, num_dim = 20)
plot_pc_variance_explained(neuron_AVK_DVA)

neuron_AVK_DVA <- reduce_dimension(neuron_AVK_DVA,
                                 reduction_method = "UMAP",
                                 preprocess_method = "PCA",
                                 umap.min_dist = 0.05,
                                 umap.n_neighbors = 15)

neuron_AVK_DVA <- cluster_cells(neuron_AVK_DVA, reduction_method = "UMAP", res = 2e-3)
table(monocle3::clusters(neuron_AVK_DVA))

plot_cells(neuron_AVK_DVA, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_AVK_DVA, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_AVK_DVA)$cluster <- monocle3::clusters(neuron_AVK_DVA)

markers <- top_markers(neuron_AVK_DVA, group_cells_by = "cluster")
markers %>% filter(cell_group == 2) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_AVK_DVA)$cluster

metadata <- colData(neuron_AVK_DVA) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 1) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_AVK_DVA)$UMAP_1 <- reducedDims(neuron_AVK_DVA)[["UMAP"]][,1]
colData(neuron_AVK_DVA)$UMAP_2 <- reducedDims(neuron_AVK_DVA)[["UMAP"]][,2]
plot.expr.UMAP(neuron_AVK_DVA, "F44E5.5", size = 0.5)

colData(neuron_AVK_DVA)$Cell.type <- as.character(colData(neuron_AVK_DVA)$cluster)
colData(neuron_AVK_DVA)$Cell.type <- dplyr::recode(
  colData(neuron_AVK_DVA)$Cell.type,
  "1" = "AVK",
  "2" = "DVA"
)

cd_new  <- colData(neuron_AVK_DVA)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))


#AWB_AWC
neuron_AWB_AWC <- final_neuron[, colData(final_neuron)$Cell.type == "AWB_AWC"]

neuron_AWB_AWC <- detect_genes(neuron_AWB_AWC)
neuron_AWB_AWC <- neuron_AWB_AWC[rowData(neuron_AWB_AWC)$num_cells_expressed > 5, ]
neuron_AWB_AWC

neuron_AWB_AWC <- preprocess_cds(neuron_AWB_AWC, num_dim = 20)
plot_pc_variance_explained(neuron_AWB_AWC)

neuron_AWB_AWC <- reduce_dimension(neuron_AWB_AWC,
                                   reduction_method = "UMAP",
                                   preprocess_method = "PCA",
                                   umap.min_dist = 0.05,
                                   umap.n_neighbors = 15)

neuron_AWB_AWC <- cluster_cells(neuron_AWB_AWC, reduction_method = "UMAP", res = 10e-3)
table(monocle3::clusters(neuron_AWB_AWC))

plot_cells(neuron_AWB_AWC, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_AWB_AWC, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_AWB_AWC)$cluster <- monocle3::clusters(neuron_AWB_AWC)

markers <- top_markers(neuron_AWB_AWC, group_cells_by = "cluster")
markers %>% filter(cell_group == 1) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_AWB_AWC)$cluster

metadata <- colData(neuron_AWB_AWC) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 1) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_AWB_AWC)$UMAP_1 <- reducedDims(neuron_AWB_AWC)[["UMAP"]][,1]
colData(neuron_AWB_AWC)$UMAP_2 <- reducedDims(neuron_AWB_AWC)[["UMAP"]][,2]
plot.expr.UMAP(neuron_AWB_AWC, "C33A12.4", size = 2)

colData(neuron_AWB_AWC)$Cell.type <- as.character(colData(neuron_AWB_AWC)$cluster)
colData(neuron_AWB_AWC)$Cell.type <- dplyr::recode(
  colData(neuron_AWB_AWC)$Cell.type,
  "1" = "AS*_Cell",
  "2" = "ASG"
)

cd_new  <- colData(neuron_AWB_AWC)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))


#CEP_PDE
neuron_CEP_PDE <- final_neuron[, colData(final_neuron)$Cell.type == "CPE_PDE"] #I misstypes this initially 

neuron_CEP_PDE <- detect_genes(neuron_CEP_PDE)
neuron_CEP_PDE <- neuron_CEP_PDE[rowData(neuron_CEP_PDE)$num_cells_expressed > 5, ]
neuron_CEP_PDE

neuron_CEP_PDE <- preprocess_cds(neuron_CEP_PDE, num_dim = 20)
plot_pc_variance_explained(neuron_CEP_PDE)

neuron_CEP_PDE <- reduce_dimension(neuron_CEP_PDE,
                                   reduction_method = "UMAP",
                                   preprocess_method = "PCA",
                                   umap.min_dist = 0.05,
                                   umap.n_neighbors = 15)

neuron_CEP_PDE <- cluster_cells(neuron_CEP_PDE, reduction_method = "UMAP", res = 10e-3)
table(monocle3::clusters(neuron_CEP_PDE))

plot_cells(neuron_CEP_PDE, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_CEP_PDE, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_CEP_PDE)$cluster <- monocle3::clusters(neuron_CEP_PDE)

markers <- top_markers(neuron_CEP_PDE, group_cells_by = "cluster")
markers %>% filter(cell_group == 1) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_CEP_PDE)$cluster

metadata <- colData(neuron_CEP_PDE) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 1) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_CEP_PDE)$UMAP_1 <- reducedDims(neuron_CEP_PDE)[["UMAP"]][,1]
colData(neuron_CEP_PDE)$UMAP_2 <- reducedDims(neuron_CEP_PDE)[["UMAP"]][,2]
plot.expr.UMAP(neuron_CEP_PDE, "C33A12.4", size = 2)

colData(neuron_CEP_PDE)$Cell.type <- as.character(colData(neuron_CEP_PDE)$cluster)
colData(neuron_CEP_PDE)$Cell.type <- dplyr::recode(
  colData(neuron_CEP_PDE)$Cell.type,
  "1" = "CEP",
  "2" = "CEP"
)

cd_new  <- colData(neuron_CEP_PDE)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))


#ASG_AWA
neuron_ASG_AWA <- final_neuron[, colData(final_neuron)$Cell.type == "ASG_AWA"] 

neuron_ASG_AWA <- detect_genes(neuron_ASG_AWA)
neuron_ASG_AWA <- neuron_ASG_AWA[rowData(neuron_ASG_AWA)$num_cells_expressed > 5, ]
neuron_ASG_AWA

neuron_ASG_AWA <- preprocess_cds(neuron_ASG_AWA, num_dim = 20)
plot_pc_variance_explained(neuron_ASG_AWA)

neuron_ASG_AWA <- reduce_dimension(neuron_ASG_AWA,
                                   reduction_method = "UMAP",
                                   preprocess_method = "PCA",
                                   umap.min_dist = 0.05,
                                   umap.n_neighbors = 15)

neuron_ASG_AWA <- cluster_cells(neuron_ASG_AWA, reduction_method = "UMAP", res = 20e-3)
table(monocle3::clusters(neuron_ASG_AWA))

plot_cells(neuron_ASG_AWA, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_ASG_AWA, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_ASG_AWA)$cluster <- monocle3::clusters(neuron_ASG_AWA)

markers <- top_markers(neuron_ASG_AWA, group_cells_by = "cluster")
markers %>% filter(cell_group == 3) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_ASG_AWA)$cluster

metadata <- colData(neuron_ASG_AWA) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 1) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_ASG_AWA)$UMAP_1 <- reducedDims(neuron_ASG_AWA)[["UMAP"]][,1]
colData(neuron_ASG_AWA)$UMAP_2 <- reducedDims(neuron_ASG_AWA)[["UMAP"]][,2]
plot.expr.UMAP(neuron_ASG_AWA, "sbt-1", size = 2)

colData(neuron_ASG_AWA)$Cell.type <- as.character(colData(neuron_ASG_AWA)$cluster)
colData(neuron_ASG_AWA)$Cell.type <- dplyr::recode(
  colData(neuron_ASG_AWA)$Cell.type,
  "1" = "Unannotated",
  "2" = "Unannotated Parental",
  "3" = "Unannotated Parental"
)

cd_new  <- colData(neuron_ASG_AWA)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))


#AS*_Cell
neuron_AS__Cell <- final_neuron[, colData(final_neuron)$Cell.type == "AS*_Cell"]

neuron_AS__Cell <- detect_genes(neuron_AS__Cell)
neuron_AS__Cell <- neuron_AS__Cell[rowData(neuron_AS__Cell)$num_cells_expressed > 5, ]
neuron_AS__Cell

neuron_AS__Cell <- preprocess_cds(neuron_AS__Cell, num_dim = 20)
plot_pc_variance_explained(neuron_AS__Cell)

neuron_AS__Cell <- reduce_dimension(neuron_AS__Cell,
                                   reduction_method = "UMAP",
                                   preprocess_method = "PCA",
                                   umap.min_dist = 0.05,
                                   umap.n_neighbors = 15)

neuron_AS__Cell <- cluster_cells(neuron_AS__Cell, reduction_method = "UMAP", res = 20e-3)
table(monocle3::clusters(neuron_AS__Cell))

plot_cells(neuron_AS__Cell, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_AS__Cell, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_AS__Cell)$cluster <- monocle3::clusters(neuron_AS__Cell)

markers <- top_markers(neuron_AS__Cell, group_cells_by = "cluster")
markers %>% filter(cell_group == 1) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_AS__Cell)$cluster

metadata <- colData(neuron_AS__Cell) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 1) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_AS__Cell)$UMAP_1 <- reducedDims(neuron_AS__Cell)[["UMAP"]][,1]
colData(neuron_AS__Cell)$UMAP_2 <- reducedDims(neuron_AS__Cell)[["UMAP"]][,2]
plot.expr.UMAP(neuron_AS__Cell, "sbt-1", size = 2)

colData(neuron_AS__Cell)$Cell.type <- as.character(colData(neuron_AS__Cell)$cluster)
colData(neuron_AS__Cell)$Cell.type <- dplyr::recode(
  colData(neuron_AS__Cell)$Cell.type,
  "1" = "A*_Cell",
  "2" = "ASH",
  "3" = "AWB",
  "4" = "ASE"
)

cd_new  <- colData(neuron_AS__Cell)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))


#ALM_PLM
neuron_ALM_PLM <- final_neuron[, colData(final_neuron)$Cell.type == "ALM_PLM"] 

neuron_ALM_PLM <- detect_genes(neuron_ALM_PLM)
neuron_ALM_PLM <- neuron_ALM_PLM[rowData(neuron_ALM_PLM)$num_cells_expressed > 5, ]
neuron_ALM_PLM

neuron_ALM_PLM <- preprocess_cds(neuron_ALM_PLM, num_dim = 20)
plot_pc_variance_explained(neuron_ALM_PLM)

neuron_ALM_PLM <- reduce_dimension(neuron_ALM_PLM,
                                    reduction_method = "UMAP",
                                    preprocess_method = "PCA",
                                    umap.min_dist = 0.05,
                                    umap.n_neighbors = 15)

neuron_ALM_PLM <- cluster_cells(neuron_ALM_PLM, reduction_method = "UMAP", res = 20e-3)
table(monocle3::clusters(neuron_ALM_PLM))

plot_cells(neuron_ALM_PLM, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_ALM_PLM, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_ALM_PLM)$cluster <- monocle3::clusters(neuron_ALM_PLM)

markers <- top_markers(neuron_ALM_PLM, group_cells_by = "cluster")
markers %>% filter(cell_group == 1) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_ALM_PLM)$cluster

metadata <- colData(neuron_ALM_PLM) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 1) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_ALM_PLM)$UMAP_1 <- reducedDims(neuron_ALM_PLM)[["UMAP"]][,1]
colData(neuron_ALM_PLM)$UMAP_2 <- reducedDims(neuron_ALM_PLM)[["UMAP"]][,2]
plot.expr.UMAP(neuron_ALM_PLM, "mec-17", size = 2)

colData(neuron_ALM_PLM)$Cell.type <- as.character(colData(neuron_ALM_PLM)$cluster)
colData(neuron_ALM_PLM)$Cell.type <- dplyr::recode(
  colData(neuron_ALM_PLM)$Cell.type,
  "1" = "Unannotatead",
  "2" = "Unannotated Parental",
  "3" = "ALM_PLM",
  "4" = "ALM_PLM",
  "5" = "PLM_ALN_parental"
)

cd_new  <- colData(neuron_ALM_PLM)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))

#ADF_AWB
neuron_ADF_AWB <- final_neuron[, colData(final_neuron)$Cell.type == "ADF_AWB"]

neuron_ADF_AWB <- detect_genes(neuron_ADF_AWB)
neuron_ADF_AWB <- neuron_ADF_AWB[rowData(neuron_ADF_AWB)$num_cells_expressed > 5, ]
neuron_ADF_AWB

neuron_ADF_AWB <- preprocess_cds(neuron_ADF_AWB, num_dim = 20)
plot_pc_variance_explained(neuron_ADF_AWB)

neuron_ADF_AWB <- reduce_dimension(neuron_ADF_AWB,
                                   reduction_method = "UMAP",
                                   preprocess_method = "PCA",
                                   umap.min_dist = 0.05,
                                   umap.n_neighbors = 15)

neuron_ADF_AWB <- cluster_cells(neuron_ADF_AWB, reduction_method = "UMAP", res = 90e-3)
table(monocle3::clusters(neuron_ADF_AWB))

plot_cells(neuron_ADF_AWB, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_ADF_AWB, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_ADF_AWB)$cluster <- monocle3::clusters(neuron_ADF_AWB)

markers <- top_markers(neuron_ADF_AWB, group_cells_by = "cluster")
markers %>% filter(cell_group == 3) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_ADF_AWB)$cluster

metadata <- colData(neuron_ADF_AWB) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 2) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_ADF_AWB)$UMAP_1 <- reducedDims(neuron_ADF_AWB)[["UMAP"]][,1]
colData(neuron_ADF_AWB)$UMAP_2 <- reducedDims(neuron_ADF_AWB)[["UMAP"]][,2]
plot.expr.UMAP(neuron_ADF_AWB, "C13G3.1", size = 2)

colData(neuron_ADF_AWB)$Cell.type <- as.character(colData(neuron_ADF_AWB)$cluster)
colData(neuron_ADF_AWB)$Cell.type <- dplyr::recode(
  colData(neuron_ADF_AWB)$Cell.type,
  "1" = "AWB",
  "2" = "Unannotated",
  "3" = "Unannotated"
)

cd_new  <- colData(neuron_ADF_AWB)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))


#A*_Cell 
neuron_A_Cell <- final_neuron[, colData(final_neuron)$Cell.type == "A*_Cell"]

neuron_A_Cell <- detect_genes(neuron_A_Cell)
neuron_A_Cell <- neuron_A_Cell[rowData(neuron_A_Cell)$num_cells_expressed > 5, ]
neuron_A_Cell

neuron_A_Cell <- preprocess_cds(neuron_A_Cell, num_dim = 20)
plot_pc_variance_explained(neuron_A_Cell)

neuron_A_Cell <- reduce_dimension(neuron_A_Cell,
                                   reduction_method = "UMAP",
                                   preprocess_method = "PCA",
                                   umap.min_dist = 0.05,
                                   umap.n_neighbors = 15)

neuron_A_Cell <- cluster_cells(neuron_A_Cell, reduction_method = "UMAP", res = 20e-3)
table(monocle3::clusters(neuron_A_Cell))

plot_cells(neuron_A_Cell, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_A_Cell, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_A_Cell)$cluster <- monocle3::clusters(neuron_A_Cell)

markers <- top_markers(neuron_A_Cell, group_cells_by = "cluster")
markers %>% filter(cell_group == 4) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_A_Cell)$cluster

metadata <- colData(neuron_A_Cell) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 3) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_A_Cell)$UMAP_1 <- reducedDims(neuron_A_Cell)[["UMAP"]][,1]
colData(neuron_A_Cell)$UMAP_2 <- reducedDims(neuron_A_Cell)[["UMAP"]][,2]
plot.expr.UMAP(neuron_A_Cell, "droe-4", size = 2)

colData(neuron_A_Cell)$Cell.type <- as.character(colData(neuron_A_Cell)$cluster)
colData(neuron_A_Cell)$Cell.type <- dplyr::recode(
  colData(neuron_A_Cell)$Cell.type,
  "1" = "IL2",
  "2" = "AWC",
  "3" = "AWC",
  "4" = "IL2",
  "5" = "ASJ",
  "6" = "ASJ"
)

cd_new  <- colData(neuron_A_Cell)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))



#AV*2
neuron_AV__2 <- final_neuron[, colData(final_neuron)$Cell.type == "AV*"]

neuron_AV__2 <- detect_genes(neuron_AV__2)
neuron_AV__2 <- neuron_AV__2[rowData(neuron_AV__2)$num_cells_expressed > 5, ]
neuron_AV__2

neuron_AV__2 <- preprocess_cds(neuron_AV__2, num_dim = 20)
plot_pc_variance_explained(neuron_AV__2)

neuron_AV__2 <- reduce_dimension(neuron_AV__2,
                                 reduction_method = "UMAP",
                                 preprocess_method = "PCA",
                                 umap.min_dist = 0.05,
                                 umap.n_neighbors = 15)

neuron_AV__2 <- cluster_cells(neuron_AV__2, reduction_method = "UMAP", res = 90e-3)
table(monocle3::clusters(neuron_AV__2))

plot_cells(neuron_AV__2, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = 1,
           label_groups_by_cluster = F)

plot_cells(neuron_AV__2, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 1)

colData(neuron_AV__2)$cluster <- monocle3::clusters(neuron_AV__2)

markers <- top_markers(neuron_AV__2, group_cells_by = "cluster")
markers %>% filter(cell_group == 1) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_AV__2)$cluster

metadata <- colData(neuron_AV__2) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 3) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_AV__2)$UMAP_1 <- reducedDims(neuron_AV__2)[["UMAP"]][,1]
colData(neuron_AV__2)$UMAP_2 <- reducedDims(neuron_AV__2)[["UMAP"]][,2]
plot.expr.UMAP(neuron_AV__2, "lite-1", size = 2)

colData(neuron_AV__2)$Cell.type <- as.character(colData(neuron_AV__2)$cluster)
colData(neuron_AV__2)$Cell.type <- dplyr::recode(
  colData(neuron_AV__2)$Cell.type,
  "1" = "Unannotated",
  "2" = "Unannotated",
  "3" = "Unannotated",
  "4" = "AVJ"
)

cd_new  <- colData(neuron_AV__2)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))

#Final Unannotated
neuron_final_unannotated <- final_neuron[, colData(final_neuron)$Cell.type == "Unannotated"]

neuron_final_unannotated <- detect_genes(neuron_final_unannotated)
neuron_final_unannotated <- neuron_final_unannotated[rowData(neuron_final_unannotated)$num_cells_expressed > 5, ]
neuron_final_unannotated

neuron_final_unannotated <- preprocess_cds(neuron_final_unannotated, num_dim = 50)
plot_pc_variance_explained(neuron_final_unannotated)

neuron_final_unannotated <- reduce_dimension(neuron_final_unannotated,
                                  reduction_method = "UMAP",
                                  preprocess_method = "PCA",
                                  umap.min_dist = 0.2,
                                  umap.n_neighbors = 40)

neuron_final_unannotated <- cluster_cells(neuron_final_unannotated, reduction_method = "UMAP", res = 3e-3)
table(monocle3::clusters(neuron_final_unannotated))

plot_cells(neuron_final_unannotated, 
           color_cells_by = "cell.subtype", #coloring it by the cell.subtype col
           group_label_size = 3, #how big text should be
           cell_size = .3,
           label_groups_by_cluster = F)

plot_cells(neuron_final_unannotated, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = .3)

colData(neuron_final_unannotated)$cluster <- monocle3::clusters(neuron_final_unannotated)

markers <- top_markers(neuron_final_unannotated, group_cells_by = "cluster")
markers %>% filter(cell_group == 46) %>% arrange(desc(marker_score)) %>% head(20)

colData(neuron_final_unannotated)$cluster

metadata <- colData(neuron_final_unannotated) %>% as.data.frame()
# Count the number of cells for each annotation
metadata %>%
  filter(cluster == 30) %>%
  count(cell.subtype, sort = TRUE) 

# Displaying egl-21 on UMAP (highlights neuronal cells we want)
colData(neuron_final_unannotated)$UMAP_1 <- reducedDims(neuron_final_unannotated)[["UMAP"]][,1]
colData(neuron_final_unannotated)$UMAP_2 <- reducedDims(neuron_final_unannotated)[["UMAP"]][,2]
plot.expr.UMAP(neuron_final_unannotated, "nlp-50", size = 2)

colData(neuron_final_unannotated)$Cell.type <- as.character(colData(neuron_final_unannotated)$cluster)
colData(neuron_final_unannotated)$Cell.type <- dplyr::recode(
  colData(neuron_final_unannotated)$Cell.type,
  "1" = "Unannotated",
  "2" = "Unannotated",
  "3" = "Unannotated",
  "4" = "Unannotated",
  "5" = "Unannotated",
  "6" = "Unannotated",
  "7" = "Unannotated",
  "8" = "Unannotated",
  "9" = "Unannotated",
  "10" = "Unannotated",
  "11" = "Unannotated",
  "12" = "Unannotated",
  "13" = "Unannotated",
  "14" = "Unannotated",
  "15" = "Unannotated",
  "16" = "Unannotated",
  "17" = "Unannotated",
  "18" = "Unannotated",
  "19" = "Unannotated",
  "20" = "Unannotated",
  "21" = "Unannotated",
  "22" = "Unannotated",
  "23" = "Unannotated",
  "24" = "Unannotated",
  "25" = "Unannotated",
  "26" = "Unannotated",
  "27" = "Unannotated",
  "28" = "Unannotated",
  "29" = "Unannotated",
  "30" = "Unannotated",
  "31" = "Unannotated",
  "32" = "Unannotated",
  "33" = "Unannotated",
  "34" = "ALM_PLM",
  "35" = "Unannotated",
  "36" = "Unannotated",
  "37" = "Unannotated",
  "38" = "Unannotated",
  "39" = "AVJ",
  "40" = "Unannotated",
  "41" = "Unannotated",
  "42" = "Unannotated",
  "43" = "Unannotated",
  "44" = "Unannotated",
  "45" = "Unannotated"
)

cd_new  <- colData(neuron_final_unannotated)
cd_final    <- colData(final_neuron)

shared_cells <- intersect(rownames(cd_new), rownames(cd_final))
cd_final[shared_cells, "Cell.type"] <- cd_new[shared_cells, "Cell.type"]
colData(final_neuron) <- cd_final

table(colData(final_neuron)$Cell.type)
(colData(final_neuron))

saveRDS(final_neuron, "~/Taylor_Lab/final_neuron_embryo_annotations.rds")
#========================================================================================
saveRDS(neuron, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/042324_L1_all_live_cells_MeOH_new_saa_P_lineage_no_doublets_or_UPR_cds.rds")

# Separating neurons from non-neuronal cells

plot.expr.UMAP(neuron, "sbt-1", size = 0.5)

plot_cells(neuron, color_cells_by = "PCA.cluster", 
           group_label_size = 4, cell_size = 0.5)

neuron <- c(51,73,57,43,69,45,50,39,64,49,60,63,47,30,6,35,23,61,78,38,74,40,2,33,3,56,19,8,16,19,80,
            55,24,48,44,59,18,31,22,65,76,76,52,77,53,75,54,12,29,10,5,58,34,25,9,62,4,26,42,11,32,27,
            21,41,61,68, 72)

ggplot(as.data.frame(colData(neuron)), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = ifelse(PCA.cluster %in% neuron, "red", "black")), show.legend = "NA") + 
  scale_color_manual(values = c("black", "red")) 

colData(neuron)$tissue <- ifelse(
  colData(neuron)$PCA.cluster %in% neuron,
  "Neuron",
  "Nonneuronal"
)

L1.neuron <- neuron[,colData(neuron)$PCA.cluster %in% neuron]

L1.neuron <- detect_genes(L1.neuron)
L1.neuron <- L1.neuron[rowData(L1.neuron)$num_cells_expressed > 5,]
L1.neuron
# 17758 features in 125504 cells

L1.neuron <- preprocess_cds(L1.neuron, num_dim = 60)
plot_pc_variance_explained(L1.neuron)

L1.neuron <- align_cds(L1.neuron, aligment_group = "Experiment", alignment_k = 5)
L1.neuron <- reduce_dimension(L1.neuron,
                                reduction_method = "UMAP",
                                preprocess_method = "Aligned",
                                umap.min_dist = 0.3,
                                umap.n_neighbors = 75)

L1.neuron <- cluster_cells(L1.neuron, reduction_method = "PCA", res = 3e-4)

plot_cells(L1.neuron, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.neuron, 
           color_cells_by = "Experiment",
           label_cell_groups = F,
           cell_size = 0.5)

plot_cells(L1.neuron, 
           color_cells_by = "hph",
           label_cell_groups = F,
           cell_size = 0.5)

colData(L1.neuron)$PCA.cluster <- monocle3::clusters(L1.neuron, reduction_method = "PCA")
colData(L1.neuron)$PCA.partition <- monocle3::partitions(L1.neuron, reduction_method = "PCA")

plot_cells(L1.neuron,
           color_cells_by = "PCA.partition",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

plot_cells(L1.neuron,
           color_cells_by = "PCA.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(L1.neuron)$UMAP_1 <- reducedDims(L1.neuron)[["UMAP"]][,1]
colData(L1.neuron)$UMAP_2 <- reducedDims(L1.neuron)[["UMAP"]][,2]

plot.expr.UMAP(L1.neuron, "sbt-1", size = 0.5)

# Cleaning up some annotations

plot.cell.type.m3(L1.neuron, "Unknown_neuron")
plot.cell.type.m3(L1.neuron, "Unknown")
plot.cell.type.m3(L1.neuron, "Unannotated")

colData(L1.neuron)$Cell.type <- ifelse(
  colData(L1.neuron)$Cell.type == "Unknown_neuron",
  as.character(colData(L1.neuron)$PCA.cluster),
  as.character(colData(L1.neuron)$Cell.type)
)

table(colData(L1.neuron)$Cell.type)

colData(L1.neuron)$Cell.type <- dplyr::recode(colData(L1.neuron)$Cell.type,
                                              "55" = "BAG",
                                              "65" = "AFD",
                                              "63" = "ASI",
                                              "34" = "RIC",
                                              "39" = "RIM",
                                              "38" = "ADL",
                                              "25" = "ASK",
                                              "35" = "ASJ",
                                              "PHso" = "ASJ",
                                              "52" = "ADF",
                                              "60" = "AWA",
                                              "31" = "ASG",
                                              "64" = "AVH",
                                              "5" = "RIA",
                                              "10" = "AIN",
                                              "17" = "AIY",
                                              "16" = "AVJ",
                                              "47" = "RIF",
                                              "62" = "URB",
                                              "67" = "hmc",
                                              "22" = "AVG",
                                              "36" = "NSM",
                                              "46" = "AVK",
                                              "49" = "AUA",
                                              "66" = "AIA",
                                              "21" = "AIB",
                                              "4" = "AIZ",
                                              "44" = "ALA")

plot_cells(L1.neuron, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.neuron,
           color_cells_by = "PCA.cluster",
           group_label_size = 4,
           cell_size = 0.5)

L1.sub.out <- L1.neuron[,colData(L1.neuron)$PCA.cluster %in% c(56,30,48,40,33,69,71,27,28,41,45,
                                                               53)]

table(colData(L1.neuron)$Cell.type)
L1.sub.out <- detect_genes(L1.sub.out)
L1.sub.out <- L1.sub.out[rowData(L1.sub.out)$num_cells_expressed > 5,]
L1.sub.out
# 12021 features in 10653 cells

L1.sub.out <- preprocess_cds(L1.sub.out, num_dim = 40)
plot_pc_variance_explained(L1.sub.out)

L1.sub.out <- align_cds(L1.sub.out, aligment_group = "Experiment", alignment_k = 5)
L1.sub.out <- reduce_dimension(L1.sub.out,
                              reduction_method = "UMAP",
                              preprocess_method = "Aligned",
                              umap.min_dist = 0.2,
                              umap.n_neighbors = 50)

L1.sub.out <- cluster_cells(L1.sub.out, reduction_method = "PCA", res = 3e-3)

plot_cells(L1.sub.out, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(L1.sub.out)$PCA.cluster <- monocle3::clusters(L1.sub.out, reduction_method = "PCA")

plot_cells(L1.sub.out,
           color_cells_by = "PCA.cluster",
           group_label_size = 4,
           cell_size = 0.5)

table(colData(L1.sub.out)$PCA.cluster,colData(L1.sub.out)$Cell.type, exclude = NULL)

colData(L1.sub.out)$UMAP_1 <- reducedDims(L1.sub.out)[["UMAP"]][,1]
colData(L1.sub.out)$UMAP_2 <- reducedDims(L1.sub.out)[["UMAP"]][,2]

plot.expr.UMAP(L1.sub.out, "sbt-1", size = 0.5)

# Cleaning up some annotations

colData(L1.sub.out)$Cell.type <- ifelse(
  colData(L1.sub.out)$PCA.cluster %in% c(2,22),
  "PQR",
  ifelse(
    colData(L1.sub.out)$PCA.cluster %in% c(5,17),
    "AQR",
    ifelse(
      colData(L1.sub.out)$PCA.cluster %in% c(13),
      "URX",
      ifelse(colData(L1.sub.out)$PCA.cluster %in% c(15),
             "AWB",
             ifelse(
               colData(L1.sub.out)$PCA.cluster %in% (19),
               "IL2_LR",
               ifelse(colData(L1.sub.out)$PCA.cluster %in% c(10),
                      "IL2_DV",
                      ifelse(
                        colData(L1.sub.out)$PCA.cluster %in% c(11),
                        "ASH",
                        ifelse(colData(L1.sub.out)$PCA.cluster %in% c(12),
                               "PHA",
                               ifelse(
                                 colData(L1.sub.out)$PCA.cluster %in% c(14),
                                 "PHB",
                                 as.character(colData(L1.sub.out)$Cell.type)
                               ))
                      )
                      )
             )
             )
    )
  )
) 

plot_cells(L1.sub.out, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.sub.out,
           color_cells_by = "PCA.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(L1.sub.out)$Cell.type <- ifelse(
  colData(L1.sub.out)$PCA.cluster %in% c(3,8),
  "RIG",
  ifelse(
    colData(L1.sub.out)$PCA.cluster %in% c(7),
    "I5",
    ifelse(
      colData(L1.sub.out)$PCA.cluster %in% c(6),
      "ALM",
      ifelse(colData(L1.sub.out)$PCA.cluster %in% c(20),
             "PLM",
             as.character(colData(L1.sub.out)$Cell.type)
                               ))
                      )
               )

plot_cells(L1.sub.out, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.sub.out,
           color_cells_by = "PCA.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(L1.sub.out)$Cell.type <- ifelse(
  colData(L1.sub.out)$PCA.cluster %in% c(16),
  "ASEL",
  ifelse(
    colData(L1.sub.out)$PCA.cluster %in% c(4,21),
    "ASER",
    ifelse(
      colData(L1.sub.out)$PCA.cluster %in% c(18),
      "PVT",
      ifelse(colData(L1.sub.out)$PCA.cluster %in% c(25),
             "FLP",
             as.character(colData(L1.sub.out)$Cell.type)
      ))
  )
)

colData(L1.sub.out)$Cell.type <- ifelse(
  colData(L1.sub.out)$PCA.cluster %in% c(26),
  "DVA",
  ifelse(
    colData(L1.sub.out)$PCA.cluster %in% c(24),
    "RMG",
    as.character(colData(L1.sub.out)$Cell.type)
  )
)

plot.cell.type.m3(L1.sub.out, "Unknown_neuron")
plot.cell.type.m3(L1.sub.out, "Unknown")
plot.cell.type.m3(L1.sub.out, "Unannotated")

table(colData(L1.sub.out)$Cell.type)

plot_cells(L1.sub.out, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.sub.out,
           color_cells_by = "PCA.cluster",
           group_label_size = 4,
           cell_size = 0.5)
L1.sub.out2 <- L1.sub.out[,colData(L1.sub.out)$PCA.cluster %in% c(1,23,20,6,9,10, 19)]
L1.sub.out2 <- detect_genes(L1.sub.out2)
L1.sub.out2 <- L1.sub.out2[rowData(L1.sub.out2)$num_cells_expressed > 5,]
L1.sub.out2
# 12021 features in 10653 cells

L1.sub.out2 <- preprocess_cds(L1.sub.out2, num_dim = 25)
plot_pc_variance_explained(L1.sub.out2)

L1.sub.out2 <- align_cds(L1.sub.out2, aligment_group = "Experiment", alignment_k = 5)
L1.sub.out2 <- reduce_dimension(L1.sub.out2,
                               reduction_method = "UMAP",
                               preprocess_method = "Aligned",
                               umap.min_dist = 0.1,
                               umap.n_neighbors = 25)

L1.sub.out2 <- cluster_cells(L1.sub.out2, reduction_method = "PCA", res = 8e-3)

plot_cells(L1.sub.out2, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(L1.sub.out2)$PCA.cluster <- monocle3::clusters(L1.sub.out2, reduction_method = "PCA")

plot_cells(L1.sub.out2,
           color_cells_by = "PCA.cluster",
           group_label_size = 4,
           cell_size = 0.5)

table(colData(L1.sub.out2)$PCA.cluster,colData(L1.sub.out2)$Cell.type, exclude = NULL)

colData(L1.sub.out2)$UMAP_1 <- reducedDims(L1.sub.out2)[["UMAP"]][,1]
colData(L1.sub.out2)$UMAP_2 <- reducedDims(L1.sub.out2)[["UMAP"]][,2]

plot.expr.UMAP(L1.sub.out2, "lin-39", size = 0.5)
plot.expr.UMAP(L1.sub.out2, "mab-5", size = 0.5)
plot.expr.UMAP(L1.sub.out2, "unc-86", size = 0.5)

colData(L1.sub.out2)$Cell.type <- ifelse(
  colData(L1.sub.out2)$PCA.cluster %in% c(12),
  "AVM",
  as.character(colData(L1.sub.out2)$Cell.type)
)

colData(L1.sub.out2)$Cell.type <- ifelse(
  colData(L1.sub.out2)$PCA.cluster %in% c(9),
  "AWC_OFF",
  ifelse(
    colData(L1.sub.out2)$PCA.cluster %in% (2),
    "AWC_ON",
  as.character(colData(L1.sub.out2)$Cell.type)
  )
)

plot.expr.UMAP(L1.sub.out2, "nlp-80", size = 0.5)
plot.expr.UMAP(L1.sub.out2, "nlp-70", size = 0.5)
plot.expr.UMAP(L1.sub.out2, "snet-1", size = 0.5)

sub2.mark <- top_markers(L1.sub.out2, group_cells_by = "PCA.cluster")

sub2.mark %>% filter(cell_group == 8) %>% arrange(desc(specificity)) %>% head(20)

colData(L1.sub.out2)$Cell.type <- ifelse(
  colData(L1.sub.out2)$PCA.cluster %in% c(4,5,7),
  "AIM",
  ifelse(
    colData(L1.sub.out2)$PCA.cluster %in% (8),
    "Doublets",
    as.character(colData(L1.sub.out2)$Cell.type)
  )
)

table(colData(L1.sub.out2)$Cell.type, colData(L1.sub.out2)$Experiment, exclude = NULL)

colData(L1.sub.out)[colnames(L1.sub.out2),]$Cell.type <- colData(L1.sub.out2)$Cell.type
colData(L1.neuron)[colnames(L1.sub.out),]$Cell.type <- colData(L1.sub.out)$Cell.type
colData(neuron)[colnames(L1.neuron),]$Cell.type <- colData(L1.neuron)$Cell.type
colData(combined_cds)[colnames(neuron),]$Cell.type <- colData(neuron)$Cell.type

plot_cells(L1.sub.out, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.neuron, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.neuron,
           color_cells_by = "PCA.cluster",
           group_label_size = 4,
           cell_size = 0.5)

plot.cell.type.m3(L1.neuron, "RIP")

colData(L1.neuron)$Cell.type <- ifelse(
  colData(L1.neuron)$PCA.cluster %in% c(5),
  "RIA",
    as.character(colData(L1.neuron)$Cell.type)
  )

sub3 <- c(29,4,71,44,59,70,59,68,54,2,9,13,32,43,23,61,29,20,37,36,51, 22, 15, 57)

ggplot(as.data.frame(colData(L1.neuron)), aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(aes(color = ifelse(PCA.cluster %in% sub3, "red", "black")),show.legend = "NA") + 
  scale_color_manual(values = c("red", "black"))

L1.sub3 <- L1.neuron[,colData(L1.neuron)$PCA.cluster %in% sub3]

L1.sub3 <- detect_genes(L1.sub3)
L1.sub3 <- L1.sub3[rowData(L1.sub3)$num_cells_expressed > 5,]
L1.sub3
# 14019 features in 39014 cells

L1.sub3 <- preprocess_cds(L1.sub3, num_dim = 55)
plot_pc_variance_explained(L1.sub3)

L1.sub3 <- align_cds(L1.sub3, aligment_group = "Experiment", alignment_k = 5)
L1.sub3 <- reduce_dimension(L1.sub3,
                                reduction_method = "UMAP",
                                preprocess_method = "Aligned",
                                umap.min_dist = 0.1,
                                umap.n_neighbors = 25)

L1.sub3 <- cluster_cells(L1.sub3, reduction_method = "UMAP", res = 3e-4)

plot_cells(L1.sub3, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(L1.sub3)$UMAP.cluster <- monocle3::clusters(L1.sub3, reduction_method = "UMAP")

plot_cells(L1.sub3,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(L1.sub3)$UMAP_1 <- reducedDims(L1.sub3)[["UMAP"]][,1]
colData(L1.sub3)$UMAP_2 <- reducedDims(L1.sub3)[["UMAP"]][,2]

sub3.mark <- top_markers(L1.sub3, group_cells_by = "UMAP.cluster", speedglm.maxiter = 100)

sub3.mark %>% filter(cell_group == 60) %>% arrange(desc(specificity)) %>% head(20)

# 62 looks like AIN, could be from the AIN cluster, but 63 and 65 look like AVG.

plot.expr.UMAP(L1.sub3, "sbt-1", size = 0.5)
plot.expr.UMAP(L1.sub3, "tab-1", size = 0.5)
plot.expr.UMAP(L1.sub3, "egl-5", size = 0.5)

colData(L1.sub3)$Cell.type <- ifelse(
  colData(L1.sub3)$UMAP.cluster %in% c(55),
  "PVR",
  ifelse(
    colData(L1.sub3)$UMAP.cluster %in% c(47),
    "LUA",
    ifelse(
      colData(L1.sub3)$UMAP.cluster %in% c(10),
      "MC",
      ifelse(colData(L1.sub3)$UMAP.cluster %in% c(1,3,12,32),
             "AIZ",
             ifelse(
               colData(L1.sub3)$UMAP.cluster %in% (53),
               "PVP",
               ifelse(colData(L1.sub3)$UMAP.cluster %in% c(52),
                      "RIR",
                      ifelse(
                        colData(L1.sub3)$UMAP.cluster %in% c(30,4),
                        "SMD",
                        ifelse(colData(L1.sub3)$UMAP.cluster %in% c(13),
                               "RMD_DV",
                               ifelse(
                                 colData(L1.sub3)$UMAP.cluster %in% c(31),
                                 "AVE",
                                 as.character(colData(L1.sub3)$Cell.type)
                               ))
                      )
               )
             )
      )
    )
  )
) 

colData(L1.sub3)$Cell.type <- ifelse(
  colData(L1.sub3)$UMAP.cluster %in% c(22),
  "RMD_LR",
  ifelse(
    colData(L1.sub3)$UMAP.cluster %in% c(42),
    "RIV",
    ifelse(
      colData(L1.sub3)$UMAP.cluster %in% c(35),
      "BDU",
      ifelse(colData(L1.sub3)$UMAP.cluster %in% c(41),
             "PVQ",
             ifelse(
               colData(L1.sub3)$UMAP.cluster %in% (15),
               "ADA",
               ifelse(colData(L1.sub3)$UMAP.cluster %in% c(48),
                      "CAN",
                      ifelse(
                        colData(L1.sub3)$UMAP.cluster %in% c(34),
                        "M2",
                        ifelse(colData(L1.sub3)$UMAP.cluster %in% c(37),
                               "RIH",
                               ifelse(
                                 colData(L1.sub3)$UMAP.cluster %in% c(33),
                                 "RME_DV",
                                 as.character(colData(L1.sub3)$Cell.type)
                               ))
                      )
               )
             )
      )
    )
  )
) 

colData(L1.sub3)$Cell.type <- ifelse(
  colData(L1.sub3)$UMAP.cluster %in% c(57),
  "RIP",
  ifelse(
    colData(L1.sub3)$UMAP.cluster %in% c(44),
    "RME_LR",
    ifelse(
      colData(L1.sub3)$UMAP.cluster %in% c(45),
      "DVC",
      ifelse(colData(L1.sub3)$UMAP.cluster %in% c(43),
             "AVL",
             ifelse(
               colData(L1.sub3)$UMAP.cluster %in% (50),
               "I4",
               ifelse(colData(L1.sub3)$UMAP.cluster %in% c(9,11),
                      "URA",
                      ifelse(
                        colData(L1.sub3)$UMAP.cluster %in% c(49),
                        "M1",
                        ifelse(colData(L1.sub3)$UMAP.cluster %in% c(64),
                               "I6",
                               ifelse(
                                 colData(L1.sub3)$UMAP.cluster %in% c(54),
                                 "I4",
                                 as.character(colData(L1.sub3)$Cell.type)
                               ))
                      )
               )
             )
      )
    )
  )
) 

plot_cells(L1.sub3, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.sub3,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(L1.neuron)[colnames(L1.sub3),]$Cell.type <- colData(L1.sub3)$Cell.type
colData(neuron)[colnames(L1.neuron),]$Cell.type <- colData(L1.neuron)$Cell.type
colData(combined_cds)[colnames(neuron),]$Cell.type <- colData(neuron)$Cell.type

saveRDS(L1.neuron, "~/Dropbox (VU Basic Sciences)/L1 Cengen manuscript/Claire data and presentations/no_intron_analysis/data/042424_all_neurons_combined_cds.rds")
saveRDS(neuron, "~/Dropbox (VU Basic Sciences)/L1 Cengen manuscript/Claire data and presentations/no_intron_analysis/data/042424_all_cells_combined_filtered_cds.rds")
saveRDS(combined_cds, "~/Dropbox (VU Basic Sciences)/L1 Cengen manuscript/Claire data and presentations/no_intron_analysis/data/042424_all_cells_combined_no_filter_cds.rds")

# Subsetting a few more clusters

L1.sub3.2 <- L1.sub3[,colData(L1.sub3)$UMAP.cluster %in% c(29, 17, 36, 40, 46, 25, 39,5,20)]

L1.sub3.2 <- detect_genes(L1.sub3.2)
L1.sub3.2 <- L1.sub3.2[rowData(L1.sub3.2)$num_cells_expressed > 5,]
L1.sub3.2
# 10258 features in 5270 cells

L1.sub3.2 <- preprocess_cds(L1.sub3.2, num_dim = 45)
plot_pc_variance_explained(L1.sub3.2)

L1.sub3.2 <- align_cds(L1.sub3.2, aligment_group = "Experiment", alignment_k = 5)
L1.sub3.2 <- reduce_dimension(L1.sub3.2,
                            reduction_method = "UMAP",
                            preprocess_method = "Aligned",
                            umap.min_dist = 0.1,
                            umap.n_neighbors = 25)

L1.sub3.2 <- cluster_cells(L1.sub3.2, reduction_method = "UMAP", res = 6e-3)

plot_cells(L1.sub3.2, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.sub3.2, 
           color_cells_by = "Experiment",
           cell_size = 0.5,
           label_cell_groups = F)

colData(L1.sub3.2)$UMAP.cluster <- monocle3::clusters(L1.sub3.2, reduction_method = "UMAP")

plot_cells(L1.sub3.2,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(L1.sub3.2)$UMAP_1 <- reducedDims(L1.sub3.2)[["UMAP"]][,1]
colData(L1.sub3.2)$UMAP_2 <- reducedDims(L1.sub3.2)[["UMAP"]][,2]

sub3.mark2 <- top_markers(L1.sub3.2, group_cells_by = "UMAP.cluster", marker_sig_test = F)
sub3.mark2$gene_short_name <- i2s(sub3.mark2$gene_id, gids)
sub3.mark2 %>% filter(cell_group == 3) %>% arrange(desc(specificity)) %>% head(20)


plot.expr.UMAP(L1.sub3.2, "tab-1", size = 0.5)
plot.expr.UMAP(L1.sub3.2, "egl-5", size = 0.5)

colData(L1.sub3.2)$Cell.type <- ifelse(
  colData(L1.sub3.2)$UMAP.cluster %in% c(3,25),
  "URY",
  ifelse(
    colData(L1.sub3.2)$UMAP.cluster %in% c(9,11),
    "SIB",
    ifelse(
      colData(L1.sub3.2)$UMAP.cluster %in% c(2,31),
      "SIA",
      ifelse(colData(L1.sub3.2)$UMAP.cluster %in% c(32),
             "SDQ",
             ifelse(
               colData(L1.sub3.2)$UMAP.cluster %in% (29),
               "ALN",
               ifelse(colData(L1.sub3.2)$UMAP.cluster %in% c(12),
                      "AVB",
                      ifelse(
                        colData(L1.sub3.2)$UMAP.cluster %in% c(30),
                        "SMB",
                        ifelse(colData(L1.sub3.2)$UMAP.cluster %in% c(16),
                               "SAA",
                               ifelse(
                                 colData(L1.sub3.2)$UMAP.cluster %in% c(13,21),
                                 "RID",
                                 ifelse(
                                   colData(L1.sub3.2)$UMAP.cluster %in% c(35,4),
                                   "AVA",
                                   ifelse(
                                     colData(L1.sub3.2)$UMAP.cluster %in% c(24),
                                     "AVD",
                                     ifelse(
                                       colData(L1.sub3.2)$UMAP.cluster %in% c(14),
                                       "PVC",
                                       ifelse(
                                         colData(L1.sub3.2)$UMAP.cluster %in% c(34),
                                         "OLL",
                                         ifelse(colData(L1.sub3.2)$UMAP.cluster %in% c(33),
                                                "OLQ",
                                                ifelse(colData(L1.sub3.2)$UMAP.cluster %in% c(36),
                                                       "IL1",
                                                       as.character(colData(L1.sub3.2)$Cell.type)
                                           )
                                         )
                                       )
                                     )
                                   )
                                 )
                                 
                          ))
                      )
               )
             )
      )
    )
  )
) 

plot_cells(L1.sub3.2, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.sub3.2,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

plot.expr.UMAP(L1.sub3.2, "sbt-1", size = 0.5)

colData(L1.sub3.2)$Cell.type <- ifelse(
  colData(L1.sub3.2)$UMAP.cluster %in% c(26),
  "Unannotated",
  ifelse(
    colData(L1.sub3.2)$UMAP.cluster %in% c(23),
    "Body_wall_muscle",
    ifelse(
      colData(L1.sub3.2)$UMAP.cluster %in% c(18,27),
      "M_lineage_progenitor",
      as.character(colData(L1.sub3.2)$Cell.type)
    )
  )
)

plot.cell.type.m3(L1.sub3.2, "Low-quality")

colData(L1.sub3.2)$Cell.type <- ifelse(
  colData(L1.sub3.2)$UMAP.cluster %in% c(22),
  "Excretory_gland_cell",
  ifelse(
    colData(L1.sub3.2)$UMAP.cluster %in% c(15),
    "Germline?",
    ifelse(
      colData(L1.sub3.2)$UMAP.cluster %in% c(20,17,10,5,6,8,1,19),
      "Unannotated",
      as.character(colData(L1.sub3.2)$Cell.type)
      )
    )
  )

# Cluster 15 has some germline-enriched genes (smc-5, smc-6, sco-1)

# Clusters 28 and 7 have ADE and CEP, which have been difficult to separate at this age. Both are there, I'm pretty sure

plot_genes_jitter(L1.sub3.2, "flp-33", grouping = "Cell.type")
plot_genes_jitter(L1.sub3.2, "nhr-67", grouping = "Cell.type")
plot_genes_jitter(L1.sub3.2, "nlp-11", grouping = "Cell.type")
plot_genes_jitter(L1.sub3.2, "ceh-9", grouping = "Cell.type")

sub4.DA <- L1.sub3.2[,colData(L1.sub3.2)$UMAP.cluster %in% c(28,7, 15, 17, 20)]

sub4.DA <- detect_genes(sub4.DA)
sub4.DA <- sub4.DA[rowData(sub4.DA)$num_cells_expressed > 5,]
sub4.DA
# 6574 features in 689 cells

sub4.DA <- preprocess_cds(sub4.DA, num_dim = 30)
plot_pc_variance_explained(sub4.DA)

sub4.DA <- align_cds(sub4.DA, aligment_group = "Experiment", alignment_k = 5)
sub4.DA <- reduce_dimension(sub4.DA,
                              reduction_method = "UMAP",
                              preprocess_method = "Aligned",
                              umap.min_dist = 0.1,
                              umap.n_neighbors = 25)

sub4.DA <- cluster_cells(sub4.DA, reduction_method = "UMAP", res = 6e-3)

plot_cells(sub4.DA, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(sub4.DA, 
           color_cells_by = "Experiment",
           cell_size = 0.5,
           label_cell_groups = F)

colData(sub4.DA)$UMAP.cluster <- monocle3::clusters(sub4.DA, reduction_method = "UMAP")

plot_cells(sub4.DA,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(sub4.DA)$UMAP_1 <- reducedDims(sub4.DA)[["UMAP"]][,1]
colData(sub4.DA)$UMAP_2 <- reducedDims(sub4.DA)[["UMAP"]][,2]

plot.expr.UMAP(sub4.DA, "flp-33", size = 0.5)
plot.expr.UMAP(sub4.DA, "nhr-67", size = 0.5)

colData(sub4.DA)$Cell.type <- ifelse(
  colData(sub4.DA)$UMAP.cluster %in% c(2,3) & colData(sub4.DA)$Cell.type != "ADE",
  "CEP",
    as.character(colData(sub4.DA)$Cell.type)
  )

colData(L1.sub3.2)[colnames(sub4.DA),]$Cell.type <- colData(sub4.DA)$Cell.type

colData(L1.sub3)[colnames(L1.sub3.2),]$Cell.type <- colData(L1.sub3.2)$Cell.type

ggplot(as.data.frame(colData(L1.neuron)), aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(aes(color = ifelse(PCA.cluster %in% sub3, "red", "black")),show.legend = "NA") + 
  scale_color_manual(values = c("red", "black"))

plot.cell.type.m3(L1.neuron, "AIN")
plot.cell.type.m3(L1.neuron, "AVG")

# possible AVG doublets - clusters 62 and 65. Cluster 63 - AIN

sub3.mark %>% filter(cell_group == 62) %>% arrange(desc(specificity)) %>% head(20)

colData(L1.sub3)$Cell.type <- ifelse(
  colData(L1.sub3)$UMAP.cluster %in% c(60),
  "RIB",
  ifelse(
    colData(L1.sub3)$UMAP.cluster %in% c(63),
    "AIN",
    ifelse(colData(L1.sub3)$UMAP.cluster %in% c(62,65),
           "Doublets",
    as.character(colData(L1.sub3)$Cell.type)
    )
  )
)

# Pharyngeal Neurons

plot_cells(L1.sub3, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.sub3,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

sub3.phar <- L1.sub3[,colData(L1.sub3)$UMAP.cluster %in% c(38,18,34,49,59,64,24,27,19,8,16,
                                                           28,61,50,10,21)]

sub3.phar <- detect_genes(sub3.phar)
sub3.phar <- sub3.phar[rowData(sub3.phar)$num_cells_expressed > 5,]
sub3.phar
# 9890 features in 7954 cells

sub3.phar <- preprocess_cds(sub3.phar, num_dim = 35)
plot_pc_variance_explained(sub3.phar)

sub3.phar <- align_cds(sub3.phar, aligment_group = "Experiment", alignment_k = 5)
sub3.phar <- reduce_dimension(sub3.phar,
                            reduction_method = "UMAP",
                            preprocess_method = "Aligned",
                            umap.min_dist = 0.1,
                            umap.n_neighbors = 25)

sub3.phar <- cluster_cells(sub3.phar, reduction_method = "UMAP", res = 6e-3)

plot_cells(sub3.phar, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(sub3.phar, 
           color_cells_by = "Experiment",
           cell_size = 0.5,
           label_cell_groups = F)

colData(sub3.phar)$UMAP.cluster <- monocle3::clusters(sub3.phar, reduction_method = "UMAP")

plot_cells(sub3.phar,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

phar.mark <- top_markers(sub3.phar, group_cells_by = "UMAP.cluster")

phar.mark %>% filter(cell_group == 34) %>% arrange(desc(marker_score)) %>% head(20)

ggplot(as.data.frame(colData(sub3.phar)), aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(aes(color = ifelse(Experiment == "unc-86_L1","black", "red")),show.legend = "NA", size = 0.5) + 
  scale_color_manual(values = c("red", "black"))

colData(sub3.phar)$UMAP_1 <- reducedDims(sub3.phar)[["UMAP"]][,1]
colData(sub3.phar)$UMAP_2 <- reducedDims(sub3.phar)[["UMAP"]][,2]

colData(sub3.phar)$Cell.type <- ifelse(
  colData(sub3.phar)$UMAP.cluster %in% c(4,25),
  "I2",
  ifelse(
    colData(sub3.phar)$UMAP.cluster %in% c(26,36,14,19),
    "I1",
    ifelse(
      colData(sub3.phar)$UMAP.cluster %in% c(34,13),
      "I1?",
      ifelse(
        colData(sub3.phar)$UMAP.cluster %in% c(42),
        "I6",
        ifelse(
          colData(sub3.phar)$UMAP.cluster %in% c(15),
          "M1",
          ifelse(colData(sub3.phar)$UMAP.cluster %in% c(9,23,5,31,18),
                 "MC",
                 ifelse(
                   colData(sub3.phar)$UMAP.cluster %in% c(38),
                   "I3",
                   ifelse(
                     colData(sub3.phar)$UMAP.cluster %in% c(41),
                     "M5",
                     ifelse(
                       colData(sub3.phar)$UMAP.cluster %in% c(6,12,10,22,17,8,11,32,27,7,20,39,3,2),
                       "MI",
                       ifelse(
                         colData(sub3.phar)$UMAP.cluster %in% c(37,40),
                         "M4",
                         ifelse(
                           colData(sub3.phar)$UMAP.cluster %in% c(24,30,21,33),
                           "NSM",
                           ifelse(
                             colData(sub3.phar)$UMAP.cluster %in% c(1),
                             "M3",
                             ifelse(
                               colData(sub3.phar)$UMAP.cluster %in% c(35,28,29),
                               "M2",
                               ifelse(
                                 colData(sub3.phar)$UMAP.cluster %in% c(16),
                                 "I4",
                                 as.character(colData(sub3.phar)$Cell.type)
                               )
                             )
                            
                           )
                         )
                       )
                     )
                   )
                 )
          )
        )
      )
   
    )
  )
)

plot_cells(sub3.phar, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(sub3.phar,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(sub3.phar)$Cell.type <- ifelse(
  colData(sub3.phar)$UMAP.cluster %in% c(17,10,8,22,11),
  "MI.UPR",
  ifelse(
    colData(sub3.phar)$UMAP.cluster %in% c(23,5),
    "MC.UPR",
           as.character(colData(sub3.phar)$Cell.type)
    )
  )

# Some possible doublets in the M4, NSM and I1? clusters.
# some possible NSM/M3, cluster 24, part of 30, 21

plot_genes_by_group(sub3.phar, 
                    markers = c("ast-1", "aqp-5", "ceh-2", "tph-1", "bas-1", 
                                "Y53G8AL.4", "degt-1", "R74.10", "flp-5", "flp-13",
                                "lgc-47", "F28F5.4", "ceh-28", "trh-1", "ins-28",
                                "clik-3", "Y60A3A.21", "del-3", "clec-179"), group_cells_by = "UMAP.cluster")

sub4.phar <- sub3.phar[,colData(sub3.phar)$UMAP.cluster %in% c(1, 37, 40, 24, 30, 21, 34, 13)]

sub4.phar <- detect_genes(sub4.phar)
sub4.phar <- sub4.phar[rowData(sub4.phar)$num_cells_expressed > 5,]
sub4.phar
# 8015 features in 1374 cells

sub4.phar <- preprocess_cds(sub4.phar, num_dim = 20)
plot_pc_variance_explained(sub4.phar)

sub4.phar <- align_cds(sub4.phar, aligment_group = "Experiment", alignment_k = 5)
sub4.phar <- reduce_dimension(sub4.phar,
                              reduction_method = "UMAP",
                              preprocess_method = "Aligned",
                              umap.min_dist = 0.08,
                              umap.n_neighbors = 15)

sub4.phar <- cluster_cells(sub4.phar, reduction_method = "UMAP", res = 5e-2)

plot_cells(sub4.phar, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(sub4.phar, 
           color_cells_by = "Experiment",
           cell_size = 0.5,
           label_cell_groups = F)

colData(sub4.phar)$UMAP.cluster <- monocle3::clusters(sub4.phar, reduction_method = "UMAP")

plot_cells(sub4.phar,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(sub4.phar)$UMAP_1 <- reducedDims(sub4.phar)[["UMAP"]][,1]
colData(sub4.phar)$UMAP_2 <- reducedDims(sub4.phar)[["UMAP"]][,2]

plot.expr.UMAP(sub4.phar, "ast-1", size = 0.5)
plot.expr.UMAP(sub4.phar, "aqp-5", size = 0.5)
plot.expr.UMAP(sub4.phar, "unc-25", size = 0.5)
plot.expr.UMAP(sub4.phar, "ceh-28", size = 0.5)
plot.expr.UMAP(sub4.phar, "ceh-2", size = 0.5)

colData(sub4.phar)$Cell.type <- ifelse(
  colData(sub4.phar)$UMAP.cluster %in% c(21,9,1,18, 23, 19),
  "Doublets",
  as.character(colData(sub4.phar)$Cell.type)
)

colData(sub4.phar)$Cell.type <- ifelse(
  colData(sub4.phar)$UMAP.cluster %in% c(4) & colData(sub4.phar)$UMAP_1 > -1.8,
  "Doublets",
  as.character(colData(sub4.phar)$Cell.type)
)

######
colData(sub3.phar)[colnames(sub4.phar),]$Cell.type <- colData(sub4.phar)$Cell.type
colData(L1.sub3)[colnames(sub3.phar),]$Cell.type <- colData(sub3.phar)$Cell.type
colData(L1.neuron)[colnames(L1.sub3),]$Cell.type <- colData(L1.sub3)$Cell.type
colData(neuron)[colnames(L1.neuron),]$Cell.type <- colData(L1.neuron)$Cell.type
colData(combined_cds)[colnames(neuron),]$Cell.type <- colData(neuron)$Cell.type

plot_cells(sub3.phar, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.sub3, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.sub3, 
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = T)

table(colData(L1.sub3)$Cell.type)

sub5 <- L1.sub3[,colData(L1.sub3)$UMAP.cluster %in% c(23, 2,7,5, 51, 56)]
sub5 <- detect_genes(sub5)
sub5 <- sub5[rowData(sub5)$num_cells_expressed > 5,]
sub5
# 10404 features in 6726 cells

sub5 <- preprocess_cds(sub5, num_dim = 20)
plot_pc_variance_explained(sub5)

sub5 <- align_cds(sub5, aligment_group = "Experiment", alignment_k = 5)
sub5 <- reduce_dimension(sub5,
                              reduction_method = "UMAP",
                              preprocess_method = "Aligned",
                              umap.min_dist = 0.08,
                              umap.n_neighbors = 15)

sub5 <- cluster_cells(sub5, reduction_method = "UMAP", res = 3e-3)

plot_cells(sub5, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(sub5)$UMAP.cluster <- monocle3::clusters(sub5, reduction_method = "UMAP")

plot_cells(sub5,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(sub5)$UMAP_1 <- reducedDims(sub5)[["UMAP"]][,1]
colData(sub5)$UMAP_2 <- reducedDims(sub5)[["UMAP"]][,2]

plot.expr.UMAP(sub5, "sbt-1", size = 0.5)
plot.expr.UMAP(sub5, "egl-21", size = 0.5)
plot.expr.UMAP(sub5, "smc-5", size = 0.5)

sub5.mark <- top_markers(sub5, group_cells_by = "UMAP.cluster", marker_sig_test = F)
sub5.mark$gene_short_name <- i2s(sub5.mark$gene_id, gids)

sub5.mark %>% filter(cell_group == 18) %>% arrange(desc(specificity)) %>% head(20)
sub5.mark %>% filter(cell_group == 19) %>% arrange(desc(specificity)) %>% head(20)

colData(sub5)$Cell.type <- ifelse(
  colData(sub5)$UMAP.cluster %in% c(18),
  "Germline",
  as.character(colData(sub5)$Cell.type)
)

colData(sub5)$Cell.type <- ifelse(
  colData(sub5)$UMAP.cluster %in% c(19),
  "Unknown_non_neuronal",
  as.character(colData(sub5)$Cell.type)
)

colData(sub5)$Cell.type <- ifelse(
  colData(sub5)$UMAP.cluster %in% c(17, 12),
  "Pharyngeal_gland_cell",
  as.character(colData(sub5)$Cell.type)
)

colData(sub5)$Cell.type <- ifelse(
  colData(sub5)$UMAP.cluster %in% c(16,9,3,15,4,11,13,5,6,14,10,7,8,1,2),
  "Unannotated",
  as.character(colData(sub5)$Cell.type)
)

plot_cells(sub5, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(sub5,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(L1.sub3)[colnames(sub5),]$Cell.type <- colData(sub5)$Cell.type

table(colData(L1.sub3)$Cell.type, exclude = NULL)

plot_cells(L1.sub3, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.sub3,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(L1.sub3)$Cell.type <- ifelse(
  colData(L1.sub3)$UMAP.cluster %in% c(6,26),
  "AVG",
  as.character(colData(L1.sub3)$Cell.type)
)

colData(L1.sub3)$Cell.type <- ifelse(
  colData(L1.sub3)$UMAP.cluster %in% c(14),
  "ALA",
  as.character(colData(L1.sub3)$Cell.type)
)

colData(L1.neuron)[colnames(L1.sub3),]$Cell.type <- colData(L1.sub3)$Cell.type
colData(neuron)[colnames(L1.neuron),]$Cell.type <- colData(L1.neuron)$Cell.type
colData(combined_cds)[colnames(neuron),]$Cell.type <- colData(neuron)$Cell.type

table(colData(L1.neuron)$Cell.type, exclude = NULL)

colData(L1.neuron)$Cell.type <- dplyr::recode(colData(L1.neuron)$Cell.type,
                                              "Low-quality" = "Low_quality",
                                              "developing_muscle"= "M_lineage_progenitor",
                                              "Pn.p_descendants" = "Pn.p_lineage",
                                              "Pn.p_descendents" = "Pn.p_lineage",
                                              "Unknown" = "Unannotated",
                                              "Unknown_non_neuronal_developing" = "Unknown_non_neuronal",
                                              "Unknown_possible_blast_cell" = "VC",
                                              "Unknown_unc-4+" = "Unannotated",
                                              "Pharyngeal_unknown" = "Unannotated")
plot_cells(L1.neuron, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

# To annotate the Progenitor clusters again, I will use UMAP clustering
# This resets the clustering from earlier in this script
L1.neuron <- cluster_cells(L1.neuron, res = 3e-4)
plot_cells(L1.neuron, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = T)

colData(L1.neuron)$UMAP.cluster <- monocle3::clusters(L1.neuron, reduction_method = "UMAP")

p <- c(62,111,6,3,43,40,53,87,4,46,10,22,39,74,97,61,34,13,27,81,7,93,38,5,24,9,12,69,85,58,30,18,71,32)

ggplot(as.data.frame(colData(L1.neuron)), aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color = ifelse(UMAP.cluster %in% p, "black", "red")), show.legend = "NA", size = 0.5) +
  scale_color_manual(values = c("red", "black"))

L1.p.sub <- L1.neuron[,colData(L1.neuron)$UMAP.cluster %in% p]
L1.p.sub <- detect_genes(L1.p.sub)
L1.p.sub <- L1.p.sub[rowData(L1.p.sub)$num_cells_expressed > 5,]
L1.p.sub
# 14765 features in 46351 cells

L1.p.sub <- preprocess_cds(L1.p.sub, num_dim = 30)
plot_pc_variance_explained(L1.p.sub)

L1.p.sub <- align_cds(L1.p.sub, aligment_group = "Experiment", alignment_k = 5)
L1.p.sub <- reduce_dimension(L1.p.sub,
                         reduction_method = "UMAP",
                         preprocess_method = "Aligned",
                         umap.min_dist = 0.1,
                         umap.n_neighbors = 25)

L1.p.sub <- cluster_cells(L1.p.sub, reduction_method = "UMAP", res = 7e-4)

plot_cells(L1.p.sub, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(L1.p.sub)$UMAP.cluster <- monocle3::clusters(L1.p.sub, reduction_method = "UMAP")

plot_cells(L1.p.sub,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(L1.p.sub)$UMAP_1 <- reducedDims(L1.p.sub)[["UMAP"]][,1]
colData(L1.p.sub)$UMAP_2 <- reducedDims(L1.p.sub)[["UMAP"]][,2]

plot.expr.UMAP(L1.p.sub, "sbt-1", size = 0.5)
plot.expr.UMAP(L1.p.sub, "egl-21", size = 0.5)
plot.expr.UMAP(L1.p.sub, "smc-5", size = 0.5)

L1.p.sub.mark <- top_markers(L1.p.sub, group_cells_by = "UMAP.cluster", marker_sig_test = F)
L1.p.sub.mark$gene_short_name <- i2s(L1.p.sub.mark$gene_id, gids)

p.sub2 <- L1.p.sub[,colData(L1.p.sub)$UMAP.cluster %in% c(49,64,59,13,45,54,62,65,55,5,31,32)]
p.sub2 <- detect_genes(p.sub2)
p.sub2 <- p.sub2[rowData(p.sub2)$num_cells_expressed > 5,]
p.sub2
# 11197 features in 6210 cells

p.sub2 <- preprocess_cds(p.sub2, num_dim = 30)
plot_pc_variance_explained(p.sub2)

p.sub2 <- align_cds(p.sub2, aligment_group = "Experiment", alignment_k = 5)
p.sub2 <- reduce_dimension(p.sub2,
                             reduction_method = "UMAP",
                             preprocess_method = "Aligned",
                             umap.min_dist = 0.1,
                             umap.n_neighbors = 25)

p.sub2 <- cluster_cells(p.sub2, reduction_method = "UMAP", res = 3e-3)

plot_cells(p.sub2, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(p.sub2)$UMAP.cluster <- monocle3::clusters(p.sub2, reduction_method = "UMAP")

plot_cells(p.sub2,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(p.sub2)$UMAP_1 <- reducedDims(p.sub2)[["UMAP"]][,1]
colData(p.sub2)$UMAP_2 <- reducedDims(p.sub2)[["UMAP"]][,2]

plot.expr.UMAP(p.sub2, "F26G1.15", size = 0.5)
plot.expr.UMAP(p.sub2, "egl-21", size = 0.5)
plot.expr.UMAP(p.sub2, "sbt-1", size = 0.5)

colData(p.sub2)$Cell.type <- ifelse(
  colData(p.sub2)$UMAP.cluster %in% c(8),
  "DD",
  ifelse(colData(p.sub2)$UMAP.cluster %in% c(18),
         "DB",
         ifelse(
           colData(p.sub2)$UMAP.cluster %in% c(7) & colData(p.sub2)$Cell.type != "DA1",
           "DA",
           ifelse(
             colData(p.sub2)$UMAP.cluster %in% c(3,4),
             "HSN",
             as.character(colData(p.sub2)$Cell.type)
           )
         )
  )
)

# To see DA1 separate, SAB subclasses, possible DB subsets, possible RMH


p.sub3 <- p.sub2[,colData(p.sub2)$UMAP.cluster %in% c(26,7,11, 18, 8)]

p.sub3 <- detect_genes(p.sub3)
p.sub3 <- p.sub3[rowData(p.sub3)$num_cells_expressed > 5,]
p.sub3
# 6940 features in 944 cells

p.sub3 <- preprocess_cds(p.sub3, num_dim = 20)
plot_pc_variance_explained(p.sub3)

p.sub3 <- align_cds(p.sub3, aligment_group = "Experiment", alignment_k = 5)
p.sub3 <- reduce_dimension(p.sub3,
                           reduction_method = "UMAP",
                           preprocess_method = "Aligned",
                           umap.min_dist = 0.1,
                           umap.n_neighbors = 15)

p.sub3 <- cluster_cells(p.sub3, reduction_method = "UMAP", res = 8e-2)

plot_cells(p.sub3, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(p.sub3)$UMAP.cluster <- monocle3::clusters(p.sub3, reduction_method = "UMAP")

plot_cells(p.sub3,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(p.sub3)$UMAP_1 <- reducedDims(p.sub3)[["UMAP"]][,1]
colData(p.sub3)$UMAP_2 <- reducedDims(p.sub3)[["UMAP"]][,2]

plot.expr.UMAP(p.sub3, "tba-9", size = 0.5)
plot.expr.UMAP(p.sub3, "elt-1", size = 0.5)
plot.expr.UMAP(p.sub3, "flp-10", size = 0.5)
plot.expr.UMAP(p.sub3, "lim-6", size = 0.5)
plot.expr.UMAP(p.sub3, "unc-30", size = 0.5)
plot.expr.UMAP(p.sub3, "flp-10", size = 0.5)
plot.expr.UMAP(p.sub3, "hlh-17", size = 0.5)
plot.expr.UMAP(p.sub3, "hlh-32", size = 0.5)

colData(p.sub3)$Cell.type <- ifelse(
  colData(p.sub3)$UMAP.cluster %in% c(16),
  "DA1",
  ifelse(colData(p.sub3)$UMAP.cluster %in% c(19),
         "RMH",
         ifelse(
           colData(p.sub3)$UMAP.cluster %in% c(8),
           "G_neuroblast",
             as.character(colData(p.sub3)$Cell.type)
           )
         )
  )

colData(p.sub3)$Cell.type <- ifelse(
  colData(p.sub3)$UMAP.cluster %in% c(5,6,10,13),
  "SAB",
  as.character(colData(p.sub3)$Cell.type)
)

colData(p.sub3)$Cell.type <- ifelse(
  colData(p.sub3)$UMAP.cluster %in% c(3) & colData(p.sub3)$UMAP_1 < 2.95 & colData(p.sub3)$UMAP_2 > -3.8,
  "DB02",
  as.character(colData(p.sub3)$Cell.type)
)

p.sub3.mark <- top_markers(p.sub3, group_cells_by = "UMAP.cluster", marker_sig_test = F)
p.sub3.mark$gene_short_name <- i2s(p.sub3.mark$gene_id, gids)

p.sub3.mark %>% filter(cell_group == 1) %>% arrange(desc(specificity)) %>% head(20)

colData(p.sub2)[colnames(p.sub3),]$Cell.type <- colData(p.sub3)$Cell.type

plot_cells(p.sub2, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(p.sub2,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

# Q.p cells

q.sub <- p.sub2[,colData(p.sub2)$UMAP.cluster %in% c(22,12,20,9)]
q.sub <- detect_genes(q.sub)
q.sub <- q.sub[rowData(q.sub)$num_cells_expressed > 5,]
q.sub
# 6296 features in 656 cells

q.sub <- preprocess_cds(q.sub, num_dim = 15)
plot_pc_variance_explained(q.sub)

q.sub <- align_cds(q.sub, aligment_group = "Experiment", alignment_k = 5)
q.sub <- reduce_dimension(q.sub,
                           reduction_method = "UMAP",
                           preprocess_method = "Aligned",
                           umap.min_dist = 0.1,
                           umap.n_neighbors = 15)

q.sub <- cluster_cells(q.sub, reduction_method = "UMAP", res = 8e-2)

plot_cells(q.sub, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(q.sub)$UMAP.cluster <- monocle3::clusters(q.sub, reduction_method = "UMAP")

plot_cells(q.sub,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(q.sub)$UMAP_1 <- reducedDims(q.sub)[["UMAP"]][,1]
colData(q.sub)$UMAP_2 <- reducedDims(q.sub)[["UMAP"]][,2]

plot.expr.UMAP(q.sub, "mec-3", size = 0.5)
plot.expr.UMAP(q.sub, "lin-39", size = 0.5)
plot.expr.UMAP(q.sub, "mab-5", size = 0.5)
plot.expr.UMAP(q.sub, "mec-17", size = 0.5)
plot.expr.UMAP(q.sub, "gcy-35", size = 0.5)
plot.expr.UMAP(q.sub, "flp-8", size = 0.5)
plot.expr.UMAP(q.sub, "cam-1", size = 0.5)
plot.expr.UMAP(q.sub, "unc-86", size = 0.5)
plot.expr.UMAP(q.sub, "F26G1.15", size = 0.5)

colData(q.sub)$Cell.type <- ifelse(
  colData(q.sub)$UMAP.cluster %in% c(8,6,12),
  "PVM",
  ifelse(colData(q.sub)$UMAP.cluster %in% c(10,7,4),
         "AVM",
         ifelse(
           colData(q.sub)$UMAP.cluster %in% c(1,11),
           "QX.pa",
           ifelse(
             colData(q.sub)$UMAP.cluster %in% c(13,5,9,2),
             "SDQ",
             as.character(colData(q.sub)$Cell.type)
           )
           
         )
  )
)

colData(q.sub)$Cell.type <- dplyr::recode(colData(q.sub)$Cell.type,
                                          "50" = "PVM",
                                          "Germline" = "AVM",
                                          "Low_quality" = "PVM")

table(colData(q.sub)$Cell.type)

colData(p.sub2)[colnames(q.sub),]$Cell.type <- colData(q.sub)$Cell.type
colData(L1.p.sub)[colnames(p.sub2),]$Cell.type <- colData(p.sub2)$Cell.type
colData(L1.neuron)[colnames(L1.p.sub),]$Cell.type <- colData(L1.p.sub)$Cell.type
colData(neuron)[colnames(L1.neuron),]$Cell.type <- colData(L1.neuron)$Cell.type
colData(combined_cds)[colnames(neuron),]$Cell.type <- colData(neuron)$Cell.type

plot_cells(p.sub2, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(p.sub2,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

# T lineage
t.sub <- p.sub2[,colData(p.sub2)$UMAP.cluster %in% c(17,10,15,24)]
t.sub <- detect_genes(t.sub)
t.sub <- t.sub[rowData(t.sub)$num_cells_expressed > 5,]
t.sub
# 6386 features in 622 cells

t.sub <- preprocess_cds(t.sub, num_dim = 15)
plot_pc_variance_explained(t.sub)

t.sub <- align_cds(t.sub, aligment_group = "Experiment", alignment_k = 5)
t.sub <- reduce_dimension(t.sub,
                          reduction_method = "UMAP",
                          preprocess_method = "PCA",
                          umap.min_dist = 0.2,
                          umap.n_neighbors = 25)

t.sub <- cluster_cells(t.sub, reduction_method = "UMAP", res = 8e-2)

plot_cells(t.sub, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(t.sub)$UMAP.cluster <- monocle3::clusters(t.sub, reduction_method = "UMAP")

plot_cells(t.sub,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(t.sub)$UMAP_1 <- reducedDims(t.sub)[["UMAP"]][,1]
colData(t.sub)$UMAP_2 <- reducedDims(t.sub)[["UMAP"]][,2]

plot.expr.UMAP(t.sub, "unc-86", size = 0.5)
plot.expr.UMAP(t.sub, "nob-1", size = 0.5)
plot.expr.UMAP(t.sub, "php-3", size = 0.5)
plot.expr.UMAP(t.sub, "lin-32", size = 0.5)
plot.expr.UMAP(t.sub, "egl-1", size = 0.5)
plot.expr.UMAP(t.sub, "lin-44", size = 0.5)
plot.expr.UMAP(t.sub, "lim-7", size = 0.5)
plot.expr.UMAP(t.sub, "ast-1", size = 0.5)
plot.expr.UMAP(t.sub, "sbt-1", size = 0.5)
plot.expr.UMAP(t.sub, "egl-21", size = 0.5)
plot.expr.UMAP(t.sub, "vab-7", size = 0.5)
plot.expr.UMAP(t.sub, "ceh-99", size = 0.5)
plot.expr.UMAP(t.sub, "sem-4", size = 0.5)
plot.expr.UMAP(t.sub, "nlp-20", size = 0.5)
plot.expr.UMAP(t.sub, "gcy-35", size = 0.5)
plot.expr.UMAP(t.sub, "ceh-14", size = 0.5)
plot.expr.UMAP(t.sub, "unc-3", size = 0.5)
plot.expr.UMAP(t.sub, "ceh-63", size = 0.5)
plot.expr.UMAP(t.sub, "tbx-2", size = 0.5)
plot.expr.UMAP(t.sub, "ahr-1", size = 0.5)

colData(t.sub)$Cell.type <- ifelse(
  colData(t.sub)$UMAP.cluster %in% c(1, 10),
  "PVW",
  ifelse(colData(t.sub)$UMAP.cluster %in% c(9,8,5),
         "TX.pppp?",
         ifelse(
           colData(t.sub)$UMAP.cluster %in% c(3,4),
           "TX.pp",
           ifelse(
             colData(t.sub)$UMAP.cluster %in% c(11),
             "PLN",
             ifelse(
               colData(t.sub)$UMAP.cluster %in% c(12),
               "PHC",
             as.character(colData(t.sub)$Cell.type)
           )
           )
         )
  )
)


colData(t.sub)$Cell.type <- dplyr::recode(colData(t.sub)$Cell.type,
                                          "VC" = "PLN",
                                          "Germline" = "TX.pppa?",
                                          "P1.aaa" = "TX.pppa?",
                                          "Pn.ap?" = "TX.ppp?",
                                          "Pn.aa?" = "TX.ppp?",
                                          "Pn.aa" = "TX.pp?",
                                          "VD" = "TX.pp?",
                                          "Germline" = "TX.pp?",
                                          "P0" = "TX.pp?",
                                          "Pn.a?" = "TX.pp?",
                                          "M_lineage_progenitor" = "TX.pp?")

colData(t.sub)$Cell.type <- dplyr::recode(colData(t.sub)$Cell.type,
                                          "TX.pppa?" = "TX.pppa",
                                          "TX.ppp?" = "TX.ppp",
                                          "TX.pp?" = "TX.pp")

colData(p.sub2)$Cell.type <- dplyr::recode(colData(p.sub2)$Cell.type,
                                          "TX.pppa?" = "TX.pppa",
                                          "TX.ppp?" = "TX.ppp",
                                          "TX.pp?" = "TX.pp")


table(colData(t.sub)$Cell.type)

colData(p.sub2)[colnames(t.sub),]$Cell.type <- colData(t.sub)$Cell.type
colData(L1.p.sub)[colnames(p.sub2),]$Cell.type <- colData(p.sub2)$Cell.type
colData(L1.neuron)[colnames(L1.p.sub),]$Cell.type <- colData(L1.p.sub)$Cell.type
colData(neuron)[colnames(L1.neuron),]$Cell.type <- colData(L1.neuron)$Cell.type
colData(combined_cds)[colnames(neuron),]$Cell.type <- colData(neuron)$Cell.type

colData(p.sub2)$Cell.type <- dplyr::recode(colData(p.sub2)$Cell.type,
                                           "Pn.aaa?" = "Body_wall_muscle",
                                           "Anal_sphincter_muscle" = "Body_wall_muscle",
                                           "Germline?" = "Body_wall_muscle",
                                           "QX.pp" = "TX.pppp?")

colData(L1.p.sub)[colnames(p.sub2),]$Cell.type <- colData(p.sub2)$Cell.type

plot_cells(p.sub2, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(p.sub2,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(p.sub2)$Cell.type <- ifelse(
  colData(p.sub2)$UMAP.cluster %in% c(27),
  "Doublets",
  as.character(colData(p.sub2)$Cell.type)
)

colData(p.sub2)$Cell.type <- ifelse(
  colData(p.sub2)$UMAP.cluster %in% c(25),
  "Unknown_non_neuronal",
  as.character(colData(p.sub2)$Cell.type)
)

colData(p.sub2)$Cell.type <- ifelse(
  colData(p.sub2)$UMAP.cluster %in% c(23),
  "Unannotated",
  as.character(colData(p.sub2)$Cell.type)
)

p.sub2.mark <- top_markers(p.sub2, group_cells_by = "UMAP.cluster", marker_sig_test = F)
p.sub2.mark$gene_short_name <- i2s(p.sub2.mark$gene_id, gids)

p.sub2.mark %>% filter(cell_group == 27) %>% arrange(desc(specificity)) %>% head(20)

colData(L1.p.sub)[colnames(p.sub2),]$Cell.type <- colData(p.sub2)$Cell.type
colData(L1.neuron)[colnames(L1.p.sub),]$Cell.type <- colData(L1.p.sub)$Cell.type
colData(neuron)[colnames(L1.neuron),]$Cell.type <- colData(L1.neuron)$Cell.type
colData(combined_cds)[colnames(neuron),]$Cell.type <- colData(neuron)$Cell.type

plot_cells(p.sub2, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(p.sub2,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

plot_cells(L1.p.sub, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.p.sub,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

p.sub4 <- L1.p.sub[,!colData(L1.p.sub)$UMAP.cluster %in% c(65,55,62,45,54,13,59,49,64)]
p.sub4 <- detect_genes(p.sub4)
p.sub4 <- p.sub4[rowData(p.sub4)$num_cells_expressed > 5,]
p.sub4
# 14474 features in 42809 cells

p.sub4 <- preprocess_cds(p.sub4, num_dim = 30)
plot_pc_variance_explained(p.sub4)

p.sub4 <- align_cds(p.sub4, aligment_group = "Experiment", alignment_k = 5)
p.sub4 <- reduce_dimension(p.sub4,
                          reduction_method = "UMAP",
                          preprocess_method = "PCA",
                          umap.min_dist = 0.2,
                          umap.n_neighbors = 25)

p.sub4 <- cluster_cells(p.sub4, reduction_method = "UMAP", res = 3e-4)

table(colData(p.sub4)$Cell.type)

colData(p.sub4)$Cell.type <- dplyr::recode(colData(p.sub4)$Cell.type,
                                           "Pn.a?" = "Pn.a",
                                           "Pn.aa?" = "Pn.aa",
                                           "Pn.ap?" = "Pn.ap",
                                           "Pn.aaa?" = "Pn.aaa",
                                           "TX.pp?" = "TX.pp",
                                           "Pn.p" = "Pn.p_lineage")

plot_cells(p.sub4, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(p.sub4)$UMAP.cluster <- monocle3::clusters(p.sub4, reduction_method = "UMAP")

plot_cells(p.sub4,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(p.sub4)$UMAP_1 <- reducedDims(p.sub4)[["UMAP"]][,1]
colData(p.sub4)$UMAP_2 <- reducedDims(p.sub4)[["UMAP"]][,2]

plot.expr.UMAP(p.sub4, "ceh-13", size = 0.5)
plot.expr.UMAP(p.sub4, "lin-39", size = 0.5)
plot.expr.UMAP(p.sub4, "mab-5", size = 0.5)
plot.expr.UMAP(p.sub4, "egl-5", size = 0.5)
plot.expr.UMAP(p.sub4, "unc-4", size = 0.5)
plot.expr.UMAP(p.sub4, "ceh-12", size = 0.5)
plot.expr.UMAP(p.sub4, "unc-25", size = 0.5)
plot.expr.UMAP(p.sub4, "egl-46", size = 0.5)
plot.expr.UMAP(p.sub4, "unc-3", size = 0.5)
plot.expr.UMAP(p.sub4, "hlh-14", size = 0.5)
plot.expr.UMAP(p.sub4, "hlh-3", size = 0.5)
plot.expr.UMAP(p.sub4, "cnd-1", size = 0.5)
plot.expr.UMAP(p.sub4, "unc-55", size = 0.5)
plot.expr.UMAP(p.sub4, "elt-1", size = 0.5)
plot.expr.UMAP(p.sub4, "egl-44", size = 0.5)
plot.expr.UMAP(p.sub4, "npr-2", size = 0.5)
plot.expr.UMAP(p.sub4, "nlp-38", size = 0.5)

p0 <- p.sub4[,colData(p.sub4)$UMAP.cluster %in% c(27,16,24,9,20,10,17,14,23,21,31,22,26)]
p0 <- detect_genes(p0)
p0 <- p0[rowData(p0)$num_cells_expressed > 5,]
p0
# 11985 features in 14110 cells

p0 <- preprocess_cds(p0, num_dim = 20)
plot_pc_variance_explained(p0)

p0 <- align_cds(p0, aligment_group = "Experiment", alignment_k = 5)
p0 <- reduce_dimension(p0,
                           reduction_method = "UMAP",
                           preprocess_method = "PCA",
                           umap.min_dist = 0.2,
                           umap.n_neighbors = 25)

p0 <- cluster_cells(p0, reduction_method = "UMAP", res = 6e-3)

table(colData(p0)$Cell.type)

plot_cells(p0, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(p0)$UMAP.cluster <- monocle3::clusters(p0, reduction_method = "UMAP")

plot_cells(p0,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(p0)$UMAP_1 <- reducedDims(p0)[["UMAP"]][,1]
colData(p0)$UMAP_2 <- reducedDims(p0)[["UMAP"]][,2]

plot.expr.UMAP(p0, "ceh-13", size = 0.5)
plot.expr.UMAP(p0, "lin-39", size = 0.5)
plot.expr.UMAP(p0, "mab-5", size = 0.5)
plot.expr.UMAP(p0, "egl-5", size = 0.5)
plot.expr.UMAP(p0, "unc-4", size = 0.5)
plot.expr.UMAP(p0, "ceh-12", size = 0.5)
plot.expr.UMAP(p0, "unc-25", size = 0.5)
plot.expr.UMAP(p0, "egl-46", size = 0.5)
plot.expr.UMAP(p0, "unc-3", size = 0.5)
plot.expr.UMAP(p0, "hlh-14", size = 0.5)
plot.expr.UMAP(p0, "hlh-3", size = 0.5)
plot.expr.UMAP(p0, "cnd-1", size = 0.5)
plot.expr.UMAP(p0, "unc-55", size = 0.5)
plot.expr.UMAP(p0, "elt-1", size = 0.5)
plot.expr.UMAP(p0, "egl-44", size = 0.5)
plot.expr.UMAP(p0, "npr-2", size = 0.5)
plot.expr.UMAP(p0, "nlp-38", size = 0.5)

colData(p0)$Cell.type <- ifelse(
  colData(p0)$UMAP.cluster %in% c(53,46,44,22),
  "VD1",
  ifelse(colData(p0)$UMAP.cluster %in% c(40,55,28,48),
         "VA1",
         ifelse(
           colData(p0)$UMAP.cluster %in% c(45,42,62,10,56,11,25,16),
           "P0.p",
           ifelse(
             colData(p0)$UMAP.cluster %in% c(27,18,1,34,29,31,12,7,26),
             "P0",
             ifelse(
               colData(p0)$UMAP.cluster %in% c(20),
               "P0.ap/P1.aap",
               ifelse(
                 colData(p0)$UMAP.cluster %in% c(59, 17, 52, 50, 32, 54),
                 "P0.aa/P1.aaa",
               as.character(colData(p0)$Cell.type)
               )
             )
           )
         )
  )
)

plot_cells(p0, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(p0)$Cell.type <- dplyr::recode(colData(p0)$Cell.type,
                                       "P0.aa" = "P0.aa/P1.aaa")

p0.sub <- p0[,!(colData(p0)$UMAP.cluster %in% c(61,58,8,38,19,36,4,37,14,41,33,51,21,2))]
p0.sub <- detect_genes(p0.sub)
p0.sub <- p0.sub[rowData(p0.sub)$num_cells_expressed > 5,]
p0.sub
# 11429 features in 10848 cells

p0.sub <- preprocess_cds(p0.sub, num_dim = 20)
plot_pc_variance_explained(p0.sub)

p0.sub <- align_cds(p0.sub, aligment_group = "Experiment", alignment_k = 5)
p0.sub <- reduce_dimension(p0.sub,
                       reduction_method = "UMAP",
                       preprocess_method = "PCA",
                       umap.min_dist = 0.2,
                       umap.n_neighbors = 25)

p0.sub <- cluster_cells(p0.sub, reduction_method = "UMAP", res = 1e-2)

table(colData(p0.sub)$Cell.type)

plot_cells(p0.sub, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(p0.sub)$UMAP.cluster <- monocle3::clusters(p0.sub, reduction_method = "UMAP")

plot_cells(p0.sub,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

plot_cells(p0.sub, 
           color_cells_by = "Experiment",
           cell_size = 0.5,
           label_cell_groups = F)

plot_cells(p0.sub, 
           color_cells_by = "hph",
           cell_size = 0.5,
           label_cell_groups = F)

colData(p0.sub)$UMAP_1 <- reducedDims(p0.sub)[["UMAP"]][,1]
colData(p0.sub)$UMAP_2 <- reducedDims(p0.sub)[["UMAP"]][,2]

plot.expr.UMAP(p0.sub, "ceh-13", size = 0.5)
plot.expr.UMAP(p0.sub, "lin-39", size = 0.5)
plot.expr.UMAP(p0.sub, "mab-5", size = 0.5)
plot.expr.UMAP(p0.sub, "egl-5", size = 0.5)
plot.expr.UMAP(p0.sub, "unc-4", size = 0.5)
plot.expr.UMAP(p0.sub, "ceh-12", size = 0.5)
plot.expr.UMAP(p0.sub, "bnc-1", size = 0.5)
plot.expr.UMAP(p0.sub, "hlh-17", size = 0.5)
plot.expr.UMAP(p0.sub, "unc-4", coexpr_gene = "unc-30", size = 0.5)
plot.expr.UMAP(p0.sub, "hlh-32", size = 0.5)
plot.expr.UMAP(p0.sub, "egl-1", size = 0.5)
plot.expr.UMAP(p0.sub, "hlh-14", size = 0.5)
plot.expr.UMAP(p0.sub, "egl-46", size = 0.5)
plot.expr.UMAP(p0.sub, "cnd-1", size = 0.5)

# I think there are some doublets of VA1/VD1 that are unc-3/unc-4/unc-30/unc-25+
plot.expr.UMAP(p0.sub, "unc-3", size = 0.5)
plot.expr.UMAP(p0.sub, "unc-17", size = 0.5)
plot.expr.UMAP(p0.sub, "unc-25", size = 0.5)
plot.expr.UMAP(p0.sub, "cdk-1", size = 0.5)

# what is different about cluster 40 compared to 44/11? I think both are P0.aa/P1.aaa, but I don't know why 
# they are segregating so differently

p0.sub.mark <- top_markers(p0.sub, group_cells_by = "UMAP.cluster", marker_sig_test = F)
p0.sub.mark$gene_short_name <- i2s(p0.sub.mark$gene_id, gids)

p0.sub.mark %>% filter(cell_group == 40) %>% arrange(desc(specificity)) %>% head(20)

DE.two.clust <- function(cds,clust1,clust2, sig = T) {
  sub.cds <- cds[,colData(cds)$UMAP.cluster %in% c(clust1,clust2)]
  sub.cds <- detect_genes(cds)
  sub.cds <- sub.cds[rowData(sub.cds)$num_cells_expressed > 5,]
  de.res <- top_markers(sub.cds, marker_sig_test = sig)
  de.res
}

c40.de <- DE.two(p0.sub, 40, 44)
# there are some differences, but unclear what they mean.

colData(p0.sub)$Cell.type <- ifelse(
  colData(p0.sub)$UMAP.cluster %in% c(59,31,26,15,12,42,56,25),
  "VB02",
  ifelse(colData(p0.sub)$UMAP.cluster %in% c(49,5),
         "VB01",
         ifelse(
           colData(p0.sub)$UMAP.cluster %in% c(14,20,33,61),
           "AVF",
           ifelse(
             colData(p0.sub)$UMAP.cluster %in% c(3),
             "P0.ap/P1.aap",
             ifelse(
               colData(p0.sub)$UMAP.cluster %in% c(46, 39, 45, 7, 62, 6, 41),
               "P0.aa/P1.aaa",
               ifelse(
                 colData(p0.sub)$UMAP.cluster %in% c(1, 17, 16, 34, 47),
                 "P0.a",
                 ifelse(
                   colData(p0.sub)$UMAP.cluster %in% c(53, 2, 8, 24, 27, 38, 19, 9, 22, 28, 57, 4, 13, 23, 51),
                   "P0",
                   ifelse(
                     colData(p0.sub)$UMAP.cluster %in% c(64),
                     "Doublets",
                     as.character(colData(p0.sub)$Cell.type)
                   )
                 )
               )
             )
               )
             )
           )
)

colData(p0.sub)$Cell.type <- ifelse(
  colData(p0.sub)$UMAP_2 < 9 & colData(p0.sub)$UMAP.cluster == 35,
  "Doublets",
  as.character(colData(p0.sub)$Cell.type)
)

plot_cells(p0.sub, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(p0.sub,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(p0)[colnames(p0.sub),]$Cell.type <- colData(p0.sub)$Cell.type   

plot_cells(p0, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(p0.sub,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(p.sub4)[colnames(p0),]$Cell.type <- colData(p0)$Cell.type

plot_cells(p.sub4, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(p.sub4,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

p.sub5 <- p.sub4[,colData(p.sub4)$UMAP.cluster %in% c(7,32,13,18,11,30,1,10,
                                                      12,8,19,14,17,3,29,31,5,25,28,2,
                                                      15,4,6)]
p.sub5 <- detect_genes(p.sub5)
p.sub5 <- p.sub5[rowData(p.sub5)$num_cells_expressed > 5,]
p.sub5
# 14108 features in 33461 cells

p.sub5 <- preprocess_cds(p.sub5, num_dim = 25)
plot_pc_variance_explained(p.sub5)

p.sub5 <- align_cds(p.sub5, aligment_group = "Experiment", alignment_k = 5)
p.sub5 <- reduce_dimension(p.sub5,
                           reduction_method = "UMAP",
                           preprocess_method = "PCA",
                           umap.min_dist = 0.2,
                           umap.n_neighbors = 25)

p.sub5 <- cluster_cells(p.sub5, reduction_method = "UMAP", res = 3e-3)

table(colData(p.sub5)$Cell.type)

plot_cells(p.sub5, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(p.sub5)$UMAP.cluster <- monocle3::clusters(p.sub5, reduction_method = "UMAP")

plot_cells(p.sub5,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

plot_cells(p.sub5, 
           color_cells_by = "Experiment",
           cell_size = 0.5,
           label_cell_groups = F)

plot_cells(p.sub5, 
           color_cells_by = "hph",
           cell_size = 0.5,
           label_cell_groups = F)

colData(p.sub5)$UMAP_1 <- reducedDims(p.sub5)[["UMAP"]][,1]
colData(p.sub5)$UMAP_2 <- reducedDims(p.sub5)[["UMAP"]][,2]

colData(p.sub5)$Cell.type <- ifelse(
 colData(p.sub5)$UMAP.cluster == 0,
  "AVF",
  as.character(colData(p.sub5)$Cell.type)
)

sub6 <- c(21,15,29,7,28,52,22,58,26,
          65,39,43,71,32,8,4,6,30,18,27,35,74,16,51)

ggplot(as.data.frame(colData(p.sub5)), aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color = ifelse(UMAP.cluster %in% sub6, "red", "black")), size = 0.5, 
             show.legend = "NA") +
  scale_color_manual(values = c("black", "red"))

p.sub6 <- p.sub5[,colData(p.sub5)$UMAP.cluster %in% sub6]

p.sub6 <- detect_genes(p.sub6)
p.sub6 <- p.sub6[rowData(p.sub6)$num_cells_expressed > 5,]
p.sub6
# 11750 features in 10754 cells

p.sub6 <- preprocess_cds(p.sub6, num_dim = 20)
plot_pc_variance_explained(p.sub6)

p.sub6 <- align_cds(p.sub6, aligment_group = "Experiment", alignment_k = 5)
p.sub6 <- reduce_dimension(p.sub6,
                           reduction_method = "UMAP",
                           preprocess_method = "PCA",
                           umap.min_dist = 0.1,
                           umap.n_neighbors = 20)

p.sub6 <- cluster_cells(p.sub6, reduction_method = "UMAP", res = 3e-3)

table(colData(p.sub6)$Cell.type)

plot_cells(p.sub6, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(p.sub6)$UMAP.cluster <- monocle3::clusters(p.sub6, reduction_method = "UMAP")

plot_cells(p.sub6,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

plot_cells(p.sub6, 
           color_cells_by = "Experiment",
           cell_size = 0.5,
           label_cell_groups = F)

plot_cells(p.sub6, 
           color_cells_by = "hph",
           cell_size = 0.5,
           label_cell_groups = F)

colData(p.sub6)$UMAP_1 <- reducedDims(p.sub6)[["UMAP"]][,1]
colData(p.sub6)$UMAP_2 <- reducedDims(p.sub6)[["UMAP"]][,2]

plot.expr.UMAP(p.sub6, "ceh-13", size = 0.5)
plot.expr.UMAP(p.sub6, "lin-39", size = 0.5)
plot.expr.UMAP(p.sub6, "mab-5", size = 0.5)
plot.expr.UMAP(p.sub6, "vab-3", size = 0.5)

plot.expr.UMAP(p.sub6, "egl-5", size = 0.5)
plot.expr.UMAP(p.sub6, "unc-4", size = 0.5)
plot.expr.UMAP(p.sub6, "ceh-12", size = 0.5)
plot.expr.UMAP(p.sub6, "bnc-1", size = 0.5)

plot.expr.UMAP(p.sub6, "egl-1", size = 0.5)
plot.expr.UMAP(p.sub6, "hlh-14", size = 0.5)
plot.expr.UMAP(p.sub6, "egl-46", size = 0.5)
plot.expr.UMAP(p.sub6, "cnd-1", size = 0.5)
plot.expr.UMAP(p.sub6, "elt-1", size = 0.5)
plot.expr.UMAP(p.sub6, "unc-55", size = 0.5)
plot.expr.UMAP(p.sub6, "unc-3", size = 0.5)

plot.expr.UMAP(p.sub6, "unc-30", size = 0.5)
plot.expr.UMAP(p.sub6, "unc-25", size = 0.5)
plot.expr.UMAP(p.sub6, "unc-47", size = 0.5)
plot.expr.UMAP(p.sub6, "unc-46", size = 0.5)
plot.expr.UMAP(p.sub6, "oig-1", size = 0.5)
plot.expr.UMAP(p.sub6, "flp-14", size = 0.5)
plot.expr.UMAP(p.sub6, "lntl-1", size = 0.5)
plot.expr.UMAP(p.sub6, "ttr-39", size = 0.5)

p.sub6.mark <- top_markers(p.sub6, group_cells_by = "UMAP.cluster", marker_sig_test = F)
p.sub6.mark$gene_short_name <- i2s(p.sub6.mark$gene_id, gids)

p.sub6.mark %>% filter(cell_group == 31) %>% arrange(desc(specificity)) %>% head(20)

colData(p.sub6)$Cell.type <- ifelse(
  colData(p.sub6)$UMAP.cluster %in% c(21),
  "VD12",
  as.character(colData(p.sub6)$Cell.type)
)

colData(p.sub6)$Cell.type <- ifelse(
  colData(p.sub6)$UMAP.cluster %in% c(24),
  "VD13",
  as.character(colData(p.sub6)$Cell.type)
)

colData(p.sub6)$Cell.type <- ifelse(
  colData(p.sub6)$UMAP.cluster %in% c(4) & colData(p.sub6)$Cell.type != "DD",
  "VD",
  as.character(colData(p.sub6)$Cell.type)
)

colData(p.sub6)$Cell.type <- ifelse(
  colData(p.sub6)$UMAP.cluster %in% c(3,16,26,19,22,2,17,14,8,23,15,1) 
  & !(colData(p.sub6)$Cell.type %in% c("DD", "VD2", "VD12")),
  "VD",
  as.character(colData(p.sub6)$Cell.type)
)

colData(p.sub6)$Cell.type <- ifelse(
  colData(p.sub6)$UMAP.cluster %in% c(25,27,7,18,29,10,20,30,13) ,
  "AS",
  as.character(colData(p.sub6)$Cell.type)
)

colData(p.sub6)$Cell.type <- ifelse(
  colData(p.sub6)$UMAP.cluster %in% c(28,12,9,5,11,6) ,
  "Pn.ap",
  as.character(colData(p.sub6)$Cell.type)
)

colData(p.sub6)$Cell.type <- ifelse(
  colData(p.sub6)$UMAP.cluster %in% c(31) ,
  "Unannotated",
  as.character(colData(p.sub6)$Cell.type)
)

plot.expr.UMAP(p.sub6, "unc-3", size = 0.5)
plot.expr.UMAP(p.sub6, "unc-30", size = 0.5)
plot.expr.UMAP(p.sub6, "lin-39", size = 0.5)
plot.expr.UMAP(p.sub6, "mab-5", size = 0.5)
plot.expr.UMAP(p.sub6, "vab-3", size = 0.5)
plot.expr.UMAP(p.sub6, "egl-5", size = 0.5)

ggplot(as.data.frame(colData(p.sub6)), aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color = ifelse(UMAP.cluster %in% c(13), "red", "black")), size = 0.5, 
             show.legend = "NA") +
  scale_color_manual(values = c("black", "red"))

plot_cells(p.sub6, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(p.sub6,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(p.sub5)[colnames(p.sub6),]$Cell.type <- colData(p.sub6)$Cell.type
colData(p.sub4)[colnames(p.sub5),]$Cell.type <- colData(p.sub5)$Cell.type
colData(L1.p.sub)[colnames(p.sub4),]$Cell.type <- colData(p.sub4)$Cell.type
colData(L1.neuron)[colnames(L1.p.sub),]$Cell.type <- colData(L1.p.sub)$Cell.type
colData(neuron)[colnames(L1.neuron),]$Cell.type <- colData(L1.neuron)$Cell.type
colData(combined_cds)[colnames(neuron),]$Cell.type <- colData(neuron)$Cell.type


sub7 <- c(25,83,82,12,34,78,86,50,64,19,80,88,44,
          24,85,69,33,53,57,13,49,42,46,76,77,
          17,66,68,11,40,36,54,47,56,14,31,87,2,62,
          41,23,81,48,79,20,73,84,61,1,67,72,3,59,
          38,55,37,60,10,63,70, 5, 75, 45)

ggplot(as.data.frame(colData(p.sub5)), aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color = ifelse(UMAP.cluster %in% sub7, "red", "black")), size = 0.5, 
             show.legend = "NA") +
  scale_color_manual(values = c("black", "red"))

p.sub7 <- p.sub5[,colData(p.sub5)$UMAP.cluster %in% sub7]

p.sub7 <- detect_genes(p.sub7)
p.sub7 <- p.sub7[rowData(p.sub7)$num_cells_expressed > 5,]
p.sub7
# 13492 features in 22535 cells

p.sub7 <- preprocess_cds(p.sub7, num_dim = 20)
plot_pc_variance_explained(p.sub7)

p.sub7 <- align_cds(p.sub7, aligment_group = "Experiment", alignment_k = 5)
p.sub7 <- reduce_dimension(p.sub7,
                           reduction_method = "UMAP",
                           preprocess_method = "PCA",
                           umap.min_dist = 0.1,
                           umap.n_neighbors = 20)

p.sub7 <- cluster_cells(p.sub7, reduction_method = "UMAP", res = 3e-3)

table(colData(p.sub7)$Cell.type)

plot_cells(p.sub7, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(p.sub7)$UMAP.cluster <- monocle3::clusters(p.sub7, reduction_method = "UMAP")

plot_cells(p.sub7,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(p.sub7)$UMAP_1 <- reducedDims(p.sub7)[["UMAP"]][,1]
colData(p.sub7)$UMAP_2 <- reducedDims(p.sub7)[["UMAP"]][,2]

plot.expr.UMAP(p.sub7, "unc-3", size = 0.5)
plot.expr.UMAP(p.sub7, "cdk-1", size = 0.5)
plot.expr.UMAP(p.sub7, "hlh-14", size = 0.5)
plot.expr.UMAP(p.sub7, "egl-46", size = 0.5)
plot.expr.UMAP(p.sub7, "hlh-3", size = 0.5)
plot.expr.UMAP(p.sub7, "lin-12", size = 0.5)
plot.expr.UMAP(p.sub7, "lin-31", size = 0.5)
plot.expr.UMAP(p.sub7, "osm-11", size = 0.5)


colData(p.sub7)$Cell.type <- ifelse(
  colData(p.sub7)$UMAP.cluster %in% c(49,45,34,73,69,65,51,5,63,60,59,20,8,75,76,72,55),
  "Pn.p_lineage",
  as.character(colData(p.sub7)$Cell.type)
)

colData(p.sub7)$Cell.type <- ifelse(
  colData(p.sub7)$UMAP.cluster %in% c(32,12,44,2,10,53),
  "VC",
  as.character(colData(p.sub7)$Cell.type)
)

p.sub7.mark <- top_markers(p.sub7, group_cells_by = "UMAP.cluster", marker_sig_test = F)
p.sub7.mark$gene_short_name <- i2s(p.sub7.mark$gene_id, gids)

p.sub7.mark %>% filter(cell_group == 5) %>% arrange(desc(specificity)) %>% head(20)
p.sub7.mark %>% filter(cell_group == 9) %>% arrange(desc(specificity)) %>% head(20)

plot_cells(p.sub7, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(p.sub7,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(p.sub7)$Cell.type <- ifelse(
  colData(p.sub7)$UMAP.cluster %in% c(44),
  "VC",
  as.character(colData(p.sub7)$Cell.type)
)

colData(p.sub7)$Cell.type <- ifelse(
  colData(p.sub7)$UMAP.cluster %in% c(17,18, 47, 52, 24, 25),
  "Low_quality",
  as.character(colData(p.sub7)$Cell.type)
)

colData(p.sub7)$Cell.type <- ifelse(
  colData(p.sub7)$UMAP.cluster %in% c(11,46,21,7),
  "Pn",
  as.character(colData(p.sub7)$Cell.type)
)

p_pa_pp <- c(21,7,62,15,61,4,66,48,14,13,19,57,31,1,56,
             11,46,26,74,42,49,34,45,73,69,65,63,5,
             60,20,59,8,75,72,76,55, 51, 39,16,71,54,68,36,30,29,6,
             53,10,2,44,32,12)

ggplot(as.data.frame(colData(p.sub7)), aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color = ifelse(UMAP.cluster %in% p_pa_pp, "red", "black")), size = 0.5, 
             show.legend = "NA") +
  scale_color_manual(values = c("black", "red"))

p.sub8 <- p.sub7[,colData(p.sub7)$UMAP.cluster %in% p_pa_pp]

p.sub8 <- detect_genes(p.sub8)
p.sub8 <- p.sub8[rowData(p.sub8)$num_cells_expressed > 5,]
p.sub8
# 12964 features in 15215 cells

p.sub8 <- preprocess_cds(p.sub8, num_dim = 20)
plot_pc_variance_explained(p.sub8)

p.sub8 <- align_cds(p.sub8, aligment_group = "Experiment", alignment_k = 5)
p.sub8 <- reduce_dimension(p.sub8,
                           reduction_method = "UMAP",
                           preprocess_method = "PCA",
                           umap.min_dist = 0.1,
                           umap.n_neighbors = 20)

p.sub8 <- cluster_cells(p.sub8, reduction_method = "UMAP", res = 3e-3)

plot_cells(p.sub8, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(p.sub8)$UMAP.cluster <- monocle3::clusters(p.sub8, reduction_method = "UMAP")

plot_cells(p.sub8,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(p.sub8)$UMAP_1 <- reducedDims(p.sub8)[["UMAP"]][,1]
colData(p.sub8)$UMAP_2 <- reducedDims(p.sub8)[["UMAP"]][,2]

ggplot(as.data.frame(colData(p.sub8)), aes(UMAP.cluster, num_genes_expressed)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, height = 1)

plot.expr.UMAP(p.sub8, "sbt-1", size = 0.5)
plot.expr.UMAP(p.sub8, "egl-21", size = 0.5)
plot.expr.UMAP(p.sub8, "egl-3", size = 0.5)
# Possible doublets of post-mitotic and mitotic cells?

plot.expr.UMAP(p.sub8, "osm-11", size = 0.5)
plot.expr.UMAP(p.sub8, "lin-12", size = 0.5)
plot.expr.UMAP(p.sub8, "lin-31", size = 0.5)

plot.expr.UMAP(p.sub8, "hlh-3", size = 0.5)
plot.expr.UMAP(p.sub8, "hlh-4", size = 0.5)
plot.expr.UMAP(p.sub8, "hlh-14", size = 0.5)
plot.expr.UMAP(p.sub8, "cnd-1", size = 0.5)
plot.expr.UMAP(p.sub8, "egl-46", size = 0.5)
plot.expr.UMAP(p.sub8, "ref-1", size = 0.5)
plot.expr.UMAP(p.sub8, "ref-2", size = 0.5)
plot.expr.UMAP(p.sub8, "ceh-6", size = 0.5)
plot.expr.UMAP(p.sub8, "ham-2", size = 0.5)
plot.expr.UMAP(p.sub8, "unc-3", size = 0.5)
plot.expr.UMAP(p.sub8, "unc-62", size = 0.5)
plot.expr.UMAP(p.sub8, "ceh-20", size = 0.5)

plot.expr.UMAP(p.sub8, "ceh-13", size = 0.5)
plot.expr.UMAP(p.sub8, "lin-39", size = 0.5)
plot.expr.UMAP(p.sub8, "mab-5", size = 0.5)

plot.expr.UMAP(p.sub8, "cdk-1", size = 0.5)
plot.expr.UMAP(p.sub8, "cki-1", size = 0.5)
plot.expr.UMAP(p.sub8, "mcm-4", size = 0.5)
plot.expr.UMAP(p.sub8, "cyb-3", size = 0.5)

p.sub8.mark <- top_markers(p.sub8, group_cells_by = "UMAP.cluster", marker_sig_test = F)
p.sub8.mark$gene_short_name <- i2s(p.sub8.mark$gene_id, gids)

p.sub8.mark %>% filter(cell_group == 10) %>% arrange(desc(specificity)) %>% head(20)
# Is this possibly a P1.p cluster? It has strong ceh-13, lin-12, osm-11, lin-31
# also strong lin-33, ref-1, seb-2

# In these UMAPS, it looks like it branches with the P.p lineage cells.
p.sub8.mark %>% filter(cell_group == 50) %>% arrange(desc(specificity)) %>% head(20)
# Maybe a glial cell?

p.sub8.mark %>% filter(cell_group == 51) %>% arrange(desc(specificity)) %>% head(20)
# Nothing particular jumping out, possible neuronal doublets

p.sub8.mark %>% filter(cell_group == 52) %>% arrange(desc(specificity)) %>% head(20)
p.sub8.mark %>% filter(cell_group == 53) %>% arrange(desc(specificity)) %>% head(20)
# Also possible RIA/AIB annd non-neuronal doublets?

p.sub8.mark %>% filter(cell_group == 47) %>% arrange(desc(specificity)) %>% head(20)
# ??

p.sub8.mark %>% filter(cell_group == 13) %>% arrange(desc(specificity)) %>% head(20)

colData(p.sub8)$Cell.type <- ifelse(
  colData(p.sub8)$UMAP.cluster %in% c(12,30,28,26,25,41,47,16,49,7,36,22,14,20,45,13),
  "Pn.p_lineage",
  as.character(colData(p.sub8)$Cell.type)
)

colData(p.sub8)$Cell.type <- ifelse(
  colData(p.sub8)$UMAP.cluster %in% c(51,52,53,50),
  "Doublets",
  as.character(colData(p.sub8)$Cell.type)
)

colData(p.sub8)$Cell.type <- ifelse(
  colData(p.sub8)$UMAP.cluster %in% c(39,10),
  "P1.p?",
  as.character(colData(p.sub8)$Cell.type)
)

plot_cells(p.sub8, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(p.sub8,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(p.sub7)[colnames(p.sub8),]$Cell.type <- colData(p.sub8)$Cell.type
colData(p.sub5)[colnames(p.sub7),]$Cell.type <- colData(p.sub7)$Cell.type
colData(p.sub4)[colnames(p.sub5),]$Cell.type <- colData(p.sub5)$Cell.type

plot_cells(p.sub5, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(p.sub4, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

# Trying without the low quality cells and with only P.a descendents

p.sub5s <- p.sub5[,!colData(p.sub5)$Cell.type %in% c("P1.p?", "Pn.p_lineage", "Doublets",
                                                      "Low_quality", "P0.aa/P1.aaa", 
                                                      "AVF")]
p.sub5s <- detect_genes(p.sub5s)
p.sub5s <- p.sub5s[rowData(p.sub5s)$num_cells_expressed > 5,]
p.sub5s
# 13248 features in 24064 cells

p.sub5s <- preprocess_cds(p.sub5s, num_dim = 15)
plot_pc_variance_explained(p.sub5s)

p.sub5s <- align_cds(p.sub5s, aligment_group = "Experiment", alignment_k = 5)
p.sub5s <- reduce_dimension(p.sub5s,
                           reduction_method = "UMAP",
                           preprocess_method = "PCA",
                           umap.min_dist = 0.1,
                           umap.n_neighbors = 20)

p.sub5s <- cluster_cells(p.sub5s, reduction_method = "UMAP", res = 3e-3)

plot_cells(p.sub5s, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(p.sub5s)$UMAP.cluster <- monocle3::clusters(p.sub5s, reduction_method = "UMAP")

plot_cells(p.sub5s,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(p.sub5s)$UMAP_1 <- reducedDims(p.sub5s)[["UMAP"]][,1]
colData(p.sub5s)$UMAP_2 <- reducedDims(p.sub5s)[["UMAP"]][,2]

plot.expr.UMAP(p.sub5s, "hlh-3", size = 0.5)
plot.expr.UMAP(p.sub5s, "hlh-2", size = 0.5)
plot.expr.UMAP(p.sub5s, "hlh-14", size = 0.5)

plot.expr.UMAP(p.sub5s, "cnd-1", size = 0.5)
plot.expr.UMAP(p.sub5s, "elt-1", size = 0.5)
plot.expr.UMAP(p.sub5s, "unc-55", size = 0.5)
plot.expr.UMAP(p.sub5s, "ref-1", size = 0.5)
plot.expr.UMAP(p.sub5s, "ref-2", size = 0.5)
plot.expr.UMAP(p.sub5s, "ceh-6", size = 0.5)
plot.expr.UMAP(p.sub5s, "ham-2", size = 0.5)
plot.expr.UMAP(p.sub5s, "unc-3", size = 0.5)
plot.expr.UMAP(p.sub5s, "unc-62", size = 0.5)
plot.expr.UMAP(p.sub5s, "ceh-20", size = 0.5)

plot.expr.UMAP(p.sub5s, "egl-1", size = 0.5)

plot.expr.UMAP(p.sub5s, "ceh-13", size = 0.5)
plot.expr.UMAP(p.sub5s, "lin-39", size = 0.5)
plot.expr.UMAP(p.sub5s, "mab-5", size = 0.5)
plot.expr.UMAP(p.sub5s, "egl-5", size = 0.5)

plot.expr.UMAP(p.sub5s, "cdk-1", size = 0.5)
plot.expr.UMAP(p.sub5s, "cki-1", size = 0.5)
plot.expr.UMAP(p.sub5s, "mcm-4", size = 0.5)
plot.expr.UMAP(p.sub5s, "cya-1", size = 0.5)
plot.expr.UMAP(p.sub5s, "cye-1", size = 0.5)
plot.expr.UMAP(p.sub5s, "cyb-1", size = 0.5)
plot.expr.UMAP(p.sub5s, "cyb-2.1", size = 0.5)
plot.expr.UMAP(p.sub5s, "cyb-2.2", size = 0.5)
plot.expr.UMAP(p.sub5s, "cyb-3", size = 0.5)

plot.expr.UMAP(p.sub5s, "cyb-3", size = 0.5)

plot.expr.UMAP(p.sub5s, "sbt-1", size = 0.5)
plot.expr.UMAP(p.sub5s, "egl-21", size = 0.5)

plot.expr.UMAP(p.sub5s, "unc-25", size = 0.5)
plot.expr.UMAP(p.sub5s, "npr-2", size = 0.5)
plot.expr.UMAP(p.sub5s, "ceh-6", size = 0.5)
plot.expr.UMAP(p.sub5s, "unc-4", size = 0.5)
plot.expr.UMAP(p.sub5s, "ceh-12", size = 0.5)

plot_cells(p.sub5s, 
           color_cells_by = "Experiment",
           group_label_size = 4,
           cell_size = 0.5,
           label_cell_groups = F)

ggplot(as.data.frame(colData(p.sub5s)), aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color = ifelse(Experiment == "P-cell_10hr_second", "red", "black")), 
             show.legend = "NA", size = 0.5) + 
  theme_classic() + 
  scale_color_manual(values = c("black", "red"))

ggplot(as.data.frame(colData(p.sub5s)), aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color = ifelse(Experiment == "P-cell_15.5hr", "red", "black")), 
             show.legend = "NA", size = 0.5) + 
  theme_classic() + 
  scale_color_manual(values = c("black", "red"))

ggplot(as.data.frame(colData(p.sub5s)), aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color = ifelse(Experiment == "P-cell_20.5hr_first", "red", "black")), 
             show.legend = "NA", size = 0.5) + 
  theme_classic() + 
  scale_color_manual(values = c("black", "red"))

ggplot(as.data.frame(colData(p.sub5s)), aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color = ifelse(Experiment == "P-cell_20.5hr_second", "red", "black")), 
             show.legend = "NA", size = 0.5) + 
  theme_classic() + 
  scale_color_manual(values = c("black", "red"))

ggplot(as.data.frame(colData(p.sub5s)), aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color = ifelse(Experiment == "L1_mec-3_nlp-13_unc-47_ceh-10", "red", "black")), 
             show.legend = "NA", size = 0.5) + 
  theme_classic() + 
  scale_color_manual(values = c("black", "red"))

colData(p.sub5s)$Cell.type <- ifelse(
  colData(p.sub5s)$UMAP.cluster %in% c(16,61,35,58,11,28,49,44,49,8),
  "Pn.ap",
  as.character(colData(p.sub5s)$Cell.type)
)

colData(p.sub5s)$Cell.type <- ifelse(
  colData(p.sub5s)$UMAP.cluster %in% c(57,14,50,54),
  "Pn.aa",
  as.character(colData(p.sub5s)$Cell.type)
)

colData(p.sub5s)$Cell.type <- ifelse(
  colData(p.sub5s)$UMAP.cluster %in% c(27,1,20,24,30,64,51,25),
  "Pn.a",
  as.character(colData(p.sub5s)$Cell.type)
)

colData(p.sub5s)$Cell.type <- ifelse(
  colData(p.sub5s)$UMAP.cluster %in% c(56),
  "P0.ap/P1.aap",
  as.character(colData(p.sub5s)$Cell.type)
)

colData(p.sub5s)$Cell.type <- ifelse(
  colData(p.sub5s)$UMAP.cluster %in% c(47) & colData(p.sub5s)$UMAP_1 < 3.5,
  "P0.ap/P1.aap",
  as.character(colData(p.sub5s)$Cell.type)
)

colData(p.sub5s)$Cell.type <- ifelse(
  colData(p.sub5s)$UMAP.cluster %in% c(47) & colData(p.sub5s)$UMAP_1 > 3.5,
  "Pn.aa",
  as.character(colData(p.sub5s)$Cell.type)
)

plot_cells(p.sub5s, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(p.sub5s,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

plot_cells(p.sub5s,
           color_cells_by = "Experiment",
           label_cell_groups = F,
           cell_size = 0.5)

colData(p.sub5)[colnames(p.sub5s),]$Cell.type <- colData(p.sub5s)$Cell.type
colData(p.sub4)[colnames(p.sub5),]$Cell.type <- colData(p.sub5)$Cell.type

# Just VA and VB from P-lineage clusters
p.sub9 <- p.sub5s[,colData(p.sub5s)$UMAP.cluster %in% c(65,66,32,59,2,
                                                        42,53,4,69,19,17,43,39)]

p.sub9 <- detect_genes(p.sub9)
p.sub9 <- p.sub9[rowData(p.sub9)$num_cells_expressed > 5,]
p.sub9
# 9423 features in 4218 cells

p.sub9 <- preprocess_cds(p.sub9, num_dim = 15)
plot_pc_variance_explained(p.sub9)

p.sub9 <- align_cds(p.sub9, aligment_group = "Experiment", alignment_k = 5)
p.sub9 <- reduce_dimension(p.sub9,
                            reduction_method = "UMAP",
                            preprocess_method = "PCA",
                            umap.min_dist = 0.1,
                            umap.n_neighbors = 15)

p.sub9 <- cluster_cells(p.sub9, reduction_method = "UMAP", res = 3e-2)

plot_cells(p.sub9, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(p.sub9)$UMAP.cluster <- monocle3::clusters(p.sub9, reduction_method = "UMAP")

plot_cells(p.sub9,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

plot_cells(p.sub9,
           color_cells_by = "Experiment",
           label_cell_groups = F,
           cell_size = 0.5)

colData(p.sub9)$UMAP_1 <- reducedDims(p.sub9)[["UMAP"]][,1]
colData(p.sub9)$UMAP_2 <- reducedDims(p.sub9)[["UMAP"]][,2]

plot.expr.UMAP(p.sub9, "hlh-3", size = 0.5)
plot.expr.UMAP(p.sub9, "hlh-2", size = 0.5)
plot.expr.UMAP(p.sub9, "hlh-14", size = 0.5)

plot.expr.UMAP(p.sub9, "cnd-1", size = 0.5)
plot.expr.UMAP(p.sub9, "elt-1", size = 0.5)

plot.expr.UMAP(p.sub9, "unc-4", size = 0.5)
plot.expr.UMAP(p.sub9, "bnc-1", size = 0.5)
plot.expr.UMAP(p.sub9, "ceh-12", size = 0.5)
plot.expr.UMAP(p.sub9, "unc-52", size = 0.5)
plot.expr.UMAP(p.sub9, "acr-14", size = 0.5)
plot.expr.UMAP(p.sub9, "eva-1", size = 0.5)

plot.expr.UMAP(p.sub9, "egl-1", size = 0.5)

plot.expr.UMAP(p.sub9, "ceh-13", size = 0.5)
plot.expr.UMAP(p.sub9, "lin-39", size = 0.5)
plot.expr.UMAP(p.sub9, "mab-5", size = 0.5)
plot.expr.UMAP(p.sub9, "egl-5", size = 0.5)
plot.expr.UMAP(p.sub9, "cdk-1", size = 0.5)
plot.expr.UMAP(p.sub9, "cyb-3", size = 0.5)
plot.expr.UMAP(p.sub9, "mcm-4", size = 0.5)
plot.expr.UMAP(p.sub9, "egl-21", size = 0.5)
plot.expr.UMAP(p.sub9, "egl-3", size = 0.5)
plot.expr.UMAP(p.sub9, "unc-4", coexpr_gene = "ceh-12", size = 0.5)

p.sub9.mark <- top_markers(p.sub9, group_cells_by = "UMAP.cluster", marker_sig_test = F)
p.sub9.mark$gene_short_name <- i2s(p.sub9.mark$gene_id, gids)

p.sub9.mark %>% filter(cell_group == 65) %>% arrange(desc(specificity)) %>% head(20)
p.sub9.mark %>% filter(cell_group == 9) %>% arrange(desc(specificity)) %>% head(20)
p.sub9.mark %>% filter(cell_group == 21) %>% arrange(desc(specificity)) %>% head(20)

# What are the hierarchies of subsets I've used for annotations?
# combined_cds > neuron > L1.neuron > L1.sub.out > L1.sub.out2

# L1.neuron > L1.sub3 > L1.sub3.2 > sub4.DA
# L1.sub3 > sub3.phar > sub4.phar
# L1.sub3 > sub5

# L1.neuron > L1.p.sub > p.sub2 > p.sub3 
# p.sub2 > q.sub
# p.sub2 > t.sub

# L1.p.sub > p.sub4 > p0 > p0.sub
# p.sub4 > p.sub5 > p.sub6
# p.sub4 > p.sub5 > p.sub5s > p.sub9
# p.sub4 > p.sub5 > p.sub7 > p.sub8
# L1.p.sub > p.sub10

colData(p.sub9)$Cell.type <- dplyr::recode(colData(p.sub9)$Cell.type,
                                           "VA_VB" = "VA_VB_Doublets")

colData(p.sub9)$Cell.type <- ifelse(
  colData(p.sub9)$UMAP.cluster %in% c(18,32,27,11,29,51),
  "VA_VB_Doublets",
  as.character(colData(p.sub9)$Cell.type)
)

colData(p.sub9)$Cell.type <- ifelse(
  colData(p.sub9)$UMAP.cluster %in% c(24,20,21,48,28,41,19,26),
  "VB",
  as.character(colData(p.sub9)$Cell.type)
)

table(colData(p.sub9)$Cell.type)
colData(p.sub9)$Cell.type <- ifelse(
  colData(p.sub9)$UMAP.cluster %in% c(36,4,10,3,15,14,50,1,52,16,34,37,38,33,45,22,8,46,31,23),
  "VA",
  as.character(colData(p.sub9)$Cell.type)
)

plot_cells(p.sub9, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(p.sub9,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(p.sub9)$Cell.type <- ifelse(
  colData(p.sub9)$UMAP.cluster %in% c(53,7,2,9,17,40,6,12,42,30,35,39,44,43,5,25,13),
  "Pn.aaa",
  as.character(colData(p.sub9)$Cell.type)
)

colData(p.sub5s)[colnames(p.sub9),]$Cell.type <- colData(p.sub9)$Cell.type

plot_cells(p.sub5s, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(p.sub5s,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(p.sub5s)$Cell.type <- ifelse(
  colData(p.sub5s)$UMAP.cluster %in% c(33,9),
  "Pn",
  as.character(colData(p.sub5s)$Cell.type)
)

colData(p.sub5s)$Cell.type <- dplyr::recode(colData(p.sub5s)$Cell.type,
                                           "VA_VB" = "VA_VB_Doublets")

colData(p.sub4)$Cell.type <- dplyr::recode(colData(p.sub4)$Cell.type,
                                            "AMso" = "VB02",
                                           "Unknown_developing" = "P0.aa/P1.aaa",
                                           "VA_VB" = "VA_VB_Doublets")

colData(L1.p.sub)[colnames(p.sub2),]$Cell.type <- colData(p.sub2)$Cell.type

colData(L1.p.sub)$Cell.type <- dplyr::recode(colData(L1.p.sub)$Cell.type,
                                           "Pn.a?" = "Pn.a",
                                           "1" = "Body_wall_muscle", # cluster 55,
                                           "URA" = "HSN")

colData(p.sub5)[colnames(p.sub5s),]$Cell.type <- colData(p.sub5s)$Cell.type
colData(p.sub4)[colnames(p.sub5),]$Cell.type <- colData(p.sub5)$Cell.type

# To resolve egl-1 positive cells from Q/T lineage and from P lineage
p.sub10 <- L1.p.sub[,colData(L1.p.sub)$UMAP.cluster %in% c(47,25,5,65,62)]
p.sub10 <- detect_genes(p.sub10)
p.sub10 <- p.sub10[rowData(p.sub10)$num_cells_expressed > 5,]
p.sub10
# 8842 features in 3016 cells

p.sub10 <- preprocess_cds(p.sub10, num_dim = 15)
plot_pc_variance_explained(p.sub10)

p.sub10 <- align_cds(p.sub10, aligment_group = "Experiment", alignment_k = 5)
p.sub10 <- reduce_dimension(p.sub10,
                           reduction_method = "UMAP",
                           preprocess_method = "PCA",
                           umap.min_dist = 0.1,
                           umap.n_neighbors = 15)

p.sub10 <- cluster_cells(p.sub10, reduction_method = "UMAP", res = 3e-2)

plot_cells(p.sub10, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(p.sub10)$UMAP.cluster <- monocle3::clusters(p.sub10, reduction_method = "UMAP")

plot_cells(p.sub10,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

plot_cells(p.sub10,
           color_cells_by = "Experiment",
           label_cell_groups = F,
           cell_size = 0.5)

colData(p.sub10)$UMAP_1 <- reducedDims(p.sub10)[["UMAP"]][,1]
colData(p.sub10)$UMAP_2 <- reducedDims(p.sub10)[["UMAP"]][,2]

plot.expr.UMAP(p.sub10, "egl-1", size = 0.5)
plot.expr.UMAP(p.sub10, "unc-86", size = 0.5)
plot.expr.UMAP(p.sub10, "unc-62", size = 0.5)
plot.expr.UMAP(p.sub10, "nob-1", size = 0.5)
plot.expr.UMAP(p.sub10, "lin-39", size = 0.5)
plot.expr.UMAP(p.sub10, "mab-5", size = 0.5)

colData(p.sub10)$Cell.type <- ifelse(
  colData(p.sub10)$UMAP.cluster %in% c(8,40),
  "TX.pppp?",
  as.character(colData(p.sub10)$Cell.type)
)

colData(L1.p.sub)[colnames(p.sub10),]$Cell.type <- colData(p.sub10)$Cell.type

table(colData(L1.p.sub)$Cell.type)
colData(L1.neuron)[colnames(L1.p.sub),]$Cell.type <- colData(L1.p.sub)$Cell.type
colData(neuron)[colnames(L1.neuron),]$Cell.type <- colData(L1.neuron)$Cell.type
colData(combined_cds)[colnames(neuron),]$Cell.type <- colData(neuron)$Cell.type

table(colData(L1.neuron)$Cell.type)
# some to clean up here as well - 11, Coelomocyte, Epidermis, Germline?, Glia, Intestine, 
# Low_quality_hypodermal, Muscle, Pharynx, RIC_PVT, Unknown_developing

colData(L1.neuron)$Cell.type <- dplyr::recode(colData(L1.neuron)$Cell.type, 
                                              "Glia" = "BAG",
                                              "Germline?" = "Germline",
                                              "Coelomocyte" = "Unannotated",
                                              "RIC_PVT" = "RIC",
                                              "Pharynx" = "AIY",
                                              "11" = "Unannotated",
                                              "Muscle" = "hmc",
                                              "Epidermis" = "Unannotated",
                                              "Intestine" = "Unannotated",
                                              "Low_quality_hypodermal" = "Low_quality",
                                              "Unknown_developing" = "Unannotated")

table(colData(L1.neuron)$Cell.type)

colData(neuron)[colnames(L1.neuron),]$Cell.type <- colData(L1.neuron)$Cell.type
colData(combined_cds)[colnames(neuron),]$Cell.type <- colData(neuron)$Cell.type

colData(L1.neuron)$Tissue <- ifelse(
  colData(L1.neuron)$Cell.type %in% c("Body_wall_muscle", "Germline", "hmc", "Pharyngeal_gland_cell",
                                      "Unknown_non_neuronal", "M_lineage_progenitor"),
  "Non-neuron",
  ifelse(
    colData(L1.neuron)$Cell.type %in% c("AIN_doublets?", "Doublets", "VA_VB_Doublets"),
    "Doublets",
    ifelse(
      colData(L1.neuron)$Cell.type %in% c("G_neuroblast", "P0", "P0.a", "P0.aa/P1.aaa", "P0.ap/P1.aap",
                                          "Pn", "Pn.a", "Pn.aa", "Pn.aaa", "Pn.ap", "QX.pa",
                                          "TX.pp", "TX.ppp", "TX.pppa", "TX.pppp?", "P0.p", "Pn.p_lineage",
                                          "P1.p?"),
      "Progenitor",
      ifelse(
        colData(L1.neuron)$Cell.type %in% c("Unannotated", "Low_quality"),
        "Unknown",
        "Neuron"
      )
    )
  )
)
table(colData(L1.neuron)$Tissue, exclude = NULL)

table(colData(neuron)$Tissue, exclude = NULL)

colData(neuron)[colnames(L1.neuron),]$Tissue <- colData(L1.neuron)$Tissue

colData(neuron)$Cell.type <- dplyr::recode(colData(neuron)$Cell.type, 
                                              "Unknown" = "Unannotated",
                                              "Unknown_non_neuronal_developing" = "Unannotated",
                                              "Unknown_developing" = "Unannotated",
                                              "Unknown_neuron" = "Unannotated",
                                              "Pn.p_descendents" = "Pn.p_lineage",
                                              "Pharyngeal_neurons" = "Unannotated", 
                                              "Pharyngeal_unknown" = "Pharynx",
                                              "Low_quality_Pn.p_descendants" = "Pn.p_lineage",
                                              "G_neuroblast?" = "G_neuroblast",
                                              "developing_muscle" = "M_lineage_progenitor",
                                              "Low-quality" = "Low_quality")

colData(combined_cds)[colnames(neuron),]$Cell.type <- colData(neuron)$Cell.type

colData(neuron)$Tissue <- ifelse(
  colData(neuron)$Cell.type %in% c("Body_wall_muscle", "Germline", "hmc", "Pharyngeal_gland_cell",
                                      "Unknown_non_neuronal", "AMsh", "AMsh_PHsh", "AMso", "AMso_PHso",
                                      "Anal_muscle", "Anal_sphincter_muscle", "Arcade_cell",
                                      "CEPsh_Glia_1_2", "Coelomocyte", "Epidermis", "Excretory_gland_cell",
                                      "Germline", "Glia", "Glia_4", "Glia_5", "Glia_hypodermis", "Glia_sheath",
                                      "Glia_socket", "Intestine", "Marginal_cell", "Muscle", "Pharyngeal_cells",
                                      "Pharyngeal_gland_cell", "Pharyngeal_intestinal_valve", "Pharyngeal_muscle",
                                      "Pharynx", "PHsh", "PHso", "Rectal_gland", "Rectal_unknown_1",
                                      "Rectal_unknown_2", "Rectal_unknown_3", "Seam_cell", "XXX",
                                      "M_lineage_progenitor"),
  "Non-neuron",
  ifelse(
    colData(neuron)$Cell.type %in% c("AIN_doublets?", "Doublets", "VA_VB_Doublets"),
    "Doublets",
    ifelse(
      colData(neuron)$Cell.type %in% c("G_neuroblast", "P0", "P0.a", "P0.aa/P1.aaa", "P0.ap/P1.aap",
                                          "Pn", "Pn.a", "Pn.aa", "Pn.aaa", "Pn.ap", "QX.pa",
                                          "TX.pp", "TX.ppp", "TX.pppa", "TX.pppp?", "P0.p", "Pn.p_lineage",
                                          "P1.p?"),
      "Progenitor",
      ifelse(
        colData(neuron)$Cell.type %in% c("Unannotated", "Low_quality", "Low_quality_glia_hypodermis",
                                              "Low_quality_hypodermal"),
        "Unknown",
        "Neuron"
      )
    )
  )
)
table(colData(neuron)$Tissue, exclude = NULL)
colData(combined_cds)[colnames(neuron),]$Tissue <- colData(neuron)$Tissue

saveRDS(L1.neuron, "~/Dropbox (VU Basic Sciences)/L1 Cengen manuscript/Claire data and presentations/no_intron_analysis/data/042424_all_neurons_combined_cds.rds")
saveRDS(neuron, "~/Dropbox (VU Basic Sciences)/L1 Cengen manuscript/Claire data and presentations/no_intron_analysis/data/042424_all_cells_combined_filtered_cds.rds")
saveRDS(combined_cds, "~/Dropbox (VU Basic Sciences)/L1 Cengen manuscript/Claire data and presentations/no_intron_analysis/data/042424_all_cells_combined_no_filter_cds.rds")

# I can now save some of the subsets

# GABA.sub to search for possible DVB cluster

gaba <- L1.neuron[,colData(L1.neuron)$Cell.type %in% c("AVL", "DD", "VD", "VD1", "VD2", "VD12", "VD12", "VD13", 
                                                       "RIB", "RME_LR", "RME_DV")]

gaba <- detect_genes(gaba)
gaba <- gaba[rowData(gaba)$num_cells_expressed > 5,]
gaba
# 10706 features in 6742 cells

gaba <- preprocess_cds(gaba, num_dim = 20)
plot_pc_variance_explained(gaba)

gaba <- align_cds(gaba, aligment_group = "Experiment", alignment_k = 5)
gaba <- reduce_dimension(gaba,
                            reduction_method = "UMAP",
                            preprocess_method = "PCA",
                            umap.min_dist = 0.1,
                            umap.n_neighbors = 15)

gaba <- cluster_cells(gaba, reduction_method = "UMAP", res = 3e-2)

plot_cells(gaba, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

colData(gaba)$UMAP.cluster <- monocle3::clusters(gaba, reduction_method = "UMAP")

plot_cells(gaba,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

plot_cells(gaba,
           color_cells_by = "Experiment",
           label_cell_groups = F,
           cell_size = 0.5)

colData(gaba)$UMAP_1 <- reducedDims(gaba)[["UMAP"]][,1]
colData(gaba)$UMAP_2 <- reducedDims(gaba)[["UMAP"]][,2]

plot.expr.UMAP(gaba, "flp-10", size = 0.5)
plot.expr.UMAP(gaba, "egl-5", size = 0.5)
plot.expr.UMAP(gaba, "unc-25", size = 0.5)
plot.expr.UMAP(gaba, "flp-28", size = 0.5)
plot.expr.UMAP(gaba, "nlp-40", size = 0.5)
plot.expr.UMAP(gaba, "flp-22", size = 0.5)
plot.expr.UMAP(gaba, "unc-47", size = 0.5)
plot.expr.UMAP(gaba, "D1007.19", size = 0.5)
plot.expr.UMAP(gaba, "T05A8.3", size = 0.5)
plot.expr.UMAP(gaba, "T02C12.5", size = 0.5)
plot.expr.UMAP(gaba, "nlp-6", size = 0.5)
plot.expr.UMAP(gaba, "flp-6", size = 0.5)
plot.expr.UMAP(gaba, "nlp-21", size = 0.5)

# All look like AVL

# what TFs are different between them?


# New filtered and neuronal datasets

L1.clean.filtered <- combined_cds[,!colData(combined_cds)$Cell.type %in% c("ADF.UPR", "AIB.UPR", "ASK.UPR", "AVJ.UPR", 
                                                               "Doublets", "Doublets_I1", "I1.UPR", "I2.UPR", "I5.UPR", 
                                                               "MC.UPR", "MI.UPR", "NSM.UPR", "PHsh_PHso_doublets?",
                                                               "RIA.UPR", "RIG.UPR", "RIS.UPR", "RIV.UPR", "RMD_LR.UPR", 
                                                               "SAA.doublets", "SIB.UPR", "SMD.UPR",
                                                               "VA_VB_Doublets", "AIN_doublets?")]

L1.clean.filtered <- detect_genes(L1.clean.filtered)
L1.clean.filtered <- L1.clean.filtered[rowData(L1.clean.filtered)$num_cells_expressed > 5,]
L1.clean.filtered
# 18570 features in 161562 cells

L1.clean.filtered <- preprocess_cds(L1.clean.filtered, num_dim = 75)
plot_pc_variance_explained(L1.clean.filtered)

L1.clean.filtered <- align_cds(L1.clean.filtered, aligment_group = "Experiment", alignment_k = 5)
L1.clean.filtered <- reduce_dimension(L1.clean.filtered,
                                reduction_method = "UMAP",
                                preprocess_method = "Aligned",
                                umap.min_dist = 0.3,
                                umap.n_neighbors = 75)

L1.clean.filtered <- cluster_cells(L1.clean.filtered, reduction_method = "UMAP", res = 3e-4)

plot_cells(L1.clean.filtered, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.clean.filtered, 
           color_cells_by = "Experiment",
           label_cell_groups = F,
           cell_size = 0.5)

plot_cells(L1.clean.filtered, 
           color_cells_by = "hph",
           label_cell_groups = F,
           cell_size = 0.5)

colData(L1.clean.filtered)$UMAP.cluster <- monocle3::clusters(L1.clean.filtered, reduction_method = "UMAP")

plot_cells(L1.clean.filtered,
           color_cells_by = "PCA.partition",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

plot_cells(L1.clean.filtered,
           color_cells_by = "PCA.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(L1.clean.filtered)$UMAP_1 <- reducedDims(L1.clean.filtered)[["UMAP"]][,1]
colData(L1.clean.filtered)$UMAP_2 <- reducedDims(L1.clean.filtered)[["UMAP"]][,2]

plot.expr.UMAP(L1.clean.filtered, "sbt-1", size = 0.5)
plot.expr.UMAP(L1.clean.filtered, "cdk-1", size = 0.5)
plot.expr.UMAP(L1.clean.filtered, "sars-1", size = 0.5)

table(colData(L1.clean.filtered)$Cell.type, exclude= NULL)
table(colData(L1.clean.filtered)$Tissue, exclude= NULL)

colData(L1.clean.filtered)$Tissue <- ifelse(
  colData(L1.clean.filtered)$Cell.type %in% c("Unannotated", "Low_quality", "Low_quality_glia_hypodermis",
                                        "Low_quality_hypodermal"),
  "Unknown",
  as.character(colData(L1.clean.filtered)$Tissue))

saveRDS(L1.clean.filtered, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/042324_L1_all_live_cells_MeOH_new_saa_P_lineage_no_doublets_or_UPR_cds.rds")

# Saving subsets
# L1.neuron > L1.sub.out > L1.sub.out2

# L1.neuron > L1.sub3 > L1.sub3.2 > sub4.DA
# L1.sub3 > sub3.phar > sub4.phar
# L1.sub3 > sub5

# L1.neuron > L1.p.sub > p.sub2 > p.sub3 
# p.sub2 > q.sub
# p.sub2 > t.sub

# L1.p.sub > p.sub4 > p0 > p0.sub
# p.sub4 > p.sub5 > p.sub6
# p.sub4 > p.sub5 > p.sub5s > p.sub9
# p.sub4 > p.sub5 > p.sub7 > p.sub8
# L1.p.sub > p.sub10

saveRDS(L1.sub.out, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_neuron_subset_1_cds")
saveRDS(L1.sub.out2, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_neuron_subset_2_cds")
saveRDS(L1.sub3, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_neuron_subset_3_cds")
saveRDS(L1.sub3.2, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_neuron_subset_3.2_cds")
saveRDS(sub4.DA, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_neuron_subset_3_dopamine_cds")
saveRDS(sub3.phar, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_neuron_subset_pharyngeal_cds")
saveRDS(sub4.phar, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_neuron_subset_pharyngeal_2_cds")
saveRDS(p.sub2, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_progenitor_newborn_subset_2_cds")
saveRDS(p.sub3, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_progenitor_newborn_subset_3_cds")
saveRDS(q.sub, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_progenitor_newborn_subset_q_lineage_cds")
saveRDS(t.sub, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_progenitor_newborn_subset_t_lineage_cds")
saveRDS(L1.p.sub, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_progenitor_newborn_subset_full_cds")
saveRDS(p.sub4, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_progenitor_newborn_subset_4_cds")
saveRDS(p0, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_progenitor_newborn_subset_P0_cds")
saveRDS(p0.sub, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_progenitor_newborn_subset_P0_sub2_cds")
saveRDS(p.sub5, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_progenitor_newborn_subset_5_cds")
saveRDS(p.sub6, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_progenitor_newborn_subset_6_cds")
saveRDS(p.sub5s, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_progenitor_newborn_subset_5_sub_cds")
saveRDS(p.sub7, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_progenitor_newborn_subset_7_cds")
saveRDS(p.sub8, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_progenitor_newborn_subset_8_cds")
saveRDS(p.sub9, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_progenitor_newborn_subset_9_cds")
saveRDS(p.sub10, "./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/050124_L1_progenitor_newborn_subset_10_cds")

L1.clean.non <- L1.clean.filtered[,!(colData(L1.clean.filtered)$Tissue %in% c("Neuron", "Progenitor"))]
# non-neuronal

L1.clean.non <- detect_genes(L1.clean.non)
L1.clean.non <- L1.clean.non[rowData(L1.clean.non)$num_cells_expressed > 5,]
L1.clean.non
# 16120 features in 49446 cells

L1.clean.non <- preprocess_cds(L1.clean.non, num_dim = 55)
plot_pc_variance_explained(L1.clean.non)

table(colData(L1.clean.non)$Cell.type)

L1.clean.non <- align_cds(L1.clean.non, aligment_group = "Experiment", alignment_k = 10)
L1.clean.non <- reduce_dimension(L1.clean.non,
                                      reduction_method = "UMAP",
                                      preprocess_method = "Aligned",
                                      umap.min_dist = 0.3,
                                      umap.n_neighbors = 75)

L1.clean.non <- cluster_cells(L1.clean.non, reduction_method = "UMAP", res = 3e-3)

plot_cells(L1.clean.non, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(L1.clean.non, 
           color_cells_by = "Experiment",
           label_cell_groups = F,
           cell_size = 0.5)

plot_cells(L1.clean.non, 
           color_cells_by = "hph",
           label_cell_groups = F,
           cell_size = 0.5)

colData(L1.clean.non)$UMAP.cluster <- monocle3::clusters(L1.clean.non, reduction_method = "UMAP")

plot_cells(L1.clean.non,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(L1.clean.non)$UMAP_1 <- reducedDims(L1.clean.non)[["UMAP"]][,1]
colData(L1.clean.non)$UMAP_2 <- reducedDims(L1.clean.non)[["UMAP"]][,2]

plot.expr.UMAP(L1.clean.non, "sbt-1", size = 0.5)
plot.expr.UMAP(L1.clean.non, "cdk-1", size = 0.5)
plot.expr.UMAP(L1.clean.non, "cup-4", size = 0.5)
plot.expr.UMAP(L1.clean.non, "pat-10", size = 0.5)
plot.expr.UMAP(L1.clean.non, "nspc-19", size = 0.5)
plot.expr.UMAP(L1.clean.non, "ceh-34", size = 0.5)

colData(L1.clean.non)$Cell.type <- ifelse(
  colData(L1.clean.non)$UMAP.cluster %in% c(115,123),
  "PHsh",
  ifelse(
    colData(L1.clean.non)$UMAP.cluster %in% c(116),
    "AMsh",
    as.character(colData(L1.clean.non)$Cell.type)
  )
)

plot.expr.UMAP(L1.clean.non, "tnt-3", size = 0.5)

# glial cells

glia.non <- L1.clean.non[,colData(L1.clean.non)$UMAP.cluster %in% c(85,87,75,119,55,
                                                                    76,103,127,123,115,116)]

glia.non <- detect_genes(glia.non)
glia.non <- glia.non[rowData(glia.non)$num_cells_expressed > 5,]
glia.non
# 10805 features in 2854 cells

glia.non <- preprocess_cds(glia.non, num_dim = 25)
plot_pc_variance_explained(glia.non)

table(colData(glia.non)$Cell.type)

glia.non <- align_cds(glia.non, aligment_group = "Experiment", alignment_k = 10)
glia.non <- reduce_dimension(glia.non,
                                 reduction_method = "UMAP",
                                 preprocess_method = "Aligned",
                                 umap.min_dist = 0.2,
                                 umap.n_neighbors = 25)

glia.non <- cluster_cells(glia.non, reduction_method = "UMAP", res = 3e-2)

plot_cells(glia.non, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           cell_size = 0.5,
           label_groups_by_cluster = F)

plot_cells(glia.non, 
           color_cells_by = "Experiment",
           label_cell_groups = F,
           cell_size = 0.5)

plot_cells(glia.non, 
           color_cells_by = "hph",
           label_cell_groups = F,
           cell_size = 0.5)

colData(glia.non)$UMAP.cluster <- monocle3::clusters(glia.non, reduction_method = "UMAP")

plot_cells(glia.non,
           color_cells_by = "UMAP.cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(glia.non)$UMAP_1 <- reducedDims(glia.non)[["UMAP"]][,1]
colData(glia.non)$UMAP_2 <- reducedDims(glia.non)[["UMAP"]][,2]

plot.expr.UMAP(glia.non, "nlp-77", size = 0.5)
plot.expr.UMAP(glia.non, "grl-12", size = 0.5)
plot.expr.UMAP(glia.non, "D2092.8", size = 0.5)
plot.expr.UMAP(glia.non, "K01D12.5", size = 0.5)
plot.expr.UMAP(glia.non, "hlh-17", size = 0.5)
plot.expr.UMAP(glia.non, "mls-2", size = 0.5)
plot.expr.UMAP(glia.non, "wrt-6", size = 0.5)
plot.expr.UMAP(glia.non, "grl-18", size = 0.5)

glia.markers <- top_markers(glia.non)
glia.markers %>% filter(cell_group == 34) %>% arrange(desc(specificity)) %>% head(20)

colData(glia.non)$Cell.type <- ifelse(
  colData(glia.non)$UMAP.cluster %in% c(34),
  "ILso",
    as.character(colData(glia.non)$Cell.type)
  )

colData(glia.non)$Cell.type <- ifelse(
  colData(glia.non)$UMAP.cluster %in% c(6,31,24,20,35,21,19),
  "AMso",
  as.character(colData(glia.non)$Cell.type)
)

colData(glia.non)$Cell.type <- ifelse(
  colData(glia.non)$UMAP.cluster %in% c(23,10,4,29),
  "PHso",
  as.character(colData(glia.non)$Cell.type)
)

colData(glia.non)$Cell.type <- ifelse(
  colData(glia.non)$UMAP.cluster %in% c(38,17),
  "Glia_socket",
  as.character(colData(glia.non)$Cell.type)
)

# K09F5.6, hlh-17, 

colData(glia.non)$Cell.type <- ifelse(
  colData(glia.non)$UMAP.cluster %in% c(22),
  "CEPsh",
  as.character(colData(glia.non)$Cell.type)
)

colData(glia.non)$Cell.type <- ifelse(
  colData(glia.non)$UMAP.cluster %in% c(1,8),
  "Glia_sheath",
  as.character(colData(glia.non)$Cell.type)
)

colData(L1.clean.non)[colnames(glia.non),]$Cell.type <- colData(glia.non)$Cell.type
colData(L1.clean.filtered)[colnames(L1.clean.non),]$Cell.type <- colData(L1.clean.non)$Cell.type

## 053024

# I want to threshold the L1 neuronal data. The rest of the non-neuronal can wait.

# Loading in packages
setwd("./Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/")
source("./Monocle3Functions.txt")
library(monocle3)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(boot)
library(wbData)
gids <- wb_load_gene_ids("WS273")

L1 <- readRDS("./CeNGEN/L1_scRNAseq/L1_analysis/single_cell_datasets/042324_L1_all_live_cells_MeOH_new_saa_P_lineage_no_doublets_or_UPR_cds.rds")
L1.n <- L1[,colData(L1)$Tissue == "Neuron"]

