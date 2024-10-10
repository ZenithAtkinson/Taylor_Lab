# scRNA-Seq Analysis Pipeline for a single dataset

#If we load it in as embryo...
embryo <- readRDS("~/Lab_Data/Original_Files/cedata-Global datasetdowncell-2019-10-22.rds")

embryo@phenoData@data[1:5,1:15]
#head(embryo@meta.data[, 1:15], 5)
# The column we want is batch
table(embryo@phenoData@data$batch)
# We needto be able to split the 300, the 400, and the 500(1) and 500(2)
W300 = embryo[,embryo@phenoData@data$batch == "Waterston_300_minutes"]

# Subsetting for Waterson 400 minutes
#W400 <- subset(embryo, subset = batch == "Waterson_400_minutes")
# This gives us ~17000, which we can split into the different barcodes
table(W300@phenoData@data$batch)
head(W300@phenoData@data)

W300$Sample = substring(W300@phenoData@data$Cell, 24, 24)

table(W300$Sample)
# Once we have these character strings, they have a constant length so we can extract out the .1.1 or whatever easily

# Setting working directory to find the files
setwd("~/Dropbox (VU Basic Sciences)/Miller Lab/10X Genomics/")

# Loading in libraries
library(DropletUtils)
library(scater)
library(wbData)
library(dplyr)

# Getting a lookup table for switching between WBGeneIDs and gene names
gids <- wb_load_gene_ids("WS273") 

# Loading in the raw gene by barcode matrix as a SingleCellExperiment object
my_path <- "~/Lab_Data/Original_Files/Murray_b02/raw_feature_bc_matrix/"

# choose a name for the SingleCellExperiment object that makes sense
sce_walkthrough <- read10xCounts(my_path)

sce_walkthrough

# preview the barcodes as column data
head(colData(sce_walkthrough)$Barcode)

head(rowData(sce_walkthrough)) # ID, Symbol, Type
head(colData(sce_walkthrough)) # sample, barcode

# make the rownames the barcodes.
# The ID row is also the barcodes but this 
# adds a dimension that makes each row LABEL the barcode, which helps with the packaged functions
rownames(colData(sce_walkthrough)) <- colData(sce_walkthrough)$Barcode
colData(sce_walkthrough)
rownames(rowData(sce_walkthrough))
# For combining the sample with other datasets, I add the following information
# to the cell metadata table, including sex, genotype, experiment (strain info),
# developmental stage

colData(sce_walkthrough)$Batch <- "Murray_b02" # This can be the "batch" value (Murray_b02, 300, etc)
colData(sce_walkthrough)$Sample <- "" # This denotes the lane/channel that it came from (specificall for 300, 400, 500)
colData(sce_walkthrough)$Sex <- "Herm"
colData(sce_walkthrough)$Stage <- "Embryo"
colData(sce_walkthrough)$Genotype <- "wt"
head(colData(sce_walkthrough))

head(rowData(sce_walkthrough))

# cleaning up the gene metadata table
rowData(sce_walkthrough) <- rowData(sce_walkthrough)[,1:2]
colnames(rowData(sce_walkthrough)) <- c("id", "gene_short_name")
rowData(sce_walkthrough)$gene_short_name <- i2s(rowData(sce_walkthrough)$id, gids)
rowData(sce_walkthrough)

# plotting the barcode by UMI plot-  this can take a while
bcrank <- barcodeRanks(counts(sce_walkthrough))

## 
# highlight this whole section and use CMD+Enter to run all at once
options(bitmapType = "cairo")
given_lower <- 500

# Now try plotting again
plot(bcrank$rank, bcrank$total, log = "xy", xlab = "Rank", ylab = "Total UMI count")
o <- order(bcrank$rank)
lines(bcrank$rank[o], bcrank$fitted[o], col = "red")
abline(h = metadata(bcrank)$knee, col = "dodgerblue", lty = 2)
abline(h = metadata(bcrank)$inflection, col = "forestgreen", lty = 2)
abline(h = given_lower, col = "orange", lty = 2)
legend("bottomleft", lty = 2, col = c("dodgerblue", "forestgreen", "orange"),
       legend = c("knee", "inflection", "Given lower"))
# # 
# Knee and Inflec points (can also get given_lower? just set to what i want)
knee_point <- metadata(bcrank)$knee
inflect_point = metadata(bcrank)$inflection

# Running emptydrops to distinguish cell-containing droplets from empty droplets
# Using droplets with fewer than 50 UMIs as the background.
# Setting 50 is subjective, you can use the plot above to decide where the 
# dropoff is. Or you can set a hard threshold like below.
# Setting the seed ensures reproducibility

# FDR = False Discovery Rate

set.seed(100)
emptydrops.out <- emptyDrops(counts(sce_walkthrough), lower = inflect_point)
head(emptydrops.out)
emptydrops.out
summary(emptydrops.out$FDR, exclude = NULL)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0     0.0     1.0     0.7     1.0     1.0  1143248 

summary(emptydrops.out$FDR <= 0.01)
#  Mode   FALSE    TRUE    NA's
#logical  19278   53054 456764

head(emptydrops.out@metadata$ambient)

# Generating a dataframe of the genes in ambient RNA from the Emptydrops calculation
ambientRNA <- emptydrops.out@metadata$ambient
dim(ambientRNA) #19962 1
head(ambientRNA)

# Create a dataframe object containing all ambient genes and their expression values
ambientRNA.df <- data.frame(row.names = rownames(ambientRNA),
                            id = rownames(ambientRNA),
                            gene_name = i2s(rownames(ambientRNA), gids),
                            expression = ambientRNA[,1])

ambientRNA.df <- ambientRNA.df %>% arrange(desc(expression))
# View the highest expressed genes in the ambient profile
head(ambientRNA.df,15)
head(ambientRNA.df$gene_name,15)


# top 15 genes here
# "ctc-3" "nduo-6"   "ctc-2"  "atp-6"  "ctc-1"  "hil-7"   "his-24"  "ctb-1"   "nduo-1" 
# "hil-2"   "ZK380.6" "T24B8.3" "crt-1"   "cav-1"   "pat-10" 

# These are the genes that show up high on the list of ambient expression

# Other diagnostic tests

# T/F operator to determine which cells are real
is.cell <- emptydrops.out$FDR <= 0.01
#is.cell
sum(is.cell, na.rm = TRUE)

# 18173 claimed cells

table(Limited=emptydrops.out$Limited, Significant=is.cell)
# Significant
# Limited FALSE  TRUE
# FALSE   19278   41359
# TRUE      0    11695

# Because Limited == TRUE and Significant == FALSE, it means the number of permutations 
# was not limiting. -- Look for 0 in bottom left

# Statistics on total UMI counts for cells called as cells

# Total = UMIs in one cell

summary(emptydrops.out[which(emptydrops.out$FDR <= 0.01),]$Total)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51    1165    2019    2781    3861   37890 

# We can add the PValue and FDR for the calculation of whether a droplet was a cell or 
# not to the cell metadata table

sce_walkthrough$PValue <- emptydrops.out$PValue
head(colData(sce_walkthrough))

sce_walkthrough$FDR <- emptydrops.out$FDR
head(colData(sce_walkthrough))


# Keeping only the barcodes that correspond to cells
# apply filtering steps above to our sce object

# We already ran this line, but we do it again to keep things ordered
is.cell <- emptydrops.out$FDR <= 0.01

# filter step
sce_walkthrough <- sce_walkthrough[, which(is.cell), drop = F]
sce_walkthrough

# Generating a new filtered gene by barcode matrix for SoupX background RNA correction
#counts(sce_walkthrough)
sce_walkthrough.filt.counts <- counts(sce_walkthrough)
#sce_walkthrough.filt.counts

# Writing to a new folder for SoupX.
# Generating new filtered matrix files for SoupX

# to make this work, I create a "SoupX" subdirectory within the results folder.
# copy the raw_feature_bc_matrix folder into that folder, and make a new
# directory named "filtered_feature_bc_matrix"

# Copy and paste the "features.tsv.gz" file from the raw_feature_bc_matrix folder into 
# the new filtered_feature_bc_matrix folder. We will generate a new barcode file 
# and a filtered matrix file from the cells detected by Emptydrops.

# Need to now transfer the matrix.mtx and barcodes.tsv files from raw_feature_bc_matrix to the filtered_feature_bc_matrix
  # mv ~/Lab_Data/Original_Files/Murray_b02/raw_feature_bc_matrix/matrix.mtx ~/Lab_Data/SoupX/Murray_b02/filtered_feature_bc_matrix/
  # the .gz extension needs to be removed. Can be done with gzip -d *.gz (or does it? There is a gzip command below)
library(Matrix)
#counts
writeMM(sce_walkthrough.filt.counts, "~/Lab_Data/SoupX/Murray_b02/filtered_feature_bc_matrix/matrix.mtx")
system("gzip ~/Lab_Data/SoupX/Murray_b02/filtered_feature_bc_matrix/matrix.mtx")
writeMM(sce_walkthrough.filt.counts, "~/Lab_Data/SoupX/Murray_b02/filtered_feature_bc_matrix/matrix.mtx")

#barcodes
fil.barcodes <- colnames(sce_walkthrough)
write.table(fil.barcodes, "~/Lab_Data/SoupX/Murray_b02/filtered_feature_bc_matrix/barcodes.tsv", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
system("gzip ~/Lab_Data/SoupX/Murray_b02/filtered_feature_bc_matrix/barcodes.tsv")
write.table(fil.barcodes, "~/Lab_Data/SoupX/Murray_b02/filtered_feature_bc_matrix/barcodes.tsv", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx# Normalizing the data


# define per-cell size factors from the library sizes(i.e., total sum of counts per cell)
# normalize because there are a wide range of UMIs contained in each droplet (cell)
sizeFactors(sce_walkthrough) <- librarySizeFactors(sce_walkthrough)

# compute log-transformed normalized expression values from a count matrix in a SCE object
sce_walkthrough <- logNormCounts(sce_walkthrough)


# SoupX works best if you have some initial clustering data, we will do that 
# in Seurat

library(Seurat)
sce_walkthrough.S <- as.Seurat(sce_walkthrough)

# Generating initial clustering with Seurat to aid in SoupX background correction
# "Normalize the count data present in a given assay."
sce_walkthrough.S <- NormalizeData(sce_walkthrough.S)

# Identifies features that are outliers on a 'mean variability plot'
sce_walkthrough.S <- FindVariableFeatures(sce_walkthrough.S, selection.method = "vst", nfeatures = 4000)

# Scales and centers features in the dataset. 
# If variables are provided in vars.to.regress, they are individually regressed against each feature, 
# and the resulting residuals are then scaled and centered.
sce_walkthrough.S <- ScaleData(sce_walkthrough.S)


# the npcs term here determines how many PCs will be calculated. I like to include at 
# least 100. The Elbow plot code below will help set this
# AA: Takes a long time
sce_walkthrough.S <- RunPCA(sce_walkthrough.S, features = VariableFeatures(sce_walkthrough.S), npcs = 75)

# We use the elbox plot of principle components (x-axis) and the variance they explain 
# (y-axis) to pick the number of PCs (or dimensions) for UMAP. 

# I pick 5 or so PCs after the curve has flattened. If the curve has not flattened,
# you can rerun the RunPCA code and increase the npcs term.
ElbowPlot(sce_walkthrough.S, ndims = 50)

# Set the dims to the number of PCs selected from the ElbowPlot

# Initialize neighbor data and then cluster 
sce_walkthrough.S <- FindNeighbors(sce_walkthrough.S, dims = 1:75)
sce_walkthrough.S <- FindClusters(sce_walkthrough.S, resolution = 1.2)


# Running UMAP dimensionality reduction
sce_walkthrough.S <- RunUMAP(sce_walkthrough.S, dims = 1:75)

# plotting the UMAP, colored by clusters. Seurat starts numbering from 0, not 1.
DimPlot(sce_walkthrough.S, reduction = "umap")


# Seurat identified 43 clusters. There is bad separation between the clusters

# You can check the expression of a gene on the UMAP by using FeaturePlot. 
# It takes the WBGene id as input, here is sbt-1, a pan-neuronal gene
FeaturePlot(sce_walkthrough.S, features = "WBGene00011392")

# you can use a function from the wbData package to plot based on name
FeaturePlot(sce_walkthrough.S, features = s2i("sbt-1", gids))
FeaturePlot(sce_walkthrough.S, features = s2i("let-858", gids))
FeaturePlot(sce_walkthrough.S, features = s2i("bnc-1", gids))
FeaturePlot(sce_walkthrough.S, features = s2i("unc-4", gids))
FeaturePlot(sce_walkthrough.S, features = s2i("flp-1", gids))
FeaturePlot(sce_walkthrough.S, features = s2i("unc-122", gids))


# FeaturePlot(sce_walkthrough.S, features = s2i("symbol", gids))

# I typically will look for expression of genes that I know are specific to 
# see if the data contain neurons/cell types I'm interested in.

# I also check the expression of genes that are high in the ambient profile 
# to see if they have very high expression in subsets of cells with low expression 
# everywhere else (a sign of possible contamination)

# Getting the umap coordinates and groupings from the Seurat object to use in SoupX
umap.coords <- sce_walkthrough.S@reductions$umap
umap.coords

Embeddings(umap.coords)[1:3,1:2]

# create a new data frame containing the umap coordinates for each cell
sce_DimRed <- as.data.frame(Embeddings(umap.coords))

# Idents = cluster assigned by Seurat 
# We will add the cluster number corresponding to each coordinate pair
sce_DimRed$Cluster <- Idents(sce_walkthrough.S)
head(sce_DimRed)

# Now to run SoupX, with soupRange = c(0, 25). This uses barcodes with fewer than 
# 25 UMIs to calculate the background RNA ("soup")
# This background determination is similar to EmptyDrops, but we 
# are going to be less stringent since we are removing gene counts

library(SoupX)
# set a path to the directory you created for SoupX correction
my_dir_S <- "../scRNA_demonstration/7596-ST/SoupX/7596-ST-1"
scl.sce <- load10X(my_dir_S, soupRange = c(0, 25))


# looking at the first five rows and columns of the expression matrix
scl.sce$toc[1:5,1:5]

# listing the soup profile. this will be similar to the ambient profile calculated
# from Emptydrops
head(scl.sce$soupProfile[order(scl.sce$soupProfile$est, decreasing = TRUE), ], n = 30)


# adding gene names
scl.sce$soupProfile$gene_short_name <- i2s(rownames(scl.sce$soupProfile), gids)

head(scl.sce$soupProfile[order(scl.sce$soupProfile$est, decreasing = TRUE), ], n = 30)

# I will often copy to here the top genes for later reference
# top genes are ribosomal RNA (rrn-3.1), mitochondrial genes (nduo-6, atp-6, ctc genes), 
# neuronal (R102.2, flp-14, snet-1, sbt-1)

# you can add the UMAP coordinates to the soupX object so you don't have to type it 
# everytime
scl.sce <- setDR(scl.sce, sce_DimRed)

# Plotting the expression of some of these genes, using the Seurat UMAP coordinates.
# The TRUE cells (green circles) are cells SoupX calculates have real expression

plotMarkerMap(scl.sce, s2i("col-140", gids))
plotMarkerMap(scl.sce, s2i("nduo-6", gids))
plotMarkerMap(scl.sce, s2i("cpn-3", gids))

plotMarkerMap(scl.sce, s2i("flp-14", gids))
plotMarkerMap(scl.sce, s2i("clik-1", gids))
plotMarkerMap(scl.sce, s2i("flp-1", gids))

# Setting a list of genes that will be used to estimate contamination.
# I use a list of non-neuronal genes (if present), like muscle genes (pat-10, cpn-3,
# clik-1, myo-3) and hypodermal genes (col-140) and genes restricted to specific subsets
# of neurons (this depends on the subsets of cells targeted for sequencing)

est.cont <- s2i(c("pat-10", "cpn-3", "flp-1", "nlp-49", "pdf-1", "clik-1",
                  "snet-1"), gids)

# This code sets which cells will be used to estimate contamination 

useToEst = estimateNonExpressingCells(scl.sce, nonExpressedGeneList = list(sce = est.cont), 
                                      clusters = setNames(sce_DimRed$Cluster, rownames(sce_DimRed)))

# plots the cells used to estimate
plotMarkerMap(scl.sce, geneSet = est.cont, DR = sce_DimRed, useToEst = useToEst)

scl.sce <- setClusters(scl.sce, sce_DimRed$Cluster)

# Calculating the contamination fraction
scl.sce <- calculateContaminationFraction(scl.sce, list(sce = est.cont), useToEst = useToEst)
# Estimated global contamination fraction of 7.81%

# SoupX also has an automated contamination estimate. It is almost always less
# than what is calculated the manual way

scl.sce = autoEstCont(scl.sce)
# 5548 genes passed tf-idf cut-off and 3278 soup quantile filter.  Taking the top 100.
# Using 4993 independent estimates of rho.
# Estimated global rho of 0.01

# The manual method estimated a higher amount of contamination. I think I will go with 
# that for now.

# Re-running to make sure the estimate of contamination is the manual one.
scl.sce <- calculateContaminationFraction(scl.sce, list(sce = est.cont), useToEst = useToEst)

# Running the actual correction. the roundToInt term rounds the corrected counts 
# to integers, which is important for some downstream modeling
soup.out <- adjustCounts(scl.sce, roundToInt = TRUE)

# checking which genes were reduced the most by the correction. 
cntSoggy = rowSums(scl.sce$toc > 0)
cntStrained = rowSums(soup.out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed

i2s(names(mostZeroed),gids)

i2s(names(tail(sort(rowSums(scl.sce$toc > soup.out)/rowSums(scl.sce$toc > 0)), n = 20)),gids)



# you can plot the change in expression for contaminating and other genes
# Expression in white cells is left alone, others are corrected based on the 
# percent of soup estimates.

plotChangeMap(scl.sce, out, s2i("col-124", gids))
plotChangeMap(scl.sce, out, s2i("col-140", gids))
plotChangeMap(scl.sce, out, s2i("pat-10", gids))
plotChangeMap(scl.sce, out, s2i("nlp-12", gids))
plotChangeMap(scl.sce, out, s2i("flp-1", gids))
plotChangeMap(scl.sce, out, s2i("flp-12", gids))

# Saving the objects. Update the file path to your proper location
saveRDS(scl.sce, "../scRNA_demonstration/7596-ST/SoupX/7596-ST-2/021522_soupX_object.rds")
saveRDS(soup.out, "../scRNA_demonstration/7596-ST/SoupX/7596-ST-2/021522_soupX_corrected_matrix.rds")

# Observe reduced counts due to SoupX filter
soup.out[1:5,1:5]
scl.sce$toc[1:5,1:5]

# Creating a SingleCellExperiment object with the corrected expression data
sce_corrected <- SingleCellExperiment(assays = list(counts = soup.out), colData = colData(sce_walkthrough))  

# setting the gene metadata
rowData(sce_corrected) <- rowData(sce_walkthrough)
head(rowData(sce_corrected))
colData(sce_corrected)
# Using scater to calculate some quality control metrics based on mitochondrial genes

# and our empirically derived set of stress-responsive genes

mt.genes <- c("WBGene00010959", "WBGene00010961", "WBGene00010966", "WBGene00010963", "WBGene00010957", "WBGene00010967", "WBGene00010958", "WBGene00010960", "WBGene00010962","WBGene00000829","WBGene00010965","WBGene00010964","WBGene00014454","WBGene00014472")
stress.genes <- c("WBGene00009692", "WBGene00009691", "WBGene00002018", "WBGene00002016", "WBGene00002017",
                  "WBGene00006959", "WBGene00002026", "WBGene00004622", "WBGene00235164", "WBGene00004131",
                  "WBGene00012813", "WBGene00000097", "WBGene00015168", "WBGene00009180", "WBGene00013161",
                  "WBGene00016250", "WBGene00003903", "WBGene00008852", "WBGene00001496", "WBGene00005663",
                  "WBGene00014472", "WBGene00003148", "WBGene00004744", "WBGene00011300", "WBGene00011564",
                  "WBGene00010470", "WBGene00019006", "WBGene00009799", "WBGene00001031", "WBGene00007980",
                  "WBGene00019983", "WBGene00003954", "WBGene00008405", "WBGene00045366", "WBGene00019457",
                  "WBGene00000502", "WBGene00018748", "WBGene00012004", "WBGene00007514", "WBGene00011303",
                  "WBGene00015176", "WBGene00019664", "WBGene00201458", "WBGene00021945", "WBGene00198236",
                  "WBGene00002019", "WBGene00002020", "WBGene00196396", "WBGene00002015")

# From scater vignette
# Do for both uncorrected and corrected to compare difference
sce <- addPerCellQC(sce_walkthrough, subsets = list(Mito = mt.genes, stress = stress.genes))
sce_corrected <- addPerCellQC(sce_corrected, subsets = list(Mito = mt.genes, stress = stress.genes))

# Summarizing the percent of mitochondrial statistics across all cells in both 
# the uncorrected and corrected datasets
summary(sce$subsets_Mito_percent)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   1.152   1.961   3.144   3.585  55.709


summary(sce_corrected$subsets_Mito_percent)
# Min.  1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.3979  1.1574  2.4996  2.8571   58.0089 

# how many cells have <20% of UMIs from mitochondrial genes
table(sce$subsets_Mito_percent < 20)
# FALSE  TRUE 
# 181   17992 

table(sce_corrected$subsets_Mito_percent < 20)
# FALSE  TRUE 
# 189   17984

# Normalizing and saving the uncorrected and corrected datasets
sizeFactors(sce) <- librarySizeFactors(sce)
sce <- logNormCounts(sce)
sce
saveRDS(sce, "../scRNA_demonstration/7596-ST/7596-ST-1/121322_L2_herm_sce_uncorrected.rds")

sizeFactors(sce_corrected) <- librarySizeFactors(sce_corrected)
sce_corrected <- logNormCounts(sce_corrected)
sce_corrected
saveRDS(sce_corrected, "../scRNA_demonstration/7596-ST/7596-ST-1/121322_L2_herm_sce_SoupX_corrected.rds")


# Convert to a monocle3 object for annotations
wt.cds <- new_cell_data_set(counts(sce_corrected),
                            cell_metadata = colData(sce_corrected),
                            gene_metadata = rowData(sce_corrected))



# Now we can subset the dataset to remove cells with > 20% of reads from mitochondrial genes

# A cell data set can be subset using [ , ] following the name of the object. Information before 
# comma refers to the genes, whereas after the comma refers to the cells. 

# You can subset based on any logical test for either genes (using the rowData dataframe) or cells
# (using the colData(dataframe)). Here we keep only those cells with subsets_Mito_percent < 20.

wt.cds <- wt.cds[, colData(wt.cds)$subsets_Mito_percent < 20]

# Now we run an automated detection of genes in each cell, and then filter out genes that are detected
# in fewer than 6 single cells.
wt.cds <- detect_genes(wt.cds)
wt.cds <- wt.cds[rowData(wt.cds)$num_cells_expressed > 5,]

wt.cds
# 15673 genes in 17983 cells.

# we can summarize various metrics or tabulate the number of cells meeting certain conditions.
colnames(colData(wt.cds))
summary(colData(wt.cds)$num_genes_expressed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  38.0   492.0   724.0   766.8   984.0  3645.0 

table(colData(wt.cds)$num_genes_expressed > 50)
# 65 17918 

# Now we are going to run PCA, using the preprocess_cds command. We need to tell the function 
# how many PCs to include in the calculation.

wt.cds <- preprocess_cds(wt.cds, num_dim = 75)
plot_pc_variance_explained(wt.cds)

# Now we will run UMAP.
wt.cds <- reduce_dimension(wt.cds, 
                         reduction_method = "UMAP",
                         umap.min_dist = 0.3,
                         umap.n_neighbors = 75)

# Storing UMAP coordinates in the cell metadata for use with custom plots later.
colData(wt.cds)$UMAP_1 <- reducedDims(wt.cds)[["UMAP"]][,1]
colData(wt.cds)$UMAP_2 <- reducedDims(wt.cds)[["UMAP"]][,2]

# Identifying clusters. The "res" argument determines the resolution. Larger numbers will detect
# more clusters. 
wt.cds <- cluster_cells(wt.cds,
                      res = 3e-4)

# Plotting the cells, colored by cluster. You can color the cells in the plot by any column in 
# the colData dataframe.

plot_cells(wt.cds, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 0.5)

# plotting by genotype.
plot_cells(wt.cds,
           color_cells_by = "Genotype",
           label_groups_by_cluster = F,
           label_cell_groups = F,
           cell_size = 0.5)

# adding cluster assignments to the colData dataframe. They are stored elsewhere by default, but storing 
# them in the colData dataframe makes it easier to work with later.

colData(wt.cds)$cluster <- monocle3::clusters(wt.cds)

# Annotating cell types and tissue
colData(wt.cds)$Cell.type <- "Unannotated"

colData(wt.cds)$Tissue <- "Unannotated"

plot_cells(wt.cds,
           color_cells_by = "Cell.type",
           label_groups_by_cluster = F,
           group_label_size = 3,
           cell_size = 0.5)

# Adding some custom functions to the environment for the current session

source("./Monocle3Functions.txt")
plot.expr.UMAP(wt.cds, "unc-4", size = 0.5)
plot.expr.UMAP(wt.cds, "bnc-1", size = 0.5)

# Alternating back and forth between the cluster plot and expression plots can allow for 
# annotation of specific clusters

plot_cells(wt.cds, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 0.5)

# Assigning specific clusters an identity
colData(wt.cds)$Cell.type <- ifelse(
  colData(wt.cds)$cluster %in% c(15,9,32,50,46),
  "VA",
  as.character(colData(wt.cds)$Cell.type)
)

plot.expr.UMAP.L2.m3(wt.cds, "flp-1", size = 0.5)

colData(wt.cds)$Cell.type <- ifelse(
  colData(wt.cds)$cluster %in% c(38,44,2,6,18),
  "AVK",
  as.character(colData(wt.cds)$Cell.type)
)

plot.expr.UMAP.L2.m3(wt.cds, "pat-10", size = 0.5)

colData(wt.cds)$Cell.type <- ifelse(
  colData(wt.cds)$cluster %in% c(19),
  "body_wall_muscle",
  as.character(colData(wt.cds)$Cell.type)
)

plot.expr.UMAP.L2.m3(wt.cds, "aman-1", size = 0.5)

colData(wt.cds)$Cell.type <- ifelse(
  colData(wt.cds)$cluster %in% c(14,31),
  "Coelomocyte",
  as.character(colData(wt.cds)$Cell.type)
)

# Distinguishing neuronal clusters from non-neuronal with the pan-neuronal markers sbt-1, egl-3, egl-21
plot.expr.UMAP.L2.m3(wt.cds, "sbt-1", size = 0.5)
plot.expr.UMAP.L2.m3(wt.cds, "egl-3", size = 0.5)
plot.expr.UMAP.L2.m3(wt.cds, "egl-21", size = 0.5)

colData(wt.cds)$Tissue <- ifelse(
  colData(wt.cds)$cluster %in% c(14,31,19,24,22,43,16,30,20,55),
  "Non-neuronal",
  "Neuron"
)

# Identifying cluster-specific markers for each cluster to use for annotation,
wt.cds.markers <- top_markers(wt.cds, genes_to_test_per_group = 35,
                            marker_sig_test = F)
head(wt.cds.markers)
wt.cds.markers$gene_name <- i2s(wt.cds.markers$gene_id, gids)
wt.cds.markers <- as.data.frame(wt.cds.markers)


plot_cells(wt.cds, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 0.5)

# Checking cluster 8 against CeNGEN data
# This code will filter the the cluster marker table for cluster 8 and arrange by genes with 
# most restricted expression (highest specificity for the cluster)
wt.cds.markers %>% filter(cell_group == 8) %>% 
  arrange(desc(specificity)) %>%
  head(20) %>% select(gene_name)


# This gene list can be copied and pasted in the heatmap on the CeNGENApp to identify the 
# corresponding neuron.

colData(wt.cds)$Cell.type <- ifelse(
  colData(wt.cds)$cluster %in% c(8),
  "AIN",
  as.character(colData(wt.cds)$Cell.type)
)


#####
# Go through this process for each cluster. 

# Below, we will go through the pipeline again for an unc-4 dataset

# At the end, we combine this wt and unc-4 data sets and do some basic comparisons
#####


########################################################################

# Starting unc-4 preprocessing

########################################################################

# ALEX NOTE: From below here, removed all the code (unneeded I beleive)
