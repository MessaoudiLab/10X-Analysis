#Use the *_filtered_feature_bc_matrix as input for Seurat for all downstream purposes
#For additional help, refer to Hemberg lab's tutorial: https://broadinstitute.github.io/2019_scWorkshop/data-wrangling-scrnaseq.html

library(Seurat)
library(cowplot)
library(dplyr)
library(Matrix)
library(gdata)

***********************************************************************************************************************************
#Step 1: Reading in the matrix
#Seurat function to read 10X data
leans_villous <- Read10X(data.dir = "Lean_finalrun/outs/filtered_feature_bc_matrix")) 


#Check the dimensions of the file
dim(leans_villous)

***********************************************************************************************************************************
##Step 2: Filtering low quality cells. 
##Droplets with low gene counts are potentially ambient RNA or failed libraries
##Droplets with unusually high gene/RNA counts are potential doublets.
##Summary Plots

counts_per_cell <- Matrix::colSums(leans_villous)
counts_per_gene <- Matrix::rowSums(leans_villous)

#Expressed genes per cell or 'complexity' of cells
genes_per_cell <- Matrix::colSums(leans_villous > 0)
#Cells expressing genes or 'rarity' of cells
cells_per_gene <- Matrix::rowSums(leans_villous > 0)
pdf("Leans_villous_counts_per_cell.pdf")
hist(log10(counts_per_cell+1),main='Counts Per Cell',col='wheat')
dev.off()

pdf("Leans_villous_counts_per_gene.pdf")
hist(log10(counts_per_gene+1),main='Counts Per Gene',col='wheat')
dev.off()

pdf("Leans_villous_genes_per_cell.pdf")
hist(log10(genes_per_cell+1),main='Counts Per Gene',col='wheat')
dev.off()

pdf("Leans_villous_cells_per_gene.pdf")
hist(log10(cells_per_gene+1),main='Counts Per Gene',col='wheat')
dev.off()

#Plot cells ranked by their number of detected genes
pdf("Leans_villous_genes_per_cell_ranked.pdf")
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
dev.off()

#Filtering data based on quality - min.cells: a gene expressed in at least these many cells; min.genes: a cell with at least these many genes
#Numbers are somewhat arbitrary. Please read papers from your field. Using settings from first trimester paper
#Data input is unnormalized raw counts.
#All downstream analyses will be of object Seurat. 

countsSeurat_lean_villous <-CreateSeuratObject(counts = counts, min.cells = 3, min.features = 500)

#Check the number of cells that pass the threshold
countsSeurat_lean_villous

#Additional filtering of damaged cells is using a metric measuring %reads coming from mitochondrial genes. Relative enrichment of mtDNA
#is indicative of leaky cells and cell death.

#Calculcate the mitochondrial and add metadata to the object
countsSeurat_lean_villous[["percent.mt"]] <- PercentageFeatureSet(countsSeurat_lean_villous, pattern = "^MT-")

#To view, add, and manipulate metadata in the object use the following command
head(countsSeurat@meta.data, 5)

pdf("Lean_villous_mtPlots.pdf")
VlnPlot(countsSeurat_lean_villous, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#Use the following plot to take a call on the gene and mito cutoffs for downstream QC
plot1 <- FeatureScatter(countsSeurat_lean_villous, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(countsSeurat_lean_villous, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("Lean_villous_nFeature_RNA.pdf")
print(plot2)
dev.off()

pdf("Lean_villous_percentMT.pdf")
print(plot1)
dev.off()

#Assessing housekeeping genes
hkgenes <- read.table("../../code/housekeeping.txt", skip = 2)
hkgenes <- as.vector(hkgenes$V1)
hkgenes.found <- which(toupper(rownames(GetAssayData(countsSeurat_lean_villous))) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(GetAssayData(countsSeurat_lean_villous)[hkgenes.found, ] > 0)
countsSeurat_lean_villous <- AddMetaData(object = countsSeurat_lean_villous, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")

#Assessing cell cycle genes

#Calculate dissociation scores


#Filter the low quality cells
countsSeurat_filtered <- subset(countsSeurat, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 20 & n.exp.hkgenes > 25)

#Save the R object for future retrieval.
saveRDS(leans_villous_filtered,file="Term Landscape/leans_villous_filted.Robj")

*************************************************************************************************************************************************
