# Code used in the study: "Continuous clonal labeling reveals uniform progenitor potential in the adult exocrine pancreas", Sophie C. Lodestijn & Tom van den Bosch, et al.
# Please direct any questions or comments related to this code to Leandro Moreno: l.ferreiramoreno@amsterdamumc.nl

library('Seurat')
library('dplyr')
library('ggplot2')
library('monocle')
library("ggpubr")

#load and process run1 EYFP-
run1_YFPminus <- Read10X(data.dir ="run1_YFPminus/filtered_feature_bc_matrix/")
run1_YFPminus_object= CreateSeuratObject(counts = run1_YFPminus)
run1_YFPminus_object[["percent.mt"]] <- PercentageFeatureSet(run1_YFPminus_object, pattern = "^mt-")
run1_YFPminus_object[["Batch"]] <- "run1"
run1_YFPminus_object[["yfp"]] <- "yfp-"

run1_YFPminus_object_subset <- subset(run1_YFPminus_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
run1_YFPminus_object_subset= SCTransform(object = run1_YFPminus_object_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))

run1_YFPminus_object_subset <- RunPCA(run1_YFPminus_object_subset, verbose = FALSE)
run1_YFPminus_object_subset <- RunUMAP(run1_YFPminus_object_subset, dims = 1:10, verbose = FALSE)
run1_YFPminus_object_subset <- FindNeighbors(run1_YFPminus_object_subset, dims = 1:30, verbose = FALSE)
run1_YFPminus_object_subset <- FindClusters(run1_YFPminus_object_subset, verbose = FALSE, resolution = 0.6)

#load and process run1 EYFP+
run1_YFPplus <- Read10X(data.dir = "run1_YFPplus/filtered_feature_bc_matrix/")
run1_YFPplus_object= CreateSeuratObject(counts = run1_YFPplus)
run1_YFPplus_object[["percent.mt"]] <- PercentageFeatureSet(run1_YFPplus_object, pattern = "^mt-")
run1_YFPplus_object[["Batch"]] <- "run1"
run1_YFPplus_object[["yfp"]] <- "yfp+"

run1_YFPplus_object_subset <- subset(run1_YFPplus_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
run1_YFPplus_object_subset= SCTransform(object = run1_YFPplus_object_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))

run1_YFPplus_object_subset <- RunPCA(run1_YFPplus_object_subset, verbose = FALSE)
run1_YFPplus_object_subset <- RunUMAP(run1_YFPplus_object_subset, dims = 1:10, verbose = FALSE)
run1_YFPplus_object_subset <- FindNeighbors(run1_YFPplus_object_subset, dims = 1:30, verbose = FALSE)
run1_YFPplus_object_subset <- FindClusters(run1_YFPplus_object_subset, verbose = FALSE, resolution = 0.9)

#load and process run2 EYFP-
run2_YFPminus <- Read10X(data.dir ="run2_YFPminus/filtered_feature_bc_matrix/")
run2_YFPminus_object= CreateSeuratObject(counts = run2_YFPminus)
run2_YFPminus_object[["percent.mt"]] <- PercentageFeatureSet(run2_YFPminus_object, pattern = "^mt-")
run2_YFPminus_object[["Batch"]] <- "run2"
run2_YFPminus_object[["yfp"]] <- "yfp_minus"

run2_YFPminus_object_subset <- subset(run2_YFPminus_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
run2_YFPminus_object_subset= SCTransform(object = run2_YFPminus_object_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))

run2_YFPminus_object_subset <- RunPCA(run2_YFPminus_object_subset, verbose = FALSE)
run2_YFPminus_object_subset <- RunUMAP(run2_YFPminus_object_subset, dims = 1:30, verbose = FALSE)
run2_YFPminus_object_subset <- FindNeighbors(run2_YFPminus_object_subset, dims = 1:30, verbose = FALSE)
run2_YFPminus_object_subset <- FindClusters(run2_YFPminus_object_subset, verbose = FALSE, resolution = 0.6)

#load and process run2 EYFP+
run2_YFPplus <- Read10X(data.dir = "run2_YFPplus/filtered_feature_bc_matrix/")
run2_YFPplus_object= CreateSeuratObject(counts = run2_YFPplus)
run2_YFPplus_object[["percent.mt"]] <- PercentageFeatureSet(run2_YFPplus_object, pattern = "^mt-")
run2_YFPplus_object[["Batch"]] <- "run2"
run2_YFPplus_object[["yfp"]] <- "yfp+"

run2_YFPplus_object_subset <- subset(run2_YFPplus_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
run2_YFPplus_object_subset= SCTransform(object = run2_YFPplus_object_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))

run2_YFPplus_object_subset <- RunPCA(run2_YFPplus_object_subset, verbose = FALSE)
run2_YFPplus_object_subset <- RunUMAP(run2_YFPplus_object_subset, dims = 1:10, verbose = FALSE)
run2_YFPplus_object_subset <- FindNeighbors(run2_YFPplus_object_subset, dims = 1:30, verbose = FALSE)
run2_YFPplus_object_subset <- FindClusters(run2_YFPplus_object_subset, verbose = FALSE, resolution = 0.9)

#load and process run3 EYFP+
run3_YFPplus <- Read10X(data.dir ="run3_YFPplus/filtered_feature_bc_matrix/")
run3_YFPplus_object= CreateSeuratObject(counts = run3_YFPplus)
run3_YFPplus_object[["percent.mt"]] <- PercentageFeatureSet(run3_YFPplus_object, pattern = "^mt-")
run3_YFPplus_object[["Batch"]] <- "run3"
run3_YFPplus_object[["yfp"]] <- "yfp+"

run3_YFPplus_object_subset <- subset(run3_YFPplus_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
run3_YFPplus_object_subset= SCTransform(object = run3_YFPplus_object_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))

run3_YFPplus_object_subset <- RunPCA(run3_YFPplus_object_subset, verbose = FALSE)
run3_YFPplus_object_subset <- RunUMAP(pancreatitis_plus_object_subset, dims = 1:30, verbose = FALSE)
run3_YFPplus_object_subset <- FindNeighbors(run3_YFPplus_object_subset, dims = 1:30, verbose = FALSE)
run3_YFPplus_object_subset <- FindClusters(run3_YFPplus_object_subset, verbose = FALSE, resolution = 0.6)

#load and process run3 EYFP-
run3_YFPminus <- Read10X(data.dir = "run3_YFPminus/filtered_feature_bc_matrix/")
run3_YFPminus_object= CreateSeuratObject(counts = run3_YFPminus)
run3_YFPminus_object[["percent.mt"]] <- PercentageFeatureSet(run3_YFPminus_object, pattern = "^mt-")
run3_YFPminus_object[["Batch"]] <- "run3"
run3_YFPminus_object[["yfp"]] <- "yfp-"

run3_YFPminus_object_subset <- subset(run3_YFPminus_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
run3_YFPminus_object_subset= SCTransform(object = run3_YFPminus_object_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))

run3_YFPminus_object_subset <- RunPCA(run3_YFPminus_object_subset, verbose = FALSE)
run3_YFPminus_object_subset <- RunUMAP(run3_YFPminus_object_subset, dims = 1:30, verbose = FALSE)
run3_YFPminus_object_subset <- FindNeighbors(run3_YFPminus_object_subset, dims = 1:30, verbose = FALSE)
run3_YFPminus_object_subset <- FindClusters(run3_YFPminus_object_subset, verbose = FALSE, resolution = 0.9)

#Merge run 1 and run2  - No batch effect detected
run1_run2_subset=merge(x = run1_YFPminus_object_subset,y = c(run1_YFPplus_object_subset, run2_YFPplus_object_subset, run2_YFPminus_object_subset))

run1_run2_subset= SCTransform(object = run1_run2_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))

run1_run2_subset <- RunPCA(run1_run2_subset, verbose = FALSE)
run1_run2_subset <- RunUMAP(run1_run2_subset, dims = 1:30, verbose = FALSE, metric = "euclidean")
run1_run2_subset <- FindNeighbors(run1_run2_subset, dims = 1:30, verbose = FALSE)
run1_run2_subset <- FindClusters(run1_run2_subset, verbose = FALSE, resolution = 0.8)

DimPlot(run1_run2_subset, label = TRUE, group.by="Batch")

#merge EYFP- and EYFP- from run3
run3_subset=merge(x = run3_YFPplus_object_subset,y = run3_YFPminus_object_subset)

run3_subset= SCTransform(object = run3_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))
run3_subset <- RunPCA(run3_subset, verbose = FALSE)
run3_subset <- RunUMAP(run3_subset, dims = 1:30, verbose = FALSE, metric = "euclidean")
run3_subset <- FindNeighbors(run3_subset, dims = 1:30, verbose = FALSE)
run3_subset <- FindClusters(run3_subset, verbose = FALSE, resolution = 0.9)

#Merge all the runs
run1_run2_run3_subset=merge(x = run3_subset,y = run1_run2_subset)

run1_run2_run3_subset= SCTransform(object = run1_run2_run3_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))
run1_run2_run3_subset <- RunPCA(run1_run2_run3_subset, verbose = FALSE)
run1_run2_run3_subset <- RunUMAP(run1_run2_run3_subset, dims = 1:30, verbose = FALSE, metric = "euclidean")
run1_run2_run3_subset <- FindNeighbors(run1_run2_run3_subset, dims = 1:30, verbose = FALSE)
run1_run2_run3_subset <- FindClusters(run1_run2_run3_subset, verbose = FALSE, resolution = 0.6)

#Confirm the presence of Batch effect between run1,2 vc run3
DimPlot(run1_run2_run3_subset, label = TRUE, group.by="Batch") + NoLegend()

#correction of Batch effet
pancreas.list <- SplitObject(run1_run2_run3_subset, split.by = "Batch")

for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
                                             nfeatures = 5000, verbose = FALSE)
}

reference.list <- pancreas.list[c("run1", "run2", "run3")]

pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30,  anchor.features = 5000) #K.filter was necessary due to the low number of cells in at least 1 dataset
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated_scaled <- ScaleData(pancreas.integrated, verbose = FALSE)
#pancreas.integrated_scaled= SCTransform(object = pancreas.integrated,vars.to.regress = c("nFeature_RNA", "percent.mt"))
pancreas.integrated_scaled <- RunPCA(pancreas.integrated_scaled, verbose = FALSE)
pancreas.integrated_scaled <- RunUMAP(pancreas.integrated_scaled, dims = 1:30, verbose = FALSE)
pancreas.integrated_scaled <- FindNeighbors(pancreas.integrated_scaled, dims = 1:30, verbose = FALSE)
pancreas.integrated_scaled <- FindClusters(pancreas.integrated_scaled, verbose = FALSE, resolution = 0.9)

#Confirm that the data is free of Batch effect
DimPlot(pancreas.integrated_scaled, label = TRUE, group.by="Batch")

DimPlot(pancreas.integrated_scaled, label = TRUE,group.by="yfp",pt.size = 0.3, cols = c("#666666", "#FFCC29")) + ylim (c(-15,15)) + xlim(c(-15,15))
