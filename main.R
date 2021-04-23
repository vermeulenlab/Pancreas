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


#Merge run 1 and run2  - No batch effect detected
run1_run2_subset=merge(x = run1_YFPminus_object_subset,y = c(run1_YFPplus_object_subset, run2_YFPplus_object_subset, run2_YFPminus_object_subset))

run1_run2_subset= SCTransform(object = run1_run2_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))

run1_run2_subset <- RunPCA(run1_run2_subset, verbose = FALSE)
run1_run2_subset <- RunUMAP(run1_run2_subset, dims = 1:30, verbose = FALSE, metric = "euclidean")
run1_run2_subset <- FindNeighbors(run1_run2_subset, dims = 1:30, verbose = FALSE)
run1_run2_subset <- FindClusters(run1_run2_subset, verbose = FALSE, resolution = 0.8)

DimPlot(run1_run2_subset, label = TRUE, group.by="Batch")
