# load environment
setwd('D:/Main')
options(stringAsFactors = F)
set.seed(0)
rm(list = ls())
library(Seurat)
library(Signac)
library(dplyr)
library(magrittr)
library(Matrix)
library(data.table)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)

if(F)
{
  Myeloid.cell.types = c('Neutrophil', 'Macrophage', 'Dendritic cell', 'Granulocyte',
                         'Mast cell', 'Hematopoietic Stem Progenitor cell')
  Erythroid.cell.types = c('Erythroid cell', 'Hematopoietic Stem Progenitor cell')
  Lymphoid.cell.types = c('NK cell', 'B cell', 'Marginal Zone B cell', 'Pre B cell', 
                          'T cell', 'Cd4+ Cd8+ Double Positive T cell', 'Alpha-Beta T cell','Gamma-Delta T cell', 'Plasmacytoid Dendritic cell','Hematopoietic Stem Progenitor cell')
  
  # Process the raw data
  atac_matrix = readMM('D:/ATAC/atac_site_filtered_del_PfcC.transfer.mtx')
  atac_matrix = atac_matrix * 1
  regions_transfer = read.table('D:/ATAC/regions_transfer.csv')
  barcodes_transfer = read.table('D:/ATAC/barcodes_transfer.csv')
  colnames(atac_matrix) = barcodes_transfer[,1]
  rownames(atac_matrix) = regions_transfer[,1]
  
  scvi.annotation.full = readRDS('scvi_annotation_refined.rds')
  Myeloid.cells = scvi.annotation.full[scvi.annotation.full$refined.anno %in% Myeloid.cell.types,]$cell
  Lymphoid.cells = scvi.annotation.full[scvi.annotation.full$refined.anno %in% Lymphoid.cell.types,]$cell
  Erythroid.cells = scvi.annotation.full[scvi.annotation.full$refined.anno %in% Erythroid.cell.types,]$cell
  
  Myeloid.mtx = atac_matrix[, Myeloid.cells]
  Lymphoid.mtx = atac_matrix[, Lymphoid.cells]
  Erythroid.mtx = atac_matrix[, Erythroid.cells]
}

annotations = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) = 'UCSC'
genome(annotations) = 'mm10'

# Create Signac seurat object and perform some analysis pipelines
SignacAnalysis = function(mtx, annotations)
{
  atac_assay = CreateChromatinAssay(mtx, min.cells = 10, genome = 'mm10', annotation = annotations, sep = c('_','_'))
  atac_object = CreateSeuratObject(atac_assay, assay = 'peaks', project = 'ATAC')
  atac_object = RunTFIDF(atac_object)
  atac_object = FindTopFeatures(atac_object, min.cutoff = 'q0')
  atac_object = RunSVD(atac_object)
  atac_object = RunUMAP(atac_object, reduction = 'lsi', dims = 2:30)
  atac_object = FindNeighbors(atac_object, reduction = 'lsi', dims = 2:30)
  atac_object = FindClusters(atac_object, algorithm = 3, resolution = 1)
  return(atac_object)
}

Myeloid.atac = SignacAnalysis(Myeloid.mtx, annotations)
Erythroid.atac = SignacAnalysis(Erythroid.mtx,annotations)
Lymphoid.atac = SignacAnalysis(Lymphoid.mtx, annotations)

# Adding Metadata to atac object
# @params: target_types: Myeloid.cell.types/Erythroid.cell.types/Lymphoid.cell.types
# @anno.file = scvi.annotation.full, key for annotation is refined.anno
AnnotationAddition = function(object, target_types, anno.file)
{
  anno.file = anno.file[anno.file$refined.anno %in% target_types,]
  metadata = object@meta.data
  metadata$cell = rownames(metadata)
  metadata = inner_join(metadata, anno.file, by = 'cell')
  object = AddMetaData(object, metadata = metadata$refined.anno, col.name = 'Cell.type')
  return(object)
}
Myeloid.atac = AnnotationAddition(Myeloid.atac, Myeloid.cell.types, scvi.annotation.full)
Erythroid.atac = AnnotationAddition(Erythroid.atac, Erythroid.cell.types, scvi.annotation.full)
Lymphoid.atac = AnnotationAddition(Lymphoid.atac, Lymphoid.cell.types, scvi.annotation.full)

# A simple function to plot Cell.type & Tissue
Dimplots = function(object, to_plot)
{
  colors = brewer.pal(10, 'Set3')
  gg = Seurat::DimPlot(object, group.by = to_plot, label = T,
                       pt.size = 1, label.size = 3.5, repel = T, cols = colors) + NoLegend() +
    ggtitle('') + 
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(fill = 'NA', color = 'black',
                                      size = 1, linetype = 'solid')) + 
    labs(x = 'UMAP1', y = 'UMAP2')
  return(gg)
}
# A simple function to save plots
SavePlots = function(gg, tosaveName, width, height)
{
  filename = tosaveName
  ggsave(gg, filename = filename, width = width, height = height)
}

SavePlots(Dimplots(Myeloid.atac, 'Cell.type'), 'Myeloid_atac_annoplot.png', 5,5)
SavePlots(Dimplots(Erythroid.atac, 'Cell.type'), 'Erythroid_atac_annoplot.png', 5,5)
SavePlots(Dimplots(Lymphoid.atac, 'Cell.type'), 'Lymphoid_atac_annoplot.png',5,5)

save(Myeloid.atac, Erythroid.atac, Lymphoid.atac, file = 'Pseudotime_ATAC_files.RData')
# Gene activity will be added in the FNN analysis file.