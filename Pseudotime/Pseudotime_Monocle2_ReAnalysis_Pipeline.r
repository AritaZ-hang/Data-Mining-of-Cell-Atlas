# load environment
options(stringAsFactors = F)
setwd('D:/Main')
rm(list = ls())
set.seed(0)
library(Seurat)
library(data.table)
library(SeuratWrappers)
library(patchwork)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(magrittr)
library(monocle)

# loading raw data
# Split the Blood cells into Myeloid, Erythroid & Lymphoid 3 categories according to the full annotation file.
if(F)
{
  Myeloid.cell.types = c('Neutrophil', 'Macrophage', 'Dendritic cell', 'Granulocyte','Mast cell', 'Hematopoietic Stem Progenitor cell')
  Erythroid.cell.types = c('Erythroid cell', 'Hematopoietic Stem Progenitor cell')
  Lymphoid.cell.types = c('NK cell', 'B cell', 'Marginal Zone B cell', 'Pre B cell', 'T cell', 'Cd4+ Cd8+ Double Positive T cell', 'Alpha-Beta T cell','Gamma-Delta T cell', 'Plasmacytoid Dendritic cell','Hematopoietic Stem Progenitor cell')
  load('D:/MCA_BatchRemove_dge/rmbatch_dge/temp_data/BoneMarrow_unintegrated.RData')
  load('D:/MCA_BatchRemove_dge/rmbatch_dge/temp_data/Spleen_uncombined.RData')
  load('D:/MCA_BatchRemove_dge/rmbatch_dge/temp_data/Thymus_unintegrated.RData')
  
  RestList = list(BoneMarrow2.filtered, BoneMarrow3.filtered, Spleen.filtered, Thymus1.filtered, Thymus2.filtered)
  rna_Bloodcells = merge(BoneMarrow1.filtered, RestList)
  cellinfo.full = fread('cellInfo.full.csv', header = T)[,2:3]
  metadata = rna_Bloodcells@meta.data
  metadata$cell = rownames(rna_Bloodcells@meta.data)
  metadata = inner_join(metadata, cellinfo.full, by = 'cell')
  rna_Bloodcells = AddMetaData(rna_Bloodcells, metadata = metadata$'Cell.type', col.name = 'Cell.type')
  Myeloid.rna = subset(rna_Bloodcells, `Cell.type` %in% Myeloid.cell.types)
  Erythroid.rna = subset(rna_Bloodcells, `Cell.type` %in% Erythroid.cell.types)
  Lymphoid.rna = subset(rna_Bloodcells, `Cell.type` %in% Lymphoid.cell.types)
  
  rm(rna_Bloodcells, BoneMarrow1.filtered, BoneMarrow2.filtered,BoneMarrow3.filtered, Thymus1.filtered, Thymus2.filtered, Spleen.filtered)
  gc()
  
}

# Aggregate Seurat Analysis functions to save time
SeuratAnalysis = function(object)
{
  object = NormalizeData(object, normalization.method = 'LogNormalize', scale.factor = 10000)
  object = FindVariableFeatures(object, selection.method = 'vst', nfeatures = 2000)
  object = ScaleData(object, rownames(object))
  object = RunPCA(object, features = VariableFeatures(object))
  object = FindNeighbors(object, dims = 1:20)
  object = FindClusters(object, resolution = 1)
  object = RunUMAP(object, dims = 1:20)
}

Myeloid.rna = SeuratAnalysis(Myeloid.rna)
Erythroid.rna = SeuratAnalysis(Erythroid.rna)
Lymphoid.rna = SeuratAnalysis(Lymphoid.rna)

# A simple function to plot Cell.type & Tissue
Dimplots = function(object, to_plot)
{
  colors = brewer.pal(10, 'Set3')
  gg = Seurat::DimPlot(object, group.by = to_plot, label = T,
               pt.size = 2, label.size = 3.5, repel = T, cols = colors) + NoLegend() +
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

SavePlots(Dimplots(Myeloid.rna, 'tissue'), 'Myeloid_rna_tissueplot.png', 5,5)
SavePlots(Dimplots(Myeloid.rna, 'Cell.type'), 'Myeloid_rna_annoplot.png', 5,5)
SavePlots(Dimplots(Erythroid.rna, 'tissue'), 'Erythroid_rna_tissueplot.png', 5,5)
SavePlots(Dimplots(Erythroid.rna, 'Cell.type'), 'Erythroid_rna_annoplot.png', 5,5)
SavePlots(Dimplots(Lymphoid.rna, 'tissue'), 'Lymphoid_rna_tissueplot.png', 5,5)
SavePlots(Dimplots(Lymphoid.rna, 'Cell.type'), 'Lymphoid_rna_annoplot.png',5,5)

# Using the data in seurat object to create monocle cds
# @params: object, a seurat scRNA-seq object
# return with a monocle cds object
CreateMonocleCds = function(object)
{
  matrix = as(as.matrix(GetAssayData(object, slot = 'counts')), 'sparseMatrix')
  feature_ann = data.frame(gene_id = rownames(matrix), gene_short_name = rownames(matrix))
  rownames(feature_ann) = rownames(matrix)
  fd = new('AnnotatedDataFrame', data = feature_ann)
  sample_ann = object@meta.data
  pd = new('AnnotatedDataFrame', data = sample_ann)
  cds = newCellDataSet(matrix, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
  return(cds)
}
Myeloid.cds = CreateMonocleCds(Myeloid.rna)
Erythroid.cds = CreateMonocleCds(Erythroid.rna)
Lymphoid.cds = CreateMonocleCds(Lymphoid.rna)

# Aggregate Monocle Pipeline in order to save some time
# @params: cds, the monocle cds object
# @params:object, a seurat scRNA-seq object
MonoclePipeline = function(cds, object)
{
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  markers = VariableFeatures(object)
  cds = setOrderingFilter(cds, ordering_genes = markers)
  cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds = orderCells(cds)
  return(cds)
}
Myeloid.cds = MonoclePipeline(Myeloid.cds, Myeloid.rna)
Erythroid.cds = MonoclePipeline(Erythroid.cds, Erythroid.rna)
Lymphoid.cds = MonoclePipeline(Lymphoid.cds, Lymphoid.rna)

# adjusting Monocle Trajectorys
# Well it does not suit functional styles...
GM_state = function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Cell.type)[,"Hematopoietic Stem Progenitor cell"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

if(F)
{
  plot_cell_trajectory(Myeloid.cds, color_by = 'Pseudotime')
  Myeloid.cds = orderCells(Myeloid.cds, root_state = GM_state(Myeloid.cds))
  Erythroid.cds = orderCells(Erythroid.cds, root_state = GM_state(Erythroid.cds))
  Lymphoid.cds = orderCells(Lymphoid.cds, root_state = GM_state(Lymphoid.cds))
  plot_cell_trajectory(Erythroid.cds, color_by = 'Pseudotime')
  
  plot_cell_trajectory(Lymphoid.cds, color_by = 'Cell.type') + theme(legend.position = 'right', legend.text = element_text(size = 8), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8), legend.title = element_text(size = 10), axis.ticks = element_blank(), axis.text = element_blank())
  plot_cell_trajectory(Lymphoid.cds, color_by = 'Pseudotime')
}



# BEAM & Branch heatmap plot 
# the operating time is restricted to cpu cores
CreateBEAMObject = function(cds, branch_point)
{
  BEAM_res = BEAM(cds, branch_point = branch_point, cores = 6)
  BEAM_res = BEAM_res[, c('gene_short_name', 'pval', 'qval')]
  return(BEAM_res)
}
# branch_point
BEAM_res.Myeloid = CreateBEAMObject(Myeloid.cds, branch_point = 1)
BEAM_res.Erythroid = CreateBEAMObject(Erythroid.cds, branch_point = 1)
BEAM_res.Lymphoid = CreateBEAMObject(Lymphoid.cds, branch_point = 2)

# save files
if(F)
{
  saveRDS(BEAM_res.Erythroid, file = 'BEAM_res.Erythroid.rds')
  saveRDS(BEAM_res.Lymphoid, file = 'BEAM_res.Lymphoid.rds')
  saveRDS(BEAM_res.Myeloid, file = 'BEAM_res.Myeloid.rds')
  
  save(Erythroid.cds, Lymphoid.cds, Myeloid.cds, file = 'Monocle2_ReAnalysis_cds.RData')
}

# plotting branchplot
if(F)
{
Myeloid.branchplot = plot_genes_branched_heatmap(Myeloid.cds[row.names(subset(BEAM_res.Myeloid, qval < 1e-4)),],
                                                 branch_point = 1,
                                                 num_clusters = 5,# how many clusters
                                                 cores = 6,
                                                 branch_labels = c('Cell fate 1', 'Cell fate 2'),
                                                 branch_colors = c('#979797','#F05662','#7990C8'), use_gene_short_name = T, show_rownames = F, return_heatmap = T)
pdf('rna_myeloid_branched_heatmap.pdf', width = 7, height = 8)
Myeloid.branchplot$ph_res
dev.off()

Erythroid.branchplot = plot_genes_branched_heatmap(Erythroid.cds[row.names(subset(BEAM_res.Erythroid, qval < 1e-4)),],
                                                 branch_point = 1,
                                                 num_clusters = 5,# how many clusters
                                                 cores = 6,
                                                 branch_labels = c('Cell fate 1', 'Cell fate 2'),
                                                 branch_colors = c('#979797','#F05662','#7990C8'), use_gene_short_name = T, show_rownames = F, return_heatmap = T)
pdf('rna_Erythroid_branched_heatmap.pdf', width = 7, height = 8)
Erythroid.branchplot$ph_res
dev.off()

Lymphoid.branchplot = plot_genes_branched_heatmap(Lymphoid.cds[row.names(subset(BEAM_res.Lymphoid, qval < 1e-4)),],
                                                   branch_point = 1,
                                                   num_clusters = 5,# how many clusters
                                                   cores = 6,
                                                   branch_labels = c('Cell fate 1', 'Cell fate 2'),
                                                   branch_colors = c('#979797','#F05662','#7990C8'), use_gene_short_name = T, show_rownames = F, return_heatmap = T)
pdf('rna_Lymphoid_branched_heatmap.pdf', width = 7, height = 8)
Lymphoid.branchplot$ph_res
dev.off()
}

Myeloid.clustered.genes = Myeloid.branchplot$annotation_row
Erythroid.clustered.genes = Erythroid.branchplot$annotation_row
Lymphoid.clustered.genes = Lymphoid.branchplot$annotation_row

write.csv(Myeloid.clustered.genes, 'Myeloid_rna_clustered_genes.csv')
write.csv(Erythroid.clustered.genes, 'Erythroid_rna_clustered_genes.csv')
write.csv(Lymphoid.clustered.genes, 'Lymphoid_rna_clustered_genes.csv')

# and these genes will undergo Gene Ontology enrichment

# save the cds pdata
Myeloid.metadata = pData(Myeloid.cds)
Erythroid.metadata = pData(Erythroid.cds)
Lymphoid.metadata = pData(Lymphoid.cds)

save(Myeloid.metadata, Erythroid.metadata, Lymphoid.metadata, file = 'Pseudotime_RNA_meta.RData')

# and now it's the end

save(Myeloid.rna, Erythroid.rna, Lymphoid.rna, file = 'Pseudotime_RNA_files.RData')