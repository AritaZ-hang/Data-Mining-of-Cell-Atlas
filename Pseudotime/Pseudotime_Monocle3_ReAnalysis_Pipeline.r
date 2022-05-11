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
library(monocle3)

load('Pseudotime_RNA_files.RData')

#Create Monocle3 cds file
MonoclePipeline = function(object)
{
  cds = as.cell_data_set(object)
  cds = estimate_size_factors(cds)
  cds = cluster_cells(cds = cds)
  cds = learn_graph(cds, use_partition = T)
  root_cells = colnames(object[, object$Cell.type == 'Hematopoietic Stem Progenitor cell'])
  cds = order_cells(cds, reduction_method = 'UMAP', root_cells = root_cells)
  #cds@rowRanges@elementMetaData@listData[['gene_short_name']] = rownames(object[['RNA']])
  return(cds)
}
Myeloid.rna.cds = MonoclePipeline(Myeloid.rna)
Erythroid.rna.cds = MonoclePipeline(Erythroid.rna)
Lymphoid.rna.cds = MonoclePipeline(Lymphoid.rna)

Myeloid.rna.cds@rowRanges@elementMetadata@listData[['gene_short_name']] = rownames(Myeloid.rna)
Erythroid.rna.cds@rowRanges@elementMetadata@listData[['gene_short_name']] = rownames(Erythroid.rna)
Lymphoid.rna.cds@rowRanges@elementMetadata@listData[['gene_short_name']] = rownames(Lymphoid.rna)

# Plotting Monocle3 Trajactories
TrajactoryPlot = function(cds)
{
  gg = plot_cells(cds, color_cells_by = 'pseudotime', cell_size = 1, label_branch_points = T, label_principal_points = F, label_leaves = F, label_roots = T, trajectory_graph_color = 'white')
  return(gg)
}
Myeloid.rna.trajactory = TrajactoryPlot(Myeloid.rna.cds)
Erythroid.rna.trajactory = TrajactoryPlot(Erythroid.rna.cds)
Lymphoid.rna.trajectory = TrajactoryPlot(Lymphoid.rna.cds)

# Gene co-expression analysis by Monocle3
MonocleModule = function(cds)
{
  Track_res = graph_test(cds, neighbor_graph = 'principal_graph', cores = 1)
  Track_res.ids = row.names(subset(Track_res[order(Track_res$morans_I, decreasing = T),],q_value < 0.05))
  module_df = find_gene_modules(cds[Track_res.ids,], resolution = 1e-2)
  group_df = tibble::tibble(cell = row.names(colData(cds)),
                            cell_group = colData(cds)$Cell.type)
  agg_mat = aggregate_gene_expression(cds, module_df, group_df)
  row.names(agg_mat) = stringr::str_c('Module', row.names(agg_mat))
  res = list(module_df, group_df, agg_mat)
  names(res) = c('module_df', 'group_df', 'agg_mat')
  return(res)
}

Myeloid.monocle3res = MonocleModule(Myeloid.rna.cds)
Erythroid.monocle3res = MonocleModule(Erythroid.rna.cds)
Lymphoid.monocle3res = MonocleModule(Lymphoid.rna.cds)

Myeloid_module_heatmap = pheatmap::pheatmap(Myeloid.monocle3res$agg_mat, scale = 'column', color = viridis(10), clustering_method = 'ward.D2', filename = 'Myeloid_module_heatmap.png', width = 5, height = 7)
Erythroid_module_heatmap = pheatmap::pheatmap(Erythroid.monocle3res$agg_mat, scale = 'column', color = viridis(10), clustering_method = 'ward.D2', filename = 'Erythroid_module_heatmap.png', width = 5, height = 7)
Lymphoid_module_heatmap = pheatmap::pheatmap(Lymphoid.monocle3res$agg_mat, scale = 'column', color = viridis(10), clustering_method = 'ward.D2', filename = 'Lymphoid_module_heatmap.png', width = 5, height = 7)

write.csv(Myeloid.monocle3res$module_df, 'Myeloid_module_df.csv')
write.csv(Erythroid.monocle3res$module_df, 'Erythroid_module_df.csv')
write.csv(Lymphoid.monocle3res$module_df, 'Lymphoid_module_df.csv')