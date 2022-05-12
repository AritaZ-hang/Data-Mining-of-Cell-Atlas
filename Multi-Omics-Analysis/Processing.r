setwd('YOURPATH')
load('Tissue_unintegrated.RData')

load('Tissue.RData')

library(Seurat)
library(Signac)
library(dplyr)
library(magrittr)
library(Matrix)
library(data.table)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(reticulate)
use_condaenv('scvi', required = T)
library(sceasy)
sc = import('scanpy', convert = F)
scvi = import('scvi', convert = F)
skmisc = import('skmisc', convert = F)
library(MySeuratWrappers)
library(cowplot)

## The merging pipeline is based on the fact that there're 3 independent samples.

##### RNA PIPELINE #####

RestList = list(Tissue2.filtered, Tissue3.filtered)
rna_Tissue = merge(Tissue1.filtered, RestList)
rna_Tissue = NormalizeData(rna_Tissue,
                               normalization.method = 'LogNormalize',
                               scale.factor = 10000)
rna_Tissue = FindVariableFeatures(rna_Tissue,
                                      selection.method = 'vst',
                                      nfeatures = 2000)
rna_Tissue = ScaleData(rna_Tissue, rownames(rna_Tissue))
rna_Tissue = RunPCA(rna_Tissue, features = VariableFeatures(rna_Tissue))

ElbowPlot(rna_Tissue, ndims = 50) # Pick 20

rna_Tissue = FindNeighbors(rna_Tissue, dims = 1:20)
rna_Tissue = FindClusters(rna_Tissue, resolution = 1)
rna_Tissue = RunUMAP(rna_Tissue, dims = 1:20)
rna_Tissue = RunTSNE(rna_Tissue, dims = 1:20)

DimPlot(rna_Tissue, label = T) + NoLegend()

rna_Tissue.markers = FindAllMarkers(rna_Tissue,
                                        only.pos = T,
                                        min.pct = 0.25, logfc.threshold = 0.25)

rna_Tissue.markers.extracted = rna_Tissue.markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

setwd('YOURPATH')
write.csv(rna_Tissue.markers.extracted, quote = F,
          file = 'rna_Tissue.markers.extracted.csv')
write.csv(rna_Tissue.markers, quote = F,
          file = 'rna_Tissue.markers.full.csv')

##### ATAC PIPELINE #####

setwd('YOURPATH')
atac_matrix = readRDS('atac_site_filtered.transfer.rds')

setwd('YOURPATH')
mouse_metadata = fread('YOURPATH/cell_metadata.tissue_freq_filtered.txt')
Tissue_cells = mouse_metadata[mouse_metadata$tissue=='Tissue']$cell
atac_matrix = atac_matrix[,Tissue_cells]
activity_scores = readRDS('activity_scores.quantitative.rds')
activity_scores = activity_scores[,Tissue_cells]
ga_names = rownames(activity_scores)
ga_names = tolower(ga_names)
library(Hmisc)
ga_names = capitalize(ga_names)
rownames(activity_scores) = ga_names

setwd('YOURPATH')
annotations = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) = 'UCSC'
genome(annotations) = 'mm10'
atac_Tissue_assay = CreateChromatinAssay(atac_matrix,
                                             min.cells = 10,
                                             genome = 'mm10',
                                             annotation = annotations,
                                             sep = c('_','_'))
atac_Tissue = CreateSeuratObject(atac_Tissue_assay,
                                     assay = 'peaks',
                                     project = 'ATAC')
atac_Tissue = RunTFIDF(atac_Tissue)


adata_atac = convertFormat(atac_Tissue, from = 'seurat',to = 'anndata',
                           main_layer = 'counts',assay = 'peaks',
                           drop_single_values = F)

## Use peakvi model to project in low dimension.

print(adata_atac)
scvi$model$PEAKVI$setup_anndata(adata_atac)
pvi = scvi$model$PEAKVI(adata_atac)
pvi$train()

scvi$model$PEAKVI$save(pvi, dir_path = 'pvi_trained2')

pvi = scvi$model$PEAKVI$load('pvi_trained2',
                             adata_atac)

latent = pvi$get_latent_representation()
latent = as.matrix(latent)
rownames(latent) = colnames(atac_Tissue)
ndims = ncol(latent)
atac_Tissue[["peakvi"]] = CreateDimReducObject(embeddings = latent, 
                                                   key = "peakvi_", 
                                                   assay = "peaks")

atac_Tissue = FindNeighbors(atac_Tissue, reduction = "peakvi", dims=1:ndims)
atac_Tissue = FindClusters(atac_Tissue, resolution = 1)

atac_Tissue = RunUMAP(atac_Tissue, reduction = "peakvi", dims=1:ndims)

atac_cluster = DimPlot(object = atac_Tissue, 
                       label = TRUE,
                       pt.size = 0.8,
                       label.size = 4,
                       repel = T) + NoLegend() + ggtitle('')+
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = 'NA', color = 'black',
                                    size = 1, linetype = 'solid')) +
  labs(x = 'UMAP1', y = 'UMAP2')
ggsave(atac_cluster, filename = 'atac_cluster.png',
       width = 6, height = 6)

atac_Tissue[['RNA']] = CreateAssayObject(counts = activity_scores)
atac_Tissue = NormalizeData(object = atac_Tissue,
                                assay = 'RNA',
                                normalization.method = 'LogNormalize',
                                scale.factor = median(atac_Tissue$nCount_RNA)
)
DefaultAssay(atac_Tissue) = 'RNA'

adata_rna = convertFormat(rna_Tissue, from="seurat", to="anndata", main_layer="counts", assay="RNA", drop_single_values=FALSE,
                          outFile = 'adata_rna.h5ad')
adata_atac_act = convertFormat(atac_Tissue, from="seurat", to="anndata", main_layer="counts", assay="RNA", drop_single_values=FALSE,outFile = 'adata_atac_act.h5ad')

##### SEE THE scANVI PIPELINE #####

df = read.table(file = 'scvi_annotation.csv', header = T,
                sep = ',') # The result file of scANVI pipeline

atac_Tissue = AddMetaData(object = atac_Tissue, metadata = df$predicted.labels, col.name="predicted.labels")

options(repr.plot.width=12, repr.plot.height=5)

atac_umap_embeddings = atac_Tissue@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(Cell.type = atac_Tissue$predicted.labels)

ATAC_Cell_type_annotation = DimPlot(object = atac_Tissue, 
                                    group.by = 'predicted.labels',
                                    label = TRUE,
                                    pt.size = 0.7,
                                    label.size = 4,
                                    repel = T) + NoLegend() + ggtitle('')+
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = 'NA', color = 'black',
                                    size = 1, linetype = 'solid')) +
  labs(x = 'UMAP1', y = 'UMAP2')
ggsave(ATAC_Cell_type_annotation, filename = 'ATAC_Cell_type_annotation.png', width = 7, height = 7,
       dpi = 600)

# draw marker plot
to_draw.markers = c('YOURMARKERS')
markers_vlnplot = VlnPlot(rna_Tissue, 
                          features = to_draw.markers, 
                          stack = T, 
                          pt.size = 0, 
                          direction = 'horizontal', 
                          x.lab = '', y.lab = '', 
                          group.by = 'Cell.type',
                          features.face = 'italic') +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6,
                                   face = 'italic'))
ggsave(markers_vlnplot, filename = 'markers_vlnplot.png',
       width = 6, height = 4)

rna_dotplot = DotPlot(rna_Tissue,
                      features = to_draw.markers,group.by = 'Cell.type', scale = T)+ coord_flip() + theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size = 12, face = 'bold'),
        axis.text.y = element_text(size = 12, face = 'bold'),
        legend.title = element_text(size = 16),
        legend.key.size = unit(5, 'mm'),
        legend.text = element_text(size = 12, face = 'bold')) + 
  scale_color_gradientn(values = seq(0,1,0.2), colors = c('#330066','#336699','#66CC66','#FFCC33'))+
  labs(x = NULL, y = NULL)+guides(size = guide_legend(order = 2))
ggsave(rna_dotplot, filename = 'Tissue_rna_dotplot.png', width = 6, height = 8)

levels(rna_Tissue$seurat_clusters) = 1:(max(as.integer((levels(rna_Tissue$seurat_clusters)))) + 1)
rna_cluster = Seurat::DimPlot(object = rna_Tissue, 
                              label = TRUE,
                              pt.size = 0.8,
                              label.size = 4,
                              group.by = 'seurat_clusters',
                              repel = T) + NoLegend() + ggtitle('')+
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = 'NA', color = 'black',
                                    size = 1, linetype = 'solid')) +
  labs(x = 'UMAP1', y = 'UMAP2') + theme_classic()
ggsave(rna_cluster, filename = 'Tissue_rna_cluster.png',
       width = 6, height = 6)

save.image(file = 'Tissue.RData')