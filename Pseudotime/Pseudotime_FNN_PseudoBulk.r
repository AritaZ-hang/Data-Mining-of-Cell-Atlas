options(stringsAsFactors = F)
setwd('D:/Main')
rm(list = ls())
set.seed(0)
library(Seurat)
library(Signac)
library(FNN)
library(Hmisc)
library(monocle)
library(dplyr)
library(chromVAR)
library(viridis)
library(data.table)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(stringr)

# some functions
# A simple function to save plots
SavePlots = function(gg, tosaveName, width, height)
{
  filename = tosaveName
  ggsave(gg, filename = filename, width = width, height = height)
}

# loading necessary files
load('Pseudotime_RNA_files.RData')
load('Pseudotime_ATAC_files.RData')
load('Pseudotime_RNA_meta.RData')

gene.activity = readRDS(file = 'D:/activity_scores.quantitative.rds')

# a simple function to add corresponding gene activity scores to query atac object, and perform basic analysis according to the tutorials of seurat
GeneActivityPipeline = function(object, gene.activity)
{
  cells = colnames(object)
  gene.activity = gene.activity[, cells]
  rownames(gene.activity) = capitalize(tolower(rownames(gene.activity)))
  object[['activity']] = CreateAssayObject(counts = gene.activity)
  DefaultAssay(object) = 'activity'
  object = NormalizeData(object)
  object = ScaleData(object, features = rownames(object))
  return(object)
}

Erythroid.atac = GeneActivityPipeline(Erythroid.atac, gene.activity)

# Create Transfer anchors in CCA L2 Space
transfer.anchors.Erythroid = FindTransferAnchors(reference = Erythroid.rna, query = Erythroid.atac, reference.assay = 'RNA', query.assay = 'activity',reduction = 'cca',k.anchor = 10, features = VariableFeatures(Erythroid.rna))

# return with a list consists of refdata matrix & imputation data matrix
ImputeCoembed = function(rna_object, atac_object, transfer.anchors, type)
{
  genes.use = VariableFeatures(rna_object)
  refdata = GetAssayData(rna_object, assay = 'RNA', slot = 'data')[genes.use,]
  imputation = TransferData(transfer.anchors, refdata, weight.reduction = atac_object[['lsi']], dims = 2:30)
  atac_object[['RNA']] = imputation
  DefaultAssay(atac_object) = 'RNA'
  rna_object$orig.ident = rep('RNA', dim(rna_object@meta.data)[1])
  coembed = merge(x = rna_object, y = atac_object)
  coembed = ScaleData(coembed, features = genes.use, do.scale = F)
  coembed = RunPCA(coembed, features = genes.use, verbose = F)
  coembed = RunUMAP(coembed, dims = 1:30)
  p1 = DimPlot(coembed, group.by = c('orig.ident')) + ggtitle('') + theme_void() + theme(legend.key.size = unit(1, 'mm'), legend.text = element_text(size = 4, face = 'bold'))
  tosaveName = paste('Psuedotime_coembedding_', type, '_plot.png', sep = '')
  ggsave(p1, filename = tosaveName, width = 3, height = 3)
  
  refdata.matrix = t(as.matrix(refdata))
  imputation.matrix = GetAssayData(imputation, slot = 'data')
  imputation.matrix = t(as.matrix(imputation.matrix))
  
  return(list(refdata.matrix, imputation.matrix))
  
  
}

Erythroid.impute.ref = ImputeCoembed(Erythroid.rna, Erythroid.atac, transfer.anchors.Erythroid, type = 'Erythroid')

#FNN & Index Matching

FnnMatching = function(impute.ref, metadata.rna, metadata.atac)
{
  meta.rna = data.frame(cell = rownames(metadata.rna),
                        Cell.type = metadata.rna$Cell.type)
  meta.atac = data.frame(cell = rownames(metadata.atac),
                         Cell.type = metadata.atac$Cell.type)
  rna.types = unique(meta.rna$Cell.type)
  atac.types = unique(meta.atac$Cell.type)
  
  cover = atac.types[atac.types %in% rna.types]
  rest = setdiff(atac.types, rna.types) # well I don't know where it belongs
  
  knn.index = data.frame()
  for(i in cover)
  {
    rna.cells = meta.rna[meta.rna$Cell.type == i,]$cell
    atac.cells = meta.atac[meta.atac$Cell.type == i,]$cell
    knn.grouping = get.knnx(data = impute.ref[[2]][atac.cells,],
                            query = impute.ref[[1]][rna.cells,],
                            algorithm = 'kd_tree', k = 1)
    temp.index = as.data.frame(knn.grouping[['nn.index']])
    rownames(temp.index) = rna.cells
    rna_index.df = data.frame(rna_cells = rna.cells,
                              index = temp.index[,1])
    atac_index.df = data.frame(atac_cells = atac.cells,
                               index = 1:length(atac.cells))
    final.index = inner_join(atac_index.df, rna_index.df, by = 'index')
    knn.index = rbind(knn.index, final.index)
  }
  if(length(rest)) # if the types of atac is different from those in RNA
  {
    for(i in rest)
    {
      atac.cells = meta.atac[meta.atac$Cell.type == i,]$cell
      knn.grouping = get.knnx(data = impute.ref[[2]],
                              query = impute.ref[[1]][rna.cells,],
                              algorithm = 'kd_tree', k = 1)
      temp.index = as.data.frame(knn.grouping[['nn.index']])
      rownames(temp.index) = rna.cells
      rna.cells = rownames(impute.ref[[1]])
      rna_index.df = data.frame(rna_cells = rna.cells,
                                index = temp.index[,1])
      atac_index.df = data.frame(atac_cells = atac.cells,
                                 index = 1:length(atac.cells))
      final.index = inner_join(atac_index.df, rna_index.df, by = 'index')
      knn.index = rbind(knn.index, final.index)
    }
  }
  return(knn.index)
}
Erythroid.matching.df = FnnMatching(Erythroid.impute.ref, Erythroid.rna@meta.data, Erythroid.atac@meta.data)

# Using Kmeans to create pseudobulk for RNA expr. As for ATAC peak count matrix, we use the result of pseudobulk to create corresponding ATAC pseudobulk.
# return with RNA pseudobulk(kmeans) & grouping(pseudogroups).

MakePseudo_RNA = function(rna.object, cell.types, k)
{
  cells = rownames(rna.object@meta.data[rna.object$Cell.type %in% cell.types,])
  rna.expr = GetAssayData(rna.object, assay = 'RNA', slot = 'counts')
  rna.expr = LogNormalize(rna.expr)
  rna.expr = as.data.frame(t(as.matrix(rna.expr)))
  
  fit = kmeans(rna.expr, iter.max = 100, centers = k)
  expr.kmeans = aggregate(rna.expr, by = list(fit$cluster), FUN = mean)
  expr.kmeans = subset(expr.kmeans, select = -c(Group.1))
  
  expr.kmeans = apply(expr.kmeans, MARGIN = 2, function(x)(x-mean(x))/sd(x))
  expr.kmeans = as.data.frame(expr.kmeans)
  
  cluster.df = data.frame(rna_cells = rownames(rna.expr),
                          pseudocell = fit$cluster)
  
  res = list(expr.kmeans, cluster.df)
  names(res) = c('kmeans', 'pseudogroups')
  return(res)
}
Erythroid.rna.pseudo = MakePseudo_RNA(Erythroid.rna, cell.types = c('Erythroid cell', 'Hematopoietic Stem Progenitor cell'), k = 50)

Erythroid.rna.kmeans = Erythroid.rna.pseudo$kmeans
Erythroid.rna.pseudogroups = Erythroid.rna.pseudo$pseudogroups

# Make pseudobulk for ATAC peak count matrix. For each pseudobulk sample, the counts of a certain peak is the sum of counts of all cells consisting in.
# return with the ATAC pseudobulk with TF-IDF normalized.

MakePseudo_ATAC = function(atac.object, matching.df, pseudogroups)
{
  atac.counts = GetAssayData(atac.object, assay = 'peaks', slot = 'counts')
  atac.counts = as.matrix(atac.counts)
  rna.groups = 1:50
  
  Pseudo.atac = c()
  atac.groups = 0
  for(i in rna.groups)
  {
    cells.in = pseudogroups[pseudogroups$pseudocell == i,]$rna_cells
    corresponding.atac = matching.df[matching.df$rna_cells %in% cells.in,]$atac_cells
    if(length(corresponding.atac) > 1)
    {
      pseudo.atac = atac.counts[, corresponding.atac]
      new.counts = rowSums(pseudo.atac)
      Pseudo.atac = cbind(Pseudo.atac, new.counts)
      atac.groups = atac.groups + 1
    }
    else if(length(corresponding.atac) == 1)
    {
      pseudo.atac = atac.counts[,corresponding.atac]
      new.counts = pseudo.atac
      Pseudo.atac = cbind(Pseudo.atac, new.counts)
      atac.groups = atac.groups + 1
    }
  }
  Pseudo.atac = as.data.frame(Pseudo.atac)
  rownames(Pseudo.atac) = rownames(atac.counts)
  colnames(Pseudo.atac) = paste("Pseudo_", 1:atac.groups, sep = '')
  
  Pseudo.atac = RunTFIDF(Pseudo.atac)
  
  return(Pseudo.atac)
}
Pseudo.atac.counts = MakePseudo_ATAC(Erythroid.atac, Erythroid.matching.df, Erythroid.rna.pseudogroups)

# compute average pseudotime for RNA/ATAC pseudobulk samples.

MeanPseudotime = function(matching.df, pseudogroups, metadata)
{
  metadata$rna_cells = rownames(metadata)
  matching = inner_join(matching.df, metadata, by = 'rna_cells')[, c('atac_cells','rna_cells', 'Pseudotime')]
  Pseudo.rna = inner_join(pseudogroups, matching, by = 'rna_cells')
  
  groups = 1:50
  rna.mean.pseudotime = c()
  for( i in groups)
  {
    pseudotime = mean(Pseudo.rna[Pseudo.rna$pseudocell == i,]$Pseudotime)
    rna.mean.pseudotime = c(rna.mean.pseudotime, pseudotime)
  }
  names(rna.mean.pseudotime) = paste('Pseudo_', groups, sep = '')
  
  return(rna.mean.pseudotime)
}
mean.pseudotime = MeanPseudotime(Erythroid.matching.df, Erythroid.rna.pseudogroups, Erythroid.metadata)

pfm = getMatrixSet(x = JASPAR2020, opts = list(collection = 'CORE', tax_group = 'vertebrates', all_versions = FALSE))

SplitPeaks = function(peaks)
{
  full= unlist(str_split(peaks, pattern = '-'))
  peaks.matrix = matrix(ncol = 3, nrow = length(full)/3)
  pairs = length(full)/3
  for(i in 1:pairs){
    peaks.matrix[i, 1] = full[i*3-2]
    peaks.matrix[i, 2] = full[i*3-1]
    peaks.matrix[i, 3] = full[i*3]
  }
  peaks.matrix = as.data.frame(peaks.matrix)
  names(peaks.matrix) = c('chr', 'start', 'end')
  peaks.matrix$start = as.numeric(peaks.matrix$start)
  peaks.matrix$end = as.numeric(peaks.matrix$end)
  return(peaks.matrix)
}
# Create TF gene - TF mapping relationship.
MotifMatchingDf = function(pfm)
{
  len = length(pfm)
  motif.id = c()
  motif.name = c()
  for(i in 1:len)
  {
    id = pfm@listData[[i]]@ID
    name = pfm@listData[[i]]@name
    motif.id = c(motif.id, id)
    motif.name = c(motif.name, name)
  }
  motif.matching.df = data.frame(id = motif.id,
                                 name = motif.name)
  return(motif.matching.df)
}
MotifSymbolSynchro = function(Motif.matching)
{
  motif.names = Motif.matching$name
  temp = motif.names
  temp = gsub('\\(var\\.[0-9]\\)', '', temp)
  temp = unique(temp)
  len = length(temp)
  MotifSymbolSynchro = data.frame()
  MotifSymbolSynchro = rbind(MotifSymbolSynchro, c('', ''))
  for(i in 1:len)
  {
    if(stringr::str_detect(temp[i], '::'))
    {
      str.list = stringr::str_split(temp[i], '::')
      strs = unlist(str.list)
      row1 = c(temp[i], strs[1])
      if((strs[1] %in% MotifSymbolSynchro[,2]) | (capitalize(tolower(strs[1])) %in% MotifSymbolSynchro[,2]) ){}
      else
        MotifSymbolSynchro = rbind(MotifSymbolSynchro, row1)
    }
    else
    {
      str = temp[i]
      row1 = c(temp[i], str)
      if((str %in% MotifSymbolSynchro[,2]) | (capitalize(tolower(str)) %in% MotifSymbolSynchro[,2])){}
      else
        MotifSymbolSynchro = rbind(MotifSymbolSynchro, row1)
    }
  }
  names(MotifSymbolSynchro) = c('motif.name', 'symbol')
  MotifSymbolSynchro$symbol = capitalize(tolower(MotifSymbolSynchro$symbol))
  MotifSymbolSynchro = unique(MotifSymbolSynchro)
  MotifSymbolSynchro = MotifSymbolSynchro[2:dim(MotifSymbolSynchro)[1],]
  return(MotifSymbolSynchro)
}
MotifQuery = function(query, motif2symbol, rna.expr)
{
  query = query[query %in% colnames(rna.expr)]
  query = query[query %in% motif2symbol$symbol]
  
  ref = motif2symbol[motif2symbol$symbol %in% query,]$motif.name
  
  MotifQuery = list(ref, query)
  names(MotifQuery) = c('ref', 'query')
  return(MotifQuery)
}
RenameDevMatrix = function(dev, motif.matching)
{
  col.names = data.frame(id = colnames(dev))
  temp = inner_join(col.names, motif.matching, by = 'id')
  new.name = temp$name
  colnames(dev) = new.name
  return(dev)
}
# compute TF activity for ATAC pseudobulks. see the vignette of chromVAR.
RecalculateChromVAR = function(atac.counts, pfm)
{
  peaks = rownames(atac.counts)
  peaks.matrix = SplitPeaks(peaks)
  peaks.grange = GRanges(seqnames = peaks.matrix$chr, 
                         IRanges(start = peaks.matrix$start,
                                 end = peaks.matrix$end))
  counts = as.matrix(atac.counts)
  
  fragment_counts = SummarizedExperiment(assays = list(counts = counts),
                                         rowRanges = peaks.grange)
  
  fragment_counts = addGCBias(fragment_counts, genome = BSgenome.Mmusculus.UCSC.mm10)
  fragment_counts = filterPeaks(fragment_counts, non_overlapping = T)
  motif.ix = matchMotifs(pfm, fragment_counts, genome = BSgenome.Mmusculus.UCSC.mm10)
  dev = computeDeviations(fragment_counts, annotations = motif.ix)
  dev.matrix = deviations(dev)
  dev.matrix = t(dev.matrix)
  dev.matrix = as.data.frame(dev.matrix)
  
  motif.matching = MotifMatchingDf(pfm)
  dev.matrix = RenameDevMatrix(dev.matrix, motif.matching)
  dev.matrix = apply(dev.matrix, MARGIN = 2, function(x)(x-mean(x))/sd(x))
  dev.matrix = as.data.frame(dev.matrix)
  
  return(dev.matrix)
}

Pseudo.atac.chromvar = RecalculateChromVAR(Pseudo.atac.counts, pfm)

# Add Pseudotime to pseudobulks, and sort them based on the ascending pseudotime.
AddPseudotime = function(object, mean.pseudotime)
{
  object = as.data.frame(object)
  object$Pseudotime = as.numeric(mean.pseudotime)
  object = object[order(object$Pseudotime),]
  return(object)
}
Pseudo.rna.final = AddPseudotime(Erythroid.rna.kmeans, mean.pseudotime)
Pseudo.atac.final = AddPseudotime(Pseudo.atac.chromvar, mean.pseudotime)

# Using Loess Regression to smooth the data
ExtractMatrix = function(kmeans, genes)
{
  extract = kmeans[, genes]
  extract = t(extract)
  
  return(extract)
}
LoessRegression = function(kmeans, genes, length)
{
  extract = t(ExtractMatrix(kmeans, genes))
  Pseudotime = kmeans$Pseudotime
  exp = c()
  for(i in genes)
  {
    temp.df = data.frame(gene.exp = extract[, i],
                         Pseudotime = Pseudotime)
    Loess = loess(temp.df$gene.exp ~ temp.df$Pseudotime, temp.df)
    res = predict(Loess, data.frame(Pseudotime = seq(0, 15, length = length)), se = T)
    loess.exp = res$fit
    exp = cbind(exp, loess.exp)
  }
  exp = as.data.frame(exp)
  names(exp) = genes
  return(exp)
}

Motif.matching = MotifMatchingDf(pfm)
motif2symbol = MotifSymbolSynchro(Motif.matching)

full_query = motif2symbol$symbol
full_query.res = MotifQuery(full_query, motif2symbol, Erythroid.rna.kmeans)
full_query.res$ref = sort(unique(full_query.res$ref))
full_query.res$query = sort(unique(full_query.res$query))


Pseudo.rna.full = ExtractMatrix(Pseudo.rna.final, genes = full_query.res$query)
Pseudo.atac.full = ExtractMatrix(Pseudo.atac.final, genes = full_query.res$ref)

Pseudo.rna.full = na.omit(Pseudo.rna.full)
Pseudo.atac.full = na.omit(Pseudo.atac.full)

final.genes = rownames(Pseudo.rna.full)

Pseudo.rna.loess = LoessRegression(Pseudo.rna.final, genes = final.genes, length = 50)
Pseudo.atac.loess = LoessRegression(Pseudo.atac.final, genes = full_query.res$ref, length = 50)

p1 = pheatmap::pheatmap(t(Pseudo.rna.loess), cluster_rows = T, cluster_cols = F, show_rownames = T, show_colnames = F, border_color = NA, scale = 'row', breaks = unique(c(seq(-2, 2, length = 100))), width = 10, height = 35, filename = 're_RNA.png')

pheatmap::pheatmap(t(Pseudo.atac.loess), cluster_rows = T, cluster_cols = F, show_rownames = T, show_colnames = F, border_color = NA, scale = 'row', breaks = unique(c(seq(-2, 2, length = 100))), width = 10, height = 35, filename = 're_ATAC.png')

# Extract some key TFs to draw a heatmap.
rna.extract = ExtractMatrix(Pseudo.rna.final, c('Spi1','Spib','Runx1','Tcf3','Bach2','Maf','Gata2','Zbtb7a','E2f4','Tal1','Gata1', 'Klf1'))
atac.extract = ExtractMatrix(Pseudo.atac.final, c('SPI1','SPIB','RUNX1','TCF3','BACH2','MAF::NFE2','GATA2','ZBTB7A','E2F4','GATA1::TAL1','GATA1::TAL1','Klf1'))

pheatmap::pheatmap((rna.extract), cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F, border_color = NA, scale = 'row', breaks = unique(c(seq(-2,2, length = 100))), width = 3, height = 6, filename = 'Kmeans_rna_Erythroid_new.png')

pheatmap::pheatmap((atac.extract), cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F, border_color = NA, scale = 'row', breaks = unique(c(seq(-2,2, length = 100))), width = 3, height = 6 ,filename = 'Kmeans_atac_Erythroid_new.png')

# compute the pearson correlation coefficient of TF genes expression and TF activity.
Correlation_matrix = function(extract, extract.atac)
{
  Correlation_matrix = matrix(ncol = length(rownames(extract)), nrow = length(rownames(extract.atac)))
  for(i in 1:length(rownames(extract)))
  {
    for(j in 1:length(rownames(extract.atac)))
    {
      Correlation_matrix[i, j] = cor(extract[i, ], extract.atac[j, ], method = 'spearman',use = 'pairwise.complete.obs')
    }
  }
  Correlation_matrix = as.data.frame(Correlation_matrix)
  names(Correlation_matrix) = rownames(extract.atac)
  rownames(Correlation_matrix) = rownames(extract)
  return(Correlation_matrix)
}

cor.matrix = Correlation_matrix(rna.extract, atac.extract)

pheatmap::pheatmap(cor.matrix, cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = T, border_color = NA, scale = 'none', breaks = unique(c(seq(-0.5, 0.5, length = 100))), width = 6, height = 6, filename = 'Kmeans_cor_Erythroid_new.png')
