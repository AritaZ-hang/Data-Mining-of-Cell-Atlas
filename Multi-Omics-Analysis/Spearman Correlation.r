# compute spearman correlation coefficient for each RNA annotation cell type and ATAC annotation cell type
# Then correct the ATAC annotation based on scc

# Extract RNA expr matrix from Seurat object

rna_expM = GetAssayData(rna_Tissue, slot = 'data')
Highly_Variable = VariableFeatures(rna_Tissue)
rna_expM = as.data.frame(rna_expM[Highly_Variable, ])

types = unique(Cell_type$Cell.type)

# compute average expression for 2000 HVGs in a certain cell type
RnaExpr_pertype = function(rna_expM, Cell_type, types)
{
  RnaExpr_pertype = c()
  for( i in types)
  {
    selected.cells = Cell_type[Cell_type$Cell.type == i,]$Cell.name
    rna_expM.extracted = rna_expM[,selected.cells]
    avg_expr = apply(rna_expM.extracted, 1, mean)
    RnaExpr_pertype = cbind(RnaExpr_pertype, avg_expr)
  }
  RnaExpr_pertype = as.data.frame(RnaExpr_pertype)
  colnames(RnaExpr_pertype) = types
  return(RnaExpr_pertype)
}
Rna_type_expr = RnaExpr_pertype(rna_expM, Cell_type, types)

# Extract ATAC gene activity matrix from Seurat object
DefaultAssay(atac_Tissue) = 'RNA'
Ga_expM = GetAssayData(atac_Tissue, slot = 'data')
rownames(Ga_expM) = toupper(rownames(Ga_expM))
Highly_Variable_extract = toupper(Highly_Variable)[toupper(Highly_Variable) %in% rownames(Ga_expM)] # change the format of GA's genes
Ga_expM = Ga_expM[Highly_Variable_extract,]


atac_types = compare_ATAC_annotations[,c('cell','pred.label')]
names(atac_types) = c('Cell.name', 'Cell.type')

# compute average gene activity for HVGs in a certain cell type
GaExpr_pertype = function(Ga_expM, atac_types, types)
{
  GaExpr_pertype = c()
  for(i in types)
  {
    selected.cells = atac_types[atac_types$Cell.type == i,]$Cell.name
    Ga_expM.extracted = Ga_expM[,selected.cells]
    avg_expr = apply(Ga_expM.extracted, 1, mean)
    GaExpr_pertype = cbind(GaExpr_pertype, avg_expr)
  }
  GaExpr_pertype = as.data.frame(GaExpr_pertype)
  colnames(GaExpr_pertype) = types
  return(GaExpr_pertype)
}

Ga_type_expr = GaExpr_pertype(Ga_expM, atac_types, types)

rownames(Rna_type_expr) = toupper(rownames(Rna_type_expr))
Rna_type_expr = Rna_type_expr[Highly_Variable_extract,]

# compute the spearman correlation coefficient for each RNA annotation cell type & ATAC annotation cell type

GetCorrelationMatrix = function(Rna_type_expr, Ga_type_expr, types)
{
  types_num = length(types)
  CorrelationMatrix = matrix(nrow = types_num, ncol = types_num)
  for(i in 1:types_num)
  {
    for(j in 1:types_num)
    {
      CorrelationMatrix[i,j] = cor(Rna_type_expr[,i], Ga_type_expr[,j], method = 'spearman')
    }
  }
  CorrelationMatrix = as.data.frame(CorrelationMatrix)
  names(CorrelationMatrix) = types
  rownames(CorrelationMatrix) = types
  return(CorrelationMatrix)
}
correlation_matrix = GetCorrelationMatrix(Rna_type_expr, Ga_type_expr, types)

# draw a heatmap for visualization

library(pheatmap)
library(viridis)
pheatmap::pheatmap(correlation_matrix,
                   color = viridis(10),
                   cluster_rows = F,
                   cluster_cols = F,
                   angle_col = 315,
                   fontsize = 6,
                   filename = 'Pheatmap.png',
                   width = 3,
                   height = 3)

write.csv(correlation_matrix, 'correlation_matrix.csv')

# set median spearman correlation as threshold
setThreshold = function(x)
{
  y = c()
  types.num = length(rownames(x))
  for(i in 1:types.num)
    y = c(y, x[,i])
  threshold = median(y, na.rm = T)
  return(threshold)
}

#atac_object: atac Seurat object
#Rna_type_expr: gene average expression per type
#Ga_expM: filtered highly variable genes' activity
#correlation_matrix: spearman correlation matrix
#Return a annotation-refined metadata(data.frame)

RefineAnnotations = function(atac_object, Rna_type_expr, Ga_expM, 
                             correlation_matrix)
{
  #create a data frame
  RefinedAnnotation = data.frame(cell = rownames(atac_object@meta.data),
                                 refined.anno = atac_object$predicted.labels)
  types = unique(atac_object$predicted.labels)
  rna_types = colnames(Rna_type_expr)
  threshold = setThreshold(correlation_matrix) # median value of the spearman correlation coefficients
  
  for(i in types)
  {
    if(correlation_matrix[i,i] < max(correlation_matrix[,i]) & !is.na(correlation_matrix[i, i]))
    {
      old.type = i
      selected.cells = RefinedAnnotation[RefinedAnnotation$refined.anno == i,]$cell
      Ga_expM.extracted = Ga_expM[, selected.cells]
      avg_expr = apply(Ga_expM.extracted, 1, mean)
      new.spearman = c()
      for(j in rna_types)
      {
        spearman = cor(avg_expr, Rna_type_expr[, j], method = 'spearman')
        new.spearman = c(new.spearman, spearman)
      }
      names(new.spearman) = rna_types
      new.spearman = new.spearman[order(new.spearman, decreasing = T)]
      if(new.spearman[1] >= threshold)
      {
        new.type = names(new.spearman)[1]
        RefinedAnnotation[RefinedAnnotation$cell %in% selected.cells,]$refined.anno = gsub(old.type,new.type,RefinedAnnotation[RefinedAnnotation$cell %in% selected.cells,]$refined.anno)
      }
      else
        RefinedAnnotation[RefinedAnnotation$cell %in% selected.cells,]$refined.anno = gsub(old.type,'Undefined',RefinedAnnotation[RefinedAnnotation$cell %in% selected.cells,]$refined.anno)
    }
  }
  return(RefinedAnnotation)
}

refined.anno = RefineAnnotations(atac_Tissue, Rna_type_expr, Ga_expM, correlation_matrix)

atac_Tissue = AddMetaData(atac_Tissue, refined.anno$refined.anno, 
                         col.name = 'refined.anno')

## compute new spearman correlation coefficient
refined.types = unique(atac_Tissue$refined.anno)
atac_types.refined = refined.anno
names(atac_types.refined) = c('Cell.name', 'Cell.type')
Ga_type_expr.refined = GaExpr_pertype(Ga_expM, atac_types.refined, refined.types)
correlation_matrix.refined = Essay_LT_CorrelationMatrix(Rna_type_expr, Ga_type_expr.refined, types, refined.types)

write.csv(correlation_matrix.refined, 'correlation_matrix_refined.csv')
write.csv(refined.anno, 'Tissue_refined_anno.csv')