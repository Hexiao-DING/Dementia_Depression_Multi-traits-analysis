#@author: Li Na, SCU, SKLB
# -------------------- 1. MAST differential analysis-------------------- 
load(file='Data/seurat_obj.Rdata')
table(seurat_obj$cell.type, seurat_obj$diagnosis)
celltype_of_interest <- ""
seurat_obj_MAST <- subset(seurat_obj, subset = cell.type == celltype_of_interest)
#seurat_obj_MAST$CDR <- colSums(seurat_obj_MAST@assays$RNA@counts > 0) / nrow(seurat_obj_MAST)
de_results <- FindMarkers(
  object = seurat_obj_MAST,
  ident.1 = "DEM",     
  ident.2 = "Control",  
  group.by = "diagnosis",  
  test.use = "MAST",        
  latent.vars = c("sex", "nFeature_RNA"),   
  logfc.threshold = 0.25,   
  min.pct = 0.1,         
  only.pos = TRUE,        
  assay = "RNA" 
)
head(de_results, 10)
library(writexl)
de_results$gene <- rownames(de_results)
write_xlsx(de_results, "MAST_AD/MAST_Endothelial_cell_DEGs.xlsx")
# -------------------- 2. hdWGCNA-------------------- 
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(scater)
Idents(seurat_obj) <- "cell.type"
#1.Create hdWGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.01, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
length(seurat_obj@misc$tutorial$wgcna_genes)
#2.Create metacell
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell.type","sample"),
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cell.type' # set the Idents of the metacell seurat object
)
#3.Dimensionality reduction and batch removal
seurat_obj <- NormalizeMetacells(seurat_obj)
metacell_obj <- GetMetacellObject(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='cell.type')
seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='harmony', dims=1:15)
print(names(seurat_obj@reductions))
DimPlotMetacells(seurat_obj,
                 group.by = 'cell.type',
                 label = TRUE) +
  umap_theme() +
  ggtitle("Cell Type")
DimPlotMetacells(seurat_obj, group.by='sample', label = T) + umap_theme() + ggtitle("Group")
#4.Extract cells for co-expression network analysis
Idents(seurat_obj) <- "cell.type"
DotPlot(seurat_obj, features = "GEMIN7", assay = "RNA")
print(Assays(seurat_obj))
print(unique(seurat_obj$cell.type))
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "Oligodendrocyte",
  group.by = 'cell.type',
  assay = 'SCT',
  slot = 'counts')
#5.Soft threshold selection
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
plot_list <- PlotSoftPowers(seurat_obj)
wrap_plots(plot_list, ncol=2)
#6.Construct an overlapping matrix
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  soft_power=5,
  setDatExpr=FALSE,
  corType = "pearson",
  networkType = "signed",
  TOMType = "signed",
  detectCutHeight = 0.995,
  minModuleSize = 50,
  mergeCutHeight = 0.2,
  tom_outdir = "AD_data/Oli_TOM",
  tom_name = 'MAC' # name of the topoligical overlap matrix written to disk
)
PlotDendrogram(seurat_obj, main='MAC hdWGCNA Dendrogram')
#7.Obtain information about all current modules
modules <- GetModules(seurat_obj)
modules_filtered <- modules %>%
  group_by(module) %>%
  mutate(module_gene_count = n()) %>%  # Add the total number of genes for each line to the corresponding module
  ungroup() %>%
  filter(module != "grey" & module_gene_count >= 20) %>%
  select(-module_gene_count) 
library(openxlsx)
write.xlsx(modules, file = "AD_data/Oli_gene_modules.xlsx", rowNames = FALSE)
#seurat_obj_filtered <- SetModules(seurat_obj, modules = modules_filtered)
#seurat_obj_filtered <- ResetModuleNames(seurat_obj_filtered)
#8.Calculate the ME value
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  scale.model.use="linear",
  assay = "SCT",
  pc_dim = 1)
hMEs <- GetMEs(seurat_obj)
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
#9.differential analysis
me_seurat <- CreateSeuratObject(counts = t(MEs))  # Need to be transposed into module × cell
me_seurat <- SetAssayData(
  me_seurat,
  layer = "data",
  new.data = t(MEs),
  assay = "RNA"
)
me_seurat@meta.data <- seurat_obj@meta.data[colnames(me_seurat), ]
Idents(me_seurat) <- "diagnosis"
de_me <- FindMarkers(
  object = me_seurat,
  ident.1 = "DEM",
  ident.2 = "Control",
  group.by = "diagnosis",  
  test.use = "wilcox", 
  logfc.threshold = 0,
  min.pct = 0,
  assay = "RNA"            
)
de_me
library(writexl)
de_me <- rownames_to_column(de_me, var = "module")
write_xlsx(de_me, "AD_data/hdWGCNA_Oligodendrocyte_neuron_DGE.xlsx")
save(seurat_obj,file="AD_data/Oli_MEseurat.Rdata")
