#@author: Li Na, SCU, SKLB
library(viridis)
library(Seurat)
library(scCustomize)
library(MAST)
setwd("")
# -------------------- DEM data --------------------
# -------------------- --------------------------------
{
  seurat_obj <- readRDS("Data/GSE303823_rna_annotated.rds")
  table(seurat_obj@meta.data[["group"]])
  seurat_obj$diagnosis <- dplyr::case_when(
    seurat_obj$group == "CTRL" ~ "Control",
    seurat_obj$group %in% c("ADD", "DLB", "PDD") ~ "DEM",
    TRUE ~ NA_character_
  )
  table(seurat_obj$group, useNA = "ifany")
  table(seurat_obj$diagnosis, useNA = "ifany")
  save(seurat_obj,file='Data/seurat_obj.Rdata')
  seurat_ad <- subset(seurat_obj, subset = group %in% c("CTRL", "ADD"))
  save(seurat_ad,file='Data/seurat_ad.Rdata')
}
# -------------------- MDD data --------------------
# -------------------- --------------------------------
# -------------------- 1. Reading data：matrix_data --------------------
folder_path <- "MDD/GSE213982"
matrix_data <- ReadMtx(
  mtx = file.path(folder_path, "matrix.mtx.gz"),
  cells = file.path(folder_path, "barcodes.tsv.gz"),
  features = file.path(folder_path, "features.tsv.gz"),
  feature.column = 1
)
save(matrix_data,file="MDD_matrix_data.Rdata")
# -------------------- 2. Parse the barcode information：cell_type--------------------
barcodes <- read.delim(file.path(folder_path, "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
full_strings <- barcodes[[1]]
split_list <- strsplit(full_strings, "\\.")
sample_id <- sapply(split_list, function(x) x[1])
table(sample_id)
cell_type_info <- sapply(split_list, function(x) paste(tail(x, 2), collapse = ".")) 
save(sample_id,cell_type_info,file="MDD_cell_type.Rdata")
# -------------------- 3. Load sample information file --------------------
library(GEOquery)
gse1 <- getGEO("GSE144136", GSEMatrix = TRUE)
gse2 <- getGEO("GSE213982", GSEMatrix = TRUE)
sample_info0 <- pData(phenoData(gse1[[1]]))
sample_info0$title <- sub("^(\\d+).*", "M\\1", sample_info0$title)
sample_info1 <- pData(phenoData(gse2[[1]]))
sample_info2 <- pData(phenoData(gse2[[2]]))
sample_info <- rbind(sample_info0[,c("title","geo_accession","group:ch1","Sex:ch1")],
                     sample_info1[,c("title","geo_accession","group:ch1","Sex:ch1")],
                     sample_info2[,c("title","geo_accession","group:ch1","Sex:ch1")])
colnames(sample_info)
table(sample_info$`group:ch1`)
table(sample_info$`Sex:ch1`)
sample_sheet <- data.frame(
  sample_id = sample_info$title,
  diagnosis = ifelse(sample_info$`group:ch1` == 'Control', 'Control', 'MDD'),
  sex = sample_info$`Sex:ch1`
)
sample_sheet <- rbind(sample_sheet, cbind(sample_id = "M24_2", sample_sheet[sample_sheet$sample_id == "M24", -1]))
save(sample_sheet,file="MDD_sample_sheet.Rdata")
table(sample_sheet$diagnosis)
# -------------------- 4. Create a Seurat object --------------------
load(file="MDD_matrix_data.Rdata")
load(file="MDD_cell_type.Rdata")
load(file="MDD_sample_sheet.Rdata")
meta_df <- data.frame(
  row.names = full_strings, 
  sample = sample_id,
  cell_type = cell_type_info,
  diagnosis = sample_sheet$diagnosis[match(sample_id, sample_sheet$sample_id)],
  sex = sample_sheet$sex[match(sample_id, sample_sheet$sample_id)],
  stringsAsFactors = FALSE
)
seurat_obj <- CreateSeuratObject(
  counts = matrix_data,
  project = "SeuratObject",
  min.cells = 3,    
  min.features = 200,  
)
seurat_obj <- AddMetaData(seurat_obj, metadata = meta_df)
# -------------------- 5. quality control --------------------
load(file="Data/预处理前seurat.Rdata")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP")
mask1 <- seurat_obj$nCount_RNA >= 0 & seurat_obj$nCount_RNA<= 50000
mask2 <- seurat_obj$nFeature_RNA >=0 & seurat_obj$nFeature_RNA <= 4000
mask3 <- seurat_obj$percent.mt <= 15
mask4<-seurat_obj$percent.rb<= 40
mask <- mask1 & mask2 & mask3 & mask4
seurat_obj <- seurat_obj[, mask]
# -------------------- 6. Normalization and Dimensionality Reduction --------------------
##LogNormalize
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of high-variable genes
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
#normalization
#all.genes <- rownames(seurat_obj)
#source('code/genes_to_scale.R')
seurat_obj <- ScaleData(seurat_obj)
#PCA
seurat_obj <- Seurat::RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))#features = VariableFeatures(object = seurat_obj),npcs = 100
# -------------------- 7. Batch processing--------------------
#harmony 
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "sample")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
# -------------------- 8. Clustering dimension reduction--------------------
#source('code/genes_to_scale.R')
seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)
JackStrawPlot(seurat_obj, dims = 1:20)
ElbowPlot(seurat_obj,ndims = 30)
seuratPC=20
seurat_obj=FindNeighbors(seurat_obj, k.param = 30, dims = 1:seuratPC, reduction = "harmony")
#?FindNeighbors
#for (res in c(0.01,0.05,0.1,1,1.5,2,2.5,3,3.5,4)) {seurat_obj=FindClusters(seurat_obj, graph.name = "RNA_snn", resolution = res, algorithm = 1)}
#apply(seurat_obj@meta.data[,grep("RNA_snn_res",colnames(seurat_obj@meta.data))],2,table)
#p2_tree=clustree(seurat_obj@meta.data, prefix = "RNA_snn_res.")
#p2_tree
# Select the resolution for dimensionality reduction
px=0.7
seurat_obj <- FindClusters(seurat_obj, resolution = px)
## MAP/tSNE
seurat_obj <- RunUMAP(seurat_obj, dims = 1:seuratPC, reduction = "harmony")
table(seurat_obj@meta.data$cell_type)
# -------------------- 9. Add cell annotation information --------------------
split_list <- strsplit(as.character(seurat_obj@meta.data$cell_type), "\\.")
raw_cell_type <- sapply(split_list, function(x) x[1])
sub_cell_type <- sapply(split_list, function(x) x[2])
map_cell_type <- function(raw_type) {
  cell_type_lower <- tolower(as.character(raw_type))
  if (grepl("^exn$", cell_type_lower)) return("Excitatory_neurons")
  if (grepl("^inn$", cell_type_lower)) return("Inhibitory_neurons")
  if (grepl("^oli$", cell_type_lower)) return("Oligodendrocytes")
  if (grepl("^ast$", cell_type_lower)) return("Astrocyte")
  if (grepl("^opc$", cell_type_lower)) return("OPCs")
  if (grepl("^end$", cell_type_lower)) return("Endothelial")
  if (grepl("^mic$", cell_type_lower)) return("Macrophage/Microglia")
  if (grepl("^mix$", cell_type_lower)) return("Mixed/Unknown")
  return(as.character(raw_type))
}
raw_cell_type  <- sapply(raw_cell_type, map_cell_type)
seurat_obj@meta.data$cell.type <- raw_cell_type
seurat_obj@meta.data$sub_cell_type <- sub_cell_type
table(seurat_obj@meta.data$sub_cell_type)
save(seurat_obj,file="seurat.Rdata")
