#@author: Li Na, SCU, SKLB

load(file='AD_data/seurat_ad.Rdata')

seurat_obj<-seurat_ad
table(seurat_obj$cell.type, seurat_obj$diagnosis)
set.seed(123)

# 1. Extract cells
pbmc <- subset(seurat_obj, subset = cell.type == "Oligodendrocyte")
pbmc <- subset(pbmc, subset = diagnosis == "DEM")
countMat = GetAssayData(pbmc, layer = "counts")
pbmc <- FindVariableFeatures(object=pbmc, selection.method="vst", nfeatures=10000)
hvgs <- VariableFeatures(pbmc)

target_gene="NEGR1"     
data=as.data.frame(countMat[unique(c(target_gene,hvgs)),])

#Select the top 10,000 cells with the highest expression values for analysis
cell_mean <- colMeans(data, na.rm = TRUE)
top_cell_idx <- order(cell_mean, decreasing = TRUE)[1:min(10000, length(cell_mean))]
data <- data[, top_cell_idx]
rm(pbmc)
rm(countMat)

#Perform virtual knockout
result <- scTenifoldKnk(countMatrix = data, 
                        gKO = target_gene,          
                        qc = TRUE,                    
                        qc_mtThreshold = 0.1,       
                        qc_minLSize = 1000,       
                        nc_nNet = 10,                
                        nc_nCells = 500             
)
head(result$diffRegulation) 
df=result$diffRegulation
df <- df[df$gene != target_gene, ]
df
write.table(df, file="allDiff.txt", sep="\t", quote=F, row.names=F)
outTab=df[df$p.adj<0.05,]
write.table(outTab, file="sigDiff.txt", sep="\t", quote=F, row.names=F)
outTab