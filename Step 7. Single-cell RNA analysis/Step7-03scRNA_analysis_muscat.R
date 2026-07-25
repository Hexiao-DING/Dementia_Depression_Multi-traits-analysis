load(file='Data/seurat_obj_noPDD.Rdata')
obj <- seurat_obj_noPDD
DimPlot(obj,reduction = 'umap',group.by = c('orig.ident','diagnosis','cell.type'))  
obj <- obj[rowSums(GetAssayData(obj,layer = 'counts') > 0) > 10, ]  
sce <- as.SingleCellExperiment(obj)

(sce <- prepSCE(sce,   
                kid = 'cell.type',    # subpopulation assignments  
                gid = 'diagnosis',       # group IDs (ctrl/stim)  
                sid = 'sample',      # sample IDs (ctrl/stim.1234)  
                drop = TRUE))        # drop all other colData columns  

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))  
names(kids) <- kids; names(sids) <- sids  

t(table(sce$cluster_id, sce$sample_id)) 
sce <- runUMAP(sce, pca = 30)  
# wrapper to prettify reduced dimension plots  
.plot_dr <- function(sce, dr, col)  
  plotReducedDim(sce, dimred = dr, colour_by = col) +  
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +  
  theme_minimal() + theme(aspect.ratio = 1)  

# downsample to max. 100 cells per cluster
cs_by_k <- split(colnames(sce), sce$cluster_id)  
cs100 <- unlist(sapply(cs_by_k, function(u)   
  sample(u, min(length(u), 100))))  

# plot t-SNE & UMAP colored by cluster & group ID  
pdf(file = '0.Dataoverview_downsample100.pdf')  
for (dr in c( "UMAP"))  
  for (col in c("cluster_id", "group_id"))  
    p = .plot_dr(sce[, cs100], dr, col)  
print(p)  
dev.off() 

sce2 <- sce[, sce$cluster_id != "Ependymal cell"]
#sce2$cluster_id <- droplevels(sce2$cluster_id)
#@ct_tab <- table(sce2$sample_id, sce2$cluster_id)
#bad_idx <- which(ct_tab < 10, arr.ind = TRUE)
#bad_pairs <- data.frame(sample_id = rownames(ct_tab)[bad_idx[, 1]],cluster_id = colnames(ct_tab)[bad_idx[, 2]],stringsAsFactors = FALSE)
#cell_pair <- paste(sce2$sample_id, sce2$cluster_id, sep = "___")
#bad_pair_names <- paste(bad_pairs$sample_id, bad_pairs$cluster_id, sep = "___")
#sce2 <- sce2[, !(cell_pair %in% bad_pair_names)]
#sce2$cluster_id <- droplevels(sce2$cluster_id)
#sce2$sample_id <- droplevels(sce2$sample_id)
#t(table(sce2$cluster_id, sce2$sample_id))

pb <- aggregateData(
  sce2,
  assay = "counts",
  fun = "sum",
  by = c("cluster_id", "sample_id")
)
# one sheet per subpopulation  
assayNames(pb)  
t(head(assay(pb)))   
(pb_mds <- pbMDS(pb))  
ggsave(filename = '1.data_aggr_celltype_sample.png',print(pb_mds),w=7,h=7)  
ggsave(filename = '1.data_aggr_celltype_sample.pdf',print(pb_mds),w=7,h=7) 

## run DS analysis##
library(UpSetR)
ei <- metadata(sce2)$experiment_info  
ei
g_n = levels(ei$group_id)  
mm <- model.matrix(~ 0 + ei$group_id)  
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))  
contrast <- makeContrasts(paste0(g_n[1],'-',g_n[2]), levels = mm)  
# run DS analysis  
res <- pbDS(pb, verbose = FALSE, design = mm, contrast = contrast,filter = 'none') 
group_ALL <- bind_rows(res$table[[1]])  
write.csv(file = paste0('2.allgene_nofilter.csv'),group_ALL)  
candidate_file <- "D:/Projects_data&code/Stage1_bioinformatics_ADandDepression/STRING/candidate_genes_TeQ.txt"
candidates <- readLines(candidate_file)
tbl_fil <- group_ALL %>% 
  dplyr::filter(gene %in% candidates)%>%  
  dplyr::filter(p_adj.loc < 0.05, abs(logFC) > 1) %>%  
  dplyr::arrange(cluster_id,desc(logFC)) %>% ungroup() 
write.csv(file = paste0('2.GroupDEG_fil.csv'),tbl_fil)  
###visualizing  
de_gs_by_k <- split(tbl_fil$gene, tbl_fil$cluster_id)  
p1 = upset(fromList(de_gs_by_k))  
p1
ggsave(filename = '3.data_aggr_celltype_sample.png',print(p1),w=7,h=7)  
ggsave(filename = '3.data_aggr_celltype_sample.pdf',print(p1),w=7,h=7)   

