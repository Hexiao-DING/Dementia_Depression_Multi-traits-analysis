#@author: Li Na, SCU, SKLB
load("Data/seurat.Rdata")
SaveH5Seurat(seurat_obj, filename = "output.h5Seurat")
Convert("output.h5Seurat", dest = "h5ad")
Idents(seurat_obj) <- "cell.type"
cell.use <- colnames(seurat_obj_noPDD)[seurat_obj_noPDD$diagnosis == "DEM"]
seurat_cellchat <- subset(seurat_obj_noPDD, cells = cell.use)
data.input <- GetAssayData(seurat_cellchat, assay = "RNA", layer = "data")
meta <- seurat_cellchat@meta.data
meta <- data.frame(
  labels = meta$cell.type,
  diagnosis = meta$diagnosis,
  group = meta$group,
  row.names = colnames(data.input)
)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat@DB <- CellChatDB.human

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
dim(cellchat@net$weight)
cellchat@net$weight
table(cellchat@idents)
net <- cellchat@net$weight
groupSize <- table(cellchat@idents)
groupSize <- as.numeric(groupSize[rownames(net)])
netVisual_circle(
  net,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE
)
cellchat@netP$pathways

pathways.show <- c("NEGR")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#------------------------cellchat_ctrl_net-----------------------
seob_list <- SplitObject(seurat_obj, split.by = "diagnosis")
cellchat_ctrl <- createCellChat(object = seob_list$Control, group.by = "cell.type")
cellchat_stim <- createCellChat(object = seob_list$MDD, group.by = "cell.type")
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
# Set up the ligand-receptor database
cellchat_ctrl@DB <- CellChatDB
cellchat_stim@DB <- CellChatDB
cellchat_ctrl <- subsetData(cellchat_ctrl)
cellchat_stim <- subsetData(cellchat_stim)
# Calculate communication probability
cellchat_ctrl <- identifyOverExpressedGenes(cellchat_ctrl)
cellchat_ctrl <- identifyOverExpressedInteractions(cellchat_ctrl)
cellchat_ctrl <- computeCommunProb(cellchat_ctrl)
cellchat_stim <- identifyOverExpressedGenes(cellchat_stim)
cellchat_stim <- identifyOverExpressedInteractions(cellchat_stim)
cellchat_stim <- computeCommunProb(cellchat_stim)
# Filter out the communication signals with low confidence levels
cellchat_ctrl <- filterCommunication(cellchat_ctrl, min.cells = 10)
cellchat_stim <- filterCommunication(cellchat_stim, min.cells = 10)
# Calculate network centrality
cellchat_ctrl <- netAnalysis_computeCentrality(cellchat_ctrl, slot.name = "net")
cellchat_stim <- netAnalysis_computeCentrality(cellchat_stim, slot.name = "net")
# Data frame
cellchat_ctrl_net <- subsetCommunication(cellchat_ctrl, slot.name = "net")
cellchat_stim_net <- subsetCommunication(cellchat_stim, slot.name = "net")
head(cellchat_stim_net)
write_xlsx(cellchat_stim_net, path = "cellchat_mdd_net.xlsx")
#-------------------Overall communication network-----------------------
# Clustering cell communication network
cellchat_ctrl <- computeCommunProbPathway(cellchat_ctrl)
cellchat_ctrl <- aggregateNet(cellchat_ctrl)
cellchat_stim <- computeCommunProbPathway(cellchat_stim)
cellchat_stim <- aggregateNet(cellchat_stim)
cellchat_ctrl_netP <- subsetCommunication(cellchat_ctrl, slot.name = "netP")
cellchat_stim_netP <- subsetCommunication(cellchat_stim, slot.name = "netP")
# netVisual_circle
netVisual_circle(cellchat_stim@net$weight, edge.width.max = 10)
# netVisual_heatmap
netVisual_heatmap(cellchat_stim, measure = "weight", color.heatmap = "Reds")
pathways.show <- c("NEGR")
netVisual_aggregate(cellchat_stim, signaling = pathways.show, layout = "circle")
netVisual_heatmap(cellchat_stim, signaling = pathways.show, color.heatmap = "Reds")
#-------------------Network centrality-----------------------
cellchat_ctrl <- netAnalysis_computeCentrality(cellchat_ctrl, slot.name = "netP")
cellchat_stim <- netAnalysis_computeCentrality(cellchat_stim, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(
  cellchat_stim, 
  pattern = "outgoing",
  width = 8,
  height = 12,
  font.size = 8,
  font.size.title = 10
)
ht2 <- netAnalysis_signalingRole_heatmap(
  cellchat_stim,
  pattern = "incoming",
  width = 8,
  height = 12,
  font.size = 8,
  font.size.title = 10
)
ht1 + ht2
#draw(ht1 + ht2, ht_gap = unit(8, "mm"))
pathway_name <- "NEGR"
gg1 <- netAnalysis_signalingRole_scatter(cellchat_stim)
gg2 <- netAnalysis_signalingRole_scatter(cellchat_stim, signaling = pathway_name)
gg1 + gg2
# -------------------Evaluation of communication mode: Similarity among cells------------------
# outgoing 和 incoming
library(NMF)
library(ggalluvial)
selectK(cellchat_stim, pattern = "outgoing")
dev.off()
# When the number of extroverted patterns reaches 6, both the Coffini coefficient and the contour coefficient suddenly decrease.
nPatterns = 6
cellchat_stim <- identifyCommunicationPatterns(cellchat_stim, pattern = "outgoing", k = nPatterns,
                                               width = 6,
                                               height = 10,
                                               font.size = 8)
netAnalysis_river(cellchat_stim, pattern = "outgoing")
netAnalysis_dot(cellchat_stim, pattern = "outgoing")
