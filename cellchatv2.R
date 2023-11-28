library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(CellChat)

SB_zebrafish <- LoadH5Seurat('C://Bioinf/Blackshaw_zebrafish_Hoang_science/zebrafish_Hoang_2023_EK_anno.h5Seurat')

SV_zebrafish_Adult <- subset(SB_zebrafish, subset = Condition == 'zfAd00')
SV_zebrafish_NMDA4 <- subset(SB_zebrafish, subset = Condition == 'zfNMDA04')
SV_zebrafish_NMDA10 <- subset(SB_zebrafish, subset = Condition == 'zfNMDA10')
SV_zebrafish_NMDA20 <- subset(SB_zebrafish, subset = Condition == 'zfNMDA20')
SV_zebrafish_NMDA36 <- subset(SB_zebrafish, subset = Condition == 'zfNMDA36')

object.list <- list(Control = cellchat.NL, LS = cellchat.LS)

Adult_cellchat <- readRDS('C://Bioinf/Blackshaw_zebrafish_Hoang_science/zebrafish_Hoang_2023_Adult_CellChat.rds')
NMDA4_cellchat <- readRDS('C://Bioinf/Blackshaw_zebrafish_Hoang_science/zebrafish_Hoang_2023_NMDA4_CellChat.rds')
NMDA10_cellchat <- readRDS('C://Bioinf/Blackshaw_zebrafish_Hoang_science/zebrafish_Hoang_2023_NMDA10_CellChat.rds')
NMDA20_cellchat <- readRDS('C://Bioinf/Blackshaw_zebrafish_Hoang_science/zebrafish_Hoang_2023_NMDA20_CellChat.rds')
NMDA36_cellchat <- readRDS('C://Bioinf/Blackshaw_zebrafish_Hoang_science/zebrafish_Hoang_2023_NMDA36_CellChat.rds')

object.list <- list(Control = Adult_cellchat, NMDA4 = NMDA4_cellchat, NMDA10 = NMDA10_cellchat, NMDA20 = NMDA20_cellchat, NMDA36 = NMDA36_cellchat)

MERGED_cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(MERGED_cellchat, show.legend = F, group = c(1:5))
gg2 <- compareInteractions(MERGED_cellchat, show.legend = F, group = c(1:5), measure = "weight")
gg1 + gg2

gg1 <- rankNet(MERGED_cellchat, mode = "comparison", comparison = c(1:5), stacked = T, do.stat = TRUE)
gg2 <- rankNet(MERGED_cellchat, mode = "comparison", comparison = c(1:5),stacked = F, do.stat = TRUE)
gg1 + gg2

library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 3, height = 8)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+4]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 3, height = 8)
draw(ht1 + ht2, ht_gap = unit(3, "cm"))


pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 3, height = 8)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+4]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 3, height = 8)
draw(ht1 + ht2, ht_gap = unit(3, "cm"))

netVisual_bubble(MERGED_cellchat, sources.use = 4, targets.use = c(1:11),  comparison = c(1:5), angle.x = 45)
#> Comparing communications on a merged object
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Control"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS and cellchat@var.features$LS.info. 
cellchat <- identifyOverExpressedGenes(MERGED_cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05, thresh.p = 1) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "Control",ligand.logFC = 0.1, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = c('NMDA4','NMDA10','NMDA20','NMDA36'),ligand.logFC = -0.05, receptor.logFC = -0.05)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1:5),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1:5),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

install.packages('wordcloud')
computeEnrichmentScore(net.down, species = 'human')
computeEnrichmentScore(net.up, species = 'human')
