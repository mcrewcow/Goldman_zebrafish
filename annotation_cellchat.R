table(SV_zebrafish@active.ident)

SV_zebrafish <- RenameIdents(SV_zebrafish, '0' = 'Transitory MG',
                             '1' = 'Cone BC',
                             '2' = 'RGC',
                             '3' = 'Resting MG',
                             '4' = 'Cones',
                             '5' = 'Cones',
                             '6' = 'Interneurons',
                             '7' = 'Rods',
                             '8' = 'GABAergic Amacrine',
                             '9' = 'Interneurons',
                             '10' = 'Bipolar',
                             '11' = 'Transitory MG',
                             '12' = 'Microglia/Macrophages',
                             '13' = 'Resting MG',
                             '14' = 'Activated MG',
                             '15' = 'Activated MG',
                             '16' = 'Pericytes',
                             '17' = 'Glycinergic Amacrine',
                             '18' = 'Microglia/Macrophages',
                             '19' = 'Activated MG',
'20' = 'GABAergic Amacrine',
'21' = 'Cones',
'22' = 'Interneurons',
'23' = 'V/E',
'24' = 'Oligodendrocytes',
'25' = 'Microglia/Macrophages',
'26' = 'Pericytes',
'27' = 'Activated MG',
'28' = 'RPE',
'29' = 'Lens',
'30' = 'Cornea',
'31' = 'Cornea',
'32' = 'Red blood cells',
'33' = 'Red blood cells')

SV_zebrafish$EK_anno_2023 <- SV_zebrafish@active.ident



DimPlot(SV_zebrafish, group.by = 'seurat_clusters', label = T, repel = T, label.box = T)

DimPlot(SV_zebrafish, group.by = 'EK_anno_2023', label = T, label.box = T)

SV_zebrafish <- subset(SV_zebrafish, subset = EK_anno_2023 == c('Red blood cells'), invert = T)

SV_zebrafish <- subset(SV_zebrafish, subset = EK_anno_2023 == c('Lens'), invert = T)

SV_zebrafish <- subset(SV_zebrafish, subset = EK_anno_2023 == c('Cornea'), invert = T)

SV_zebrafish1 <- ProcessInt(SV_zebrafish)

SaveH5Seurat(SV_zebrafish1, 'C://Bioinf/Blackshaw_zebrafish_Hoang_science/zebrafish_Hoang_2023_EK_anno.h5Seurat')
saveRDS(SV_zebrafish1, 'C://Bioinf/Blackshaw_zebrafish_Hoang_science/zebrafish_Hoang_2023_EK_anno.rds')

SV_zebrafish1$EK_anno_2023 <- droplevels(SV_zebrafish1$EK_anno_2023, exclude = setdiff(levels(SV_zebrafish1$EK_anno_2023), unique(SV_zebrafish1$EK_anno_2023)))
SV_zebrafish1@active.ident <- droplevels(SV_zebrafish1@active.ident, exclude = setdiff(levels(SV_zebrafish1@active.ident), unique(SV_zebrafish1@active.ident)))

SV_zebrafish_Adult <- subset(SV_zebrafish1, subset = Condition == 'zfAd00')
SV_zebrafish_NMDA4 <- subset(SV_zebrafish1, subset = Condition == 'zfNMDA04')
SV_zebrafish_NMDA10 <- subset(SV_zebrafish1, subset = Condition == 'zfNMDA10')
SV_zebrafish_NMDA20 <- subset(SV_zebrafish1, subset = Condition == 'zfNMDA20')
SV_zebrafish_NMDA36 <- subset(SV_zebrafish1, subset = Condition == 'zfNMDA36')

library(CellChat)
cellchat <- createCellChat(object = SV_zebrafish_Adult, group.by = "EK_anno_2023", assay = 'RNA')
CellChatDB <- CellChatDB.zebrafish
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = TRUE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

ht1 <- netAnalysis_signalingRole_heatmap(SV_zebrafish_Adult_cellchat, pattern = "outgoing", width = 10)
ht2 <- netAnalysis_signalingRole_heatmap(SV_zebrafish_Adult_cellchat, pattern = "incoming", width = 10)
ht1 + ht2


SV_zebrafish_Adult_cellchat <- cellchat

cellchat <- createCellChat(object = SV_zebrafish_NMDA4, group.by = "EK_anno_2023", assay = 'RNA')
CellChatDB <- CellChatDB.zebrafish
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = TRUE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

ht1 <- netAnalysis_signalingRole_heatmap(SV_zebrafish_NMDA4_cellchat, pattern = "outgoing", width = 10)
ht2 <- netAnalysis_signalingRole_heatmap(SV_zebrafish_NMDA4_cellchat, pattern = "incoming", width = 10)
ht1 + ht2

SV_zebrafish_NMDA4_cellchat <- cellchat

cellchat <- createCellChat(object = SV_zebrafish_NMDA10, group.by = "EK_anno_2023", assay = 'RNA')
CellChatDB <- CellChatDB.zebrafish
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = TRUE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

ht1 <- netAnalysis_signalingRole_heatmap(SV_zebrafish_NMDA10_cellchat, pattern = "outgoing", width = 10)
ht2 <- netAnalysis_signalingRole_heatmap(SV_zebrafish_NMDA10_cellchat, pattern = "incoming", width = 10)
ht1 + ht2

SV_zebrafish_NMDA10_cellchat <- cellchat

cellchat <- createCellChat(object = SV_zebrafish_NMDA20, group.by = "EK_anno_2023", assay = 'RNA')
CellChatDB <- CellChatDB.zebrafish
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = TRUE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

ht1 <- netAnalysis_signalingRole_heatmap(SV_zebrafish_NMDA20_cellchat, pattern = "outgoing", width = 10, height = 12)
ht2 <- netAnalysis_signalingRole_heatmap(SV_zebrafish_NMDA20_cellchat, pattern = "incoming", width = 10, height = 12)
ht1 + ht2

SV_zebrafish_NMDA20_cellchat <- cellchat

cellchat <- createCellChat(object = SV_zebrafish_NMDA36, group.by = "EK_anno_2023", assay = 'RNA')
CellChatDB <- CellChatDB.zebrafish
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = TRUE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

ht1 <- netAnalysis_signalingRole_heatmap(SV_zebrafish_NMDA36_cellchat, pattern = "outgoing", width = 10)
ht2 <- netAnalysis_signalingRole_heatmap(SV_zebrafish_NMDA36_cellchat, pattern = "incoming", width = 10)
ht1 + ht2

SV_zebrafish_NMDA36_cellchat <- cellchat

saveRDS(SV_zebrafish_Adult_cellchat, 'C://Bioinf/Blackshaw_zebrafish_Hoang_science/zebrafish_Hoang_2023_Adult_CellChat.rds')
saveRDS(SV_zebrafish_NMDA4_cellchat, 'C://Bioinf/Blackshaw_zebrafish_Hoang_science/zebrafish_Hoang_2023_NMDA4_CellChat.rds')
saveRDS(SV_zebrafish_NMDA10_cellchat, 'C://Bioinf/Blackshaw_zebrafish_Hoang_science/zebrafish_Hoang_2023_NMDA10_CellChat.rds')
saveRDS(SV_zebrafish_NMDA20_cellchat, 'C://Bioinf/Blackshaw_zebrafish_Hoang_science/zebrafish_Hoang_2023_NMDA20_CellChat.rds')
saveRDS(SV_zebrafish_NMDA36_cellchat, 'C://Bioinf/Blackshaw_zebrafish_Hoang_science/zebrafish_Hoang_2023_NMDA36_CellChat.rds')

