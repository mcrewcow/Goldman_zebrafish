pbmc <- RenameIdents(pbmc, '1' = 'Resting MG', '2' = 'Activated MG', '3' = "Progenitors", '4' = 'Rods', '5' = 'Cone  BC', '6' = 'BC',
                     '7' = 'GABAergic AC', '8' = 'Glycinergic AC', '9' = 'Cones', '10' = 'RGC', '11' = 'Microglia',
                     '12' = 'V/E cells', '13' = 'HC', '14' = 'RPE', '15' = 'Pericytes')
pbmc$SB_anno  <- pbmc@active.ident

ProcessInt <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mito',"percent.ribo","S.Score","G2M.Score"))
  data.integrated <- RunPCA(data.integrated, npcs = 75, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:75)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:75)
}

pbmc <- Seurat.utils::RenameGenesSeurat(pbmc, newnames = genes$Symbol)
remove_rows <- which(table(rownames(pbmc)) > 1)
pbmc <- pbmc[-remove_rows, ]

pbmc[["RNA"]]@meta.features <- data.frame(row.names = rownames(pbmc[["RNA"]]))
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ProcessInt(pbmc)

zfAd00 <- subset(pbmc, subset = Condition == 'zfAd00')
zfNMDA04  <- subset(pbmc, subset = Condition == 'zfNMDA04')
zfNMDA10  <- subset(pbmc, subset = Condition == 'zfNMDA10')
zfNMDA20 <- subset(pbmc, subset = Condition == 'zfNMDA20')
zfNMDA36  <- subset(pbmc, subset = Condition == 'zfNMDA36')

ProcessSeu <- function(Seurat){
  Seurat <- RunPCA(Seurat, npcs = 75)
  Seurat <- FindNeighbors(Seurat, dims = 1:75)
  Seurat <- FindClusters(Seurat, resolution = 1)
  Seurat <- RunUMAP(Seurat, dims = 1:75)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:75)
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}

counts <- Matrix::readMM('C://Users/rodri/Downloads/Zebrafish_LD_count_matrix.mtx')
genes <- readr::read_tsv('C://Users/rodri/Downloads/Zebrafish_gene_feature.tsv', )
features <- readr::read_tsv('C://Users/rodri/Downloads/Zebrafish_LD_cell_feature.tsv')
colnames(counts) <- features$`#Barcode`
rownames(counts) <- genes$Symbol
LD <- CreateSeuratObject(counts = counts, project = "hippo_EK", min.cells = 3, min.features = 200)
LD$Condition <- 'LD'

ProcessSeu <- function(Seurat){
  Seurat <- RunPCA(Seurat, npcs = 75)
  Seurat <- FindNeighbors(Seurat, dims = 1:75)
  Seurat <- FindClusters(Seurat, resolution = 1)
  Seurat <- RunUMAP(Seurat, dims = 1:75)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:75)
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}
LD <- FindVariableFeatures(LD, selection.method = "vst", nfeatures = 2000)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "www", host = "dec2021.archive.ensembl.org")
  mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "www", host = "dec2021.archive.ensembl.org")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)

LD$percent.mito <- features$Percentage.of.mitochondrial.genes
LD$percent.ribo <- features$Percentage.of.ribosomal.protein.genes
LD <- CellCycleScoring(LD, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = FALSE)

LD <- ScaleData(LD, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mito',"percent.ribo","S.Score","G2M.Score"))
LD <- ProcessSeu(LD)

zfAd00 <- ProcessSeu(zfAd00)
zfNMDA04 <- ProcessSeu(zfNMDA04)
zfNMDA10 <- ProcessSeu(zfNMDA10)
zfNMDA20 <- ProcessSeu(zfNMDA20)
zfNMDA36 <- ProcessSeu(zfNMDA36)

integration_list <- list(zfAd00, zfNMDA04, zfNMDA10, zfNMDA20, zfNMDA36)
features <- SelectIntegrationFeatures(object.list = integration_list)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)
SB_zebrafish_integrated <- IntegrateData(anchorset = data.anchors)
SB_zebrafish_integrated <- ProcessInt(SB_zebrafish_integrated)

data.combined.markers <- FindAllMarkers(SB_zebrafish_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

data.combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
data.combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(SB_zebrafish_integrated, features = top10$gene) + NoLegend()


SaveH5Seurat(SB_zebrafish_integrated, 'C://Bioinf/Blackshaw_zebrafish_EK_reintegrated.h5Seurat')
