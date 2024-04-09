
# Set Libraries
library(dplyr)
library(Seurat)
library(ggplot2)

# Load datasets
zf102219.data <- Read10X(data.dir = 'V:/Ellie/scRNASeq/zf_102219_filtered/')
zf102219 <- CreateSeuratObject(counts = zf102219.data, min.cells = 3, min.features = 200, project = "nuc1")
zf120519.data <- Read10X(data.dir = 'V:/Ellie/scRNASeq/zf_120519_filtered/')
zf120519 <- CreateSeuratObject(counts = zf120519.data, min.cells = 3, min.features = 200, project = "nuc2")
zf013020.data <- Read10X(data.dir = 'V:/Ellie/scRNASeq/zf_013020_filtered/')
zf013020 <- CreateSeuratObject(counts = zf013020.data, min.cells = 3, min.features = 200, project = "nuc3")

# Merge datasets
comboRaw <- merge(zf102219, y = c(zf120519, zf013020), add.cell.ids = c("nuc1","nuc2","nuc3"), project = "comboCombo")

# Filter cells with unique feature counts >2500 or <200 and mt counts >1.5%
comboRaw[["percent.mt"]] <- PercentageFeatureSet(comboRaw, pattern = "^MT")
VlnPlot(comboRaw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
combo <- subset(comboRaw, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 1.5)
unique(sapply(X = strsplit(colnames(combo), split = "_"), FUN = "[", 1))

#Plot QC results
table(combo$orig.ident)
VlnPlot(combo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = .5) 

# Log-normalize whole dataset
combo <- NormalizeData(combo, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features of whole dataset
combo <- FindVariableFeatures(combo, selection.method = "vst", nfeatures = 2000)

# Scale whole dataset through linear transformation
all.genes <- rownames(combo)
combo <- ScaleData(combo, features = all.genes)

# Run PCA and find significant PCs of whole dataset
combo <- RunPCA(combo, npcs = 100)
DimPlot(combo, reduction = "pca")
ElbowPlot(combo, ndims = 100)
combo <- JackStraw(combo, num.replicate = 100, dims = 60)
combo <- ScoreJackStraw(combo, reduction = 'pca', dims = 1:60)
JackStrawPlot(combo, dims = 1:60)

# Determine PC cutoff 
goodDim = 55 

# Find nearest neighbors in graph and cluster whole dataset
combo <- FindNeighbors(combo, reduction = 'pca', dims = 1:goodDim)
combo <- FindClusters(combo, resolution = .1, n.start = 10) #used to be .1

# Visualize all clusters with UMAP
combo <- RunUMAP(combo, dims = 1:goodDim, min.dist = 0.75)
DimPlot(combo, reduction = "umap", pt.size = 1, label.size = 7, label = FALSE, group.by = "orig.ident")

# Plot signature markers to find general HVC cell types
FeaturePlot(combo, features = c("SYT1","SLC17A6","GAD2","SLC32A1","DCX","NR2E1","SLC15A2","CSF1R","PLP1","PDGFRA","SOX4","SPEF2","HBAD"), reduction = 'umap', label = FALSE, pt.size = .6, ncol = 3)
axis.title = element_text(size=20,face="bold")

# Extract GABAergic clusters 
comboInts <- subset(combo, subset = seurat_clusters == '4' | seurat_clusters == '6' | seurat_clusters == '14' | seurat_clusters == '16' | seurat_clusters == '18')

# Run PCA and find significant PCs of GABAergic subset
comboInts <- RunPCA(comboInts, npcs = 100)
DimPlot(comboInts, reduction = "pca")
ElbowPlot(comboInts, ndims = 100)
comboInts <- JackStraw(comboInts, num.replicate = 100, dims = 50)
comboInts <- ScoreJackStraw(comboInts, reduction = 'pca', dims = 1:50)
JackStrawPlot(comboInts, dims = 1:50)
intDims = 15

# Find nearest neighbors in graph and cluster GABAergic subset
comboInts <- FindNeighbors(comboInts, reduction = 'pca', dims = 1:intDims,do.plot = TRUE)
comboInts <- FindClusters(comboInts, resolution = .25, n.start = 10) #used to be .03 after .01
comboInts <- RunUMAP(comboInts, dims = 1:intDims, min.dist = 0.75)
DimPlot(comboInts, reduction = "umap", label = TRUE, pt.size = 1, label.size = 15)   #group.by = "orig.ident"

# Calculate silhouette metric to tweak cluster resolution above
library(cluster, quietly = TRUE)
reduction <- "pca"
dist.matrix <- dist(x = Embeddings(object = comboInts[[reduction]])[, 1:intDims])
clusters <- comboInts$seurat_clusters
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
comboInts$sil <- sil[, 3]
#windows()
plot(sil)

# List top differentiating markers of GABAergic clusters
comboInts.markers <- FindAllMarkers(comboInts, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
comboInts.markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_logFC) %>% View()

# Create heatmap using top 8 differentiating markers of GABAergic clusters
topList = c("MAF","NPY","TRHDE","RASL11A","PRSS12","ALKAL1","NKD1","SST","NXPH1","LOC115491088","LOC115493828","CYGB","ERBB4","CEMIP","LOC115496037","PPFIBP2","LTK","LOC115495115","TTLL5","RUNX2","ADSS1","IFITM10","UTRN","PVALB","ADARB2","PENK","PROX1","CRH","CNR1","NFIA","SCG3","NR2F2","CRTAC1","LOC100218415","CCN3", "ANO1","RPRM","QRFPR","CDH7","LOC100223410","CXCL12","PTN","CRHR2","CCK","NFIA","RGS16","TP53I11","CALB2","ZIC1","LHX8","ASIC4","KHDRBS2","IQSEC3","ETV1","SORCS1","LAPTM4B","MEIS2","PBX3","TMEM272","SIX3","TSHZ1","MEIS1","LOC105759036","ARC")
windows()
lev = levels(comboInts) <- c(1, 2, 5, 3, 4, 6, 7, 0)
comboInts$groupings = factor(comboInts$seurat_clusters,levels(comboInts))
DoHeatmap(comboInts, features = topList) + scale_fill_gradientn(colors = c("dark blue", "black", "yellow"))

# Visualize canonical markers that differentiate interneuron classes across species
shortList = c("MAF","SST","NPY","NXPH1","TTLL5","PVALB","ADARB2","PROX1","PENK","LAMP5","CALB2","CCK","RELN","LHX8","MEIS2","FOXP2")
windows()
#VlnPlot(comboInts, features = shortList, ncol = 2, group.by = 'seurat_clusters') 
DotPlot(comboInts, features = shortList, group.by = 'groupings',dot.scale = 10) + coord_flip()
list <- AverageExpression(comboInts, features = shortList,return.seurat = TRUE)

# Characterize MGE-derived cells with significant markers
mge <- subset(comboInts, subset = seurat_clusters == c('1','2','5'))
sum(comboInts$seurat_clusters == 1 | 2 | 5)
FeaturePlot(comboInts, features = c("PVALB","SST"),reduction = "umap",pt.size = 1, label.size = 15, cols = c("lightgrey", "red"))
#,"CHODL","CALB2",'NPY','GABRG1','GLRA3',"NOS1",'VIPR2','NKX2-1','IRES2','CPNE5','TAC1','MME','TPBG','ISLR','GPR149','SEMA3A'), ) 
colq = c('PVALB','SST')
DotPlot(mge, features = colq,  cols = c("lightgrey", "red"), scale.by = 'radius', dot.scale = 10, scale.min = 20)

# Compare datasets to other amniotes with Tosches et al. top markers
toschesList = c('PLAU','PVALB', 'ETV1', 'BCAN', 'LMO3', 'GABRB2', 'CNTNAP4', 'TRPS1', 'BTBD3', 'S100A1', 'KCNC1', 'KCNC2', 'ENO2', 'SOX6', 'LHX6', 'ARX', 'BCL11A', 'SATB1', 'SST', 'NPY', 'HTR1A', 'CACNG3', 'CTNND2', 'EFNA5', 'ELFN1', 'CBLN4', 'CALB1', 'DLX1', 'ZEB2', 'RELN', 'ID2', 'NDNF', 'RGS12', 'TNFAIP8L3', 'PENK', 'LRRTM4', 'SLC1A2', 'NR2F2', 'NPY1R', 'NPAS1', 'ADARB2', 'CNR1', 'SALL1', 'ZBTB16', 'HTR3A', 'NR2E1', 'PROX1', 'NFIX', 'GRP', 'SP8', 'VIP', 'TRH', 'MEIS2', 'PBX3', 'TSHZ1', 'SIX3', 'FOXP1', 'FOXP2', 'ZFHX4')
DoHeatmap(comboInts, features = toschesList) + scale_fill_gradientn(colors = c("dark blue", "black", "yellow"))


#saveRDS(combo, file = paste(dir,"combo.rds"))
#saveRDS(comboRaw, file = paste(dir, "comboRaw.rds"))
