seurat_object <- Read10X(data.dir = "/location/")
seurat_object <- CreateSeuratObject(counts = seurat_object, project = "name", min.cells = 3, min.features = 200)
seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^MT-")

VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
#nFeature_RNA: genes per cell; nCount_RNA: molecules detected per cell
seurat_object <- subset(x = seurat_object, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 10)
seurat_object <- SCTransform(seurat_object, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_object <- RunPCA(object = seurat_object, assay = "SCT", npcs = 50)
DimPlot(seurat_object, reduction = "pca")
seurat_object <- RunUMAP(object = seurat_object, assay = "SCT", return.model = TRUE, dims = 1:50)
seurat_object <- FindNeighbors(object = seurat_object, assay = "SCT", dims = 1:50)
seurat_object <- FindClusters(object = seurat_object, resolution = .5)

DimPlot(seurat_object, label = T)
FeaturePlot(seurat_object, features = kir_nk_marker)
