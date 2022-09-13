multi <- readRDS("multiome_file.rds")

anchors <- FindTransferAnchors(
    reference = multi,
    query = seurat_object,
    k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50,recompute.residuals = FALSE)
seurat_object <- MapQuery(
    anchorset = anchors, 
    query = seurat_object,
    reference = multi,
    reference.reduction = "spca",
    reduction.model = "wnn.umap")

Idents(seurat_object) <- "predicted.celltype"
DimPlot(seurat_object, reduction = 'ref.umap', label=T)
p1 <- DimPlot(multi, reduction = 'umap', label = T)
p2 <- DimPlot(seurat_object, reduction = 'ref.umap', label = T)
p1+p2
