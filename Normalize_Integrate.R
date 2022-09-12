#Split merged Seurat object into a list of objects (i.e. each list item is one patient dataset)
seurat_object.list <- SplitObject(all_TIL, split.by = "orig.ident")
#Normalize each dataset via scTransform (fast and accurate, may want to run it on the cluster)
for (i in 1:length(seurat_object.list)) {
    seurat_object.list[[i]] <- SCTransform(seurat_object.list[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
}

# Save variable features for later integration:
TIL.features <- SelectIntegrationFeatures(object.list = seurat_object.list, nfeatures = 3000)
# Merge data and run PCA:
all_TIL.sct <- merge(seurat_object.list[[1]], y = seurat_object.list[2:length(seurat_object.list)], project = "all_TIL", merge.data = TRUE)

#Computationally intensive, thus made on cluster:
VariableFeatures(all_TIL.sct) <- TIL.features
all_TIL.sct <- RunPCA(object = all_TIL.sct, assay = "SCT", npcs = 50)
DimPlot(all_TIL.sct, reduction = "pca", group.by = 'orig.ident')#Sanity check of PCA
#If integrating datasets, harmony corrects for batch effects:
all_TIL.sct.h <- RunHarmony(object = all_TIL.sct, 
                                    assay.use = "SCT",
                                   reduction = "pca",
                                    dims.use = 1:50,
                                    group.by.vars = "orig.ident",
                                    plot_convergence = TRUE)
all_TIL.sct.h <- RunUMAP(object = all_TIL.sct.h, assay = "SCT", reduction = "harmony", return.model = TRUE, dims = 1:50)
all_TIL.sct.h <- FindNeighbors(object = all_TIL.sct.h, assay = "SCT", reduction = "harmony", dims = 1:50)
all_TIL.sct.h <- FindClusters(object = all_TIL.sct.h, resolution = .5)

#Check results of neighbourhood analysis & clustering:
DimPlot(all_TIL.sct.h, label = T)
