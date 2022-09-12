#Optional: 
#Deletion of genes:
geneTCRBCR <- function(Data = 10X.data){
   TCRBCR.features <- c(grep('^TRAV', rownames(Data)),
                       grep('^TRAJ', rownames(Data)),
                       grep('^TRBV', rownames(Data)),
                       grep('^TRBD', rownames(Data)),
                     grep('^TRBJ', rownames(Data)))
    return(TCRBCR.features)
}

for(i in 1:length(datasets)){
features_to_delete <- geneTCRBCR(Data = datasets[[i]])
datasets[[i]] <- datasets[[i]][-features_to_delete,]
}

bs833_d0 <- CreateSeuratObject(counts = datasets[[i]], project = "bs833_d0", min.cells = 3, min.features = 200)

#if multiple datasets:
all_TIL <- merge(bs833_d0, y = c(dataset1,dataset2),
                          add.cell.ids = c("bs833_noMHC", "dataset1", 
                                    "dataset2"),
                          project = "all_TIL")
  
#Add mitchondrial info:
all_TIL[["percent.mt"]] <- PercentageFeatureSet(object = all_TIL, pattern = "^MT-")
head(all_TIL@meta.data, 5)

#Check QC measures:
#VlnPlot(all_TIL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(object = all_TIL, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(object = all_TIL, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1, plot2))
#nFeature_RNA: genes per cell; nCount_RNA: molecules detected per cell

#Subset cells that fulfill requirements:
all_TIL[[1]] <- subset(x = all_TIL[[1]], subset = nFeature_RNA > 150 & nFeature_RNA < 4500 & percent.mt < 10)


