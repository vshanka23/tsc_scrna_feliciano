##Load libraries and set max memory limits for R objects

library(Seurat)
library(dplyr)
library(rlang)
library(gtools)
options(future.globals.maxSize= 53687091200)

##Setup a variable to capture sample names from list of files in the working directory
sample_set <- gsub(".txt","",list.files()[grep(list.files(),pattern = ".txt")])

##Use the sample list to load the counts matrices into R/Seurat
for (x in sample_set) {
data_var <- paste("S",x,"_data",sep="")
data_dir_loc <- paste(x,"/outs/filtered_feature_bc_matrix/",sep = "")
assign(data_var,Read10X(data.dir=data_dir_loc))
}

##Convert counts data to Seurat single cell objects
for (x in sample_set) {
data_var <- get0(paste("S",x,"_data",sep=""))
project <- gsub("([2][4][7-9][1-8]_)","",x)
object <- paste("S",x,sep="")
assign(object,CreateSeuratObject(counts=data_var,project=project,min.cells=5))
}

##Remove dummy variable
rm(data_var)

##Set explanatory factors
S2476_FLFL$stim <- "FLFL"
S2478_FLFL$stim <- "FLFL"
S2484_FLFL$stim <- "FLFL"
S2486_FLFL$stim <- "FLFL"
S2488_WTWT$stim <- "WTWT"
S2491_WTWT$stim <- "WTWT"
S2492_WTWT$stim <- "WTWT"
S2493_WTWT$stim <- "WTWT"

S2476_FLFL$sample_id <- "2476_FLFL"
S2478_FLFL$sample_id <- "2478_FLFL"
S2484_FLFL$sample_id <- "2484_FLFL"
S2486_FLFL$sample_id <- "2486_FLFL"
S2488_WTWT$sample_id <- "2488_WTWT"
S2491_WTWT$sample_id <- "2491_WTWT"
S2492_WTWT$sample_id <- "2492_WTWT"
S2493_WTWT$sample_id <- "2493_WTWT"

##Check quality of sample using violin plots
VlnPlot(get0(paste("S",sample_set[1],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(get0(paste("S",sample_set[2],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(get0(paste("S",sample_set[3],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(get0(paste("S",sample_set[4],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(get0(paste("S",sample_set[5],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(get0(paste("S",sample_set[6],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(get0(paste("S",sample_set[7],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(get0(paste("S",sample_set[8],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)

##Apply standard Seurat filters to remove low quality cells (based on feature counts)
for (x in sample_set) {
sample <- paste("S",x,sep="")
assign(sample,subset(get0(sample),subset=nFeature_RNA > 300 & nFeature_RNA < 2500))
}

##Check quality again post filtering
VlnPlot(get0(paste("S",sample_set[1],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(get0(paste("S",sample_set[2],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(get0(paste("S",sample_set[3],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(get0(paste("S",sample_set[3],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(get0(paste("S",sample_set[4],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(get0(paste("S",sample_set[5],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(get0(paste("S",sample_set[6],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(get0(paste("S",sample_set[7],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(get0(paste("S",sample_set[8],sep="")),features=c("nFeature_RNA","nCount_RNA"), ncol=2)

for (x in sample_set) {
sample <- paste("S",x,sep="")
assign(sample,SCTransform(get0(sample),verbose = FALSE,return.only.var.genes = FALSE))
}

##Start the integration procedure
Integration_feature_set <-SelectIntegrationFeatures(object.list = c(S2476_FLFL,S2478_FLFL,S2484_FLFL,S2486_FLFL,S2488_WTWT,S2491_WTWT,S2492_WTWT,S2493_WTWT),nfeatures=1500)

Integration_list <- PrepSCTIntegration(object.list = c(S2476_FLFL,S2478_FLFL,S2484_FLFL,S2486_FLFL,S2488_WTWT,S2491_WTWT,S2492_WTWT,S2493_WTWT),anchor.features = Integration_feature_set, verbose=FALSE)

##Set the mode of integration to SCT (check 2017 paper for this procedure)
Integration_anchors <- FindIntegrationAnchors(object.list = Integration_list, normalization.method = "SCT", anchor.features = Integration_feature_set, verbose = FALSE)

##Normalize using SCT
Integrated <- IntegrateData(anchorset = Integration_anchors, normalization.method = "SCT", verbose = FALSE)

##Run dimensionality reduction procedures
Integrated <- RunPCA(object = Integrated, verbose = FALSE)

Integrated <- FindNeighbors(Integrated, reduction = "pca", dims = 1:10)

Integrated <- RunUMAP(Integrated, reduction = "pca", dims= 1:10)

##Iterate with different resolutions to parametrize granularity and identify cluster number plateau. Range of 0.4 to 2.0 is from Seurat documentation
for (x in seq(from=0.4, to=2.0, by=0.1)) {
cluster <- FindClusters(Integrated, resolution = x)
}

##0.9 identified. 1.0 actually took cluster number down from 24 to 21, so 0.9 is a local max.
Integrated <- FindClusters(Integrated, resolution =0.9)

##Add celltype IDs to a variable
Integrated$celltype <- Idents(Integrated)

##Add UMAP colored by different expanatory variables to plotting object
plots <- DimPlot(Integrated, reduction= "umap", group.by = c("stim","celltype","sample_id"),combine=FALSE)

##Plot and check (looks like integration worked after forcing cells number to 65K)
plots

##Create a new label variable with celltype IDs and group
Integrated$celltype.stim <- paste(Idents(Integrated),Integrated$stim,sep = "_")

##Prepare for marker identification and differential expression
##Seurat team recommends using SCT for integration but RNA mode for marker and DE analysis, so switching to RNA
DefaultAssay(Integrated) <- "RNA"
Integrated <- NormalizeData(Integrated,verbose=FALSE)

##Identify markers of all cluster based on which gene have enriched expression in each cluster compared the rest of the dataset
Integrated_markers_all <- FindAllMarkers(Integrated,min.pct = 0.1, logfc.threshold = 0.25, only.pos = TRUE)

##Select top 5 markers and write to a CSV file
Integrated_cluster_markers <- Integrated_markers_all %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC) %>% print(n=5*24)
write.csv(Integrated_cluster_markers,"scTSC_cluster_markers_top5.csv")

##Select top 10 markers and write to a CSV file
Integrated_cluster_markers <- Integrated_markers_all %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC) %>% print(n=10*24)
write.csv(Integrated_cluster_markers,"scTSC_cluster_markers_top10.csv")

##Switch to Integrated assay type to take advantage of pearson residuals generated from SCT for DE
DefaultAssay(Integrated) <- "integrated"
Idents(Integrated) <- "celltype.stim"

##Set factor labels
group <- c("WTWT","FLFL")
cluster <- c(0:(length(levels(Integrated$celltype))-1))

##Correct pearson residuals by back-regressing based on minimum median counts as covariats 
Integrated <- PrepSCTFindMarkers(Integrated)

##Run for loop to calculate DE per cluster
for (i in 0:(length(levels(Integrated$celltype))-1)) {
DE <- paste("DE","C",i,sep="_")
assign(DE,FindMarkers(Integrated,ident.1 = Index_FLFL[i+1], ident.2 = Index_WTWT[i+1],test.use = "MAST", assay = "SCT"))
}

##Consolidate DE results into a list of objects and write to xlsx
DE_list <- mget(mixedsort(ls(pattern="DE_C")))
write.xlsx(DE_list,file="DE_SCT_default-params.xlsx", rowNames=TRUE)

##Create a table consisting of cell numbers per cluster per identity
cell_numbers_per_cluster <- table(Integrated@active.ident, Integrated@meta.data$orig.ident)
write.xlsx(cell_numbers_per_cluster,"cluster_composition.xlsx",rowNames=TRUE)
