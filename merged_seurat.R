library(Seurat)
library(SeuratObject)
library(dplyr)
library(patchwork)
library(ggplot2)

if(!file.exists(outdir)){
    dir.create(outdir)
}

pdf(paste(outdir, "young_aged_cluster.pdf", sep = "/"))

aged.data   <- Read10X(data.dir = aged_dir)
young.data <- Read10X(data.dir = young_dir)

aged   <- CreateSeuratObject(counts = aged.data,   project = "aged",   min.cells = 3, min.features = 200)
young <- CreateSeuratObject(counts = young.data, project = "young", min.cells = 3, min.features = 200)

merged <- merge( aged, y = c(young), add.cell.ids = c("aged","young"))
table(merged$orig.ident)


merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^mt-")
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

save(merged, file = paste(outdir, "merged.RData", sep = "/"))

plot1 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

merged<- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 15)

merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)

top20 <- head(VariableFeatures(merged), 20)

plot1 <- VariableFeaturePlot(merged)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2

all.genes <- rownames(merged)
merged <- ScaleData(merged, features = all.genes)

merged <- RunPCA(merged )
print(merged[["pca"]], dims = 1 : 30, nfeatures = 10)
VizDimLoadings(merged, dims = 1 : 30, reduction = "pca")

DimPlot(merged, reduction = "pca")
DimHeatmap(merged, dims = 1 : 30, cells = 500, balanced = TRUE)

merged <- JackStraw(merged, num.replicate = 100)
merged <- ScoreJackStraw(merged, dims = 1 : 30)
JackStrawPlot(merged, dims = 1 : 30)

ElbowPlot(merged, ndims = 30)

merged <- FindNeighbors(merged, dims = 1 : dims)
merged <- FindClusters(merged, resolution = res)
table(Idents(merged))

table(merged$seurat_clusters, merged$orig.ident)
prop.table(table(Idents(merged), merged$orig.ident), margin = 2)
prop.table(table(merged$seurat_clusters))


cluster_num <- length(levels(merged))
new.cluster.ids <- c("Wnt4+", "ETCLC", "Col15a1+", "oMYB", "oMPs", "ITC", "Sars-", "IMC-1", "IMC-2") 
names(new.cluster.ids) <- levels(merged)
merged <- RenameIdents(merged, new.cluster.ids)


merged <- RunTSNE(merged, dims = 1 : dims, label = TRUE)
head(merged@reductions$tsne@cell.embeddings)
tsnep1 <- DimPlot(merged, reduction = "tsne",group.by = "orig.ident")
tsnep1

tsnep2 <- DimPlot(merged, reduction = "tsne")
tsnep2

FeaturePlot(merged, features = c("Smim41", "Enpep", "Cyp11a1", "Cyp17a1", "Wt1", "Clec1a", "Actg2", "Lmod1", "Cadm4", "F2r", "Top2a", "Mki67"), combine= F, reduction = "tsne",  pt.size = 0.65)

head(merged@reductions$tsne@cell.embeddings)
tSNEp1 <- DimPlot(merged, reduction = "tsne")
tSNEp1 
table(merged$seurat_clusters, merged$orig.ident)


merged <- RunUMAP(merged, dims = 1 : dims, label = TRUE)
head(merged@reductions$umap@cell.embeddings) 
UMAPp2 <- DimPlot(merged, reduction = "umap")
UMAPp2

merged.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(merged.markers, file = paste(outdir, "young_aged_markers.txt", sep = "/"), sep = "\t", quote = F)

save(merged, file = paste(outdir, "merged.RData", sep = "/"))

merged[['cell_type']] <- merged@active.ident
levels(merged)
merged_sub <- subset(merged, idents = c("Wnt4+", "ETCLC", "Col15a1+", "ITC", "oMYB", "oMPs"))
merged_sub
levels(merged_sub) <- c("oMYB", "Wnt4+", "Col15a1+", "ETCLC", "ITC", "oMPs")
levels(merged_sub)


tSNEp_sub <- DimPlot(merged_sub, reduction = "umap", pt.size =1.5, cols = c("chocolate4", "lightskyblue", "peachpuff3", "seashell4", "hotpink", "dodgerblue"))
tSNEp_sub

prop.table(table(Idents(merged_sub), merged_sub$orig.ident), margin = 2)

merged_sub.markers <- FindAllMarkers(merged_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(merged_sub.markers, file = paste(outdir, "y2_o1_sub_markers.txt", sep = "/"), sep = "\t", quote = F)
save(merged_sub, file = paste(outdir, "merged_sub.RData", sep = "/"))




