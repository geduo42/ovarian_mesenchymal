library(Seurat)
library(SeuratObject)
library(dplyr)
library(patchwork)
library(ggplot2)


if(!file.exists(outdir)){
    dir.create(outdir)
}

pdf(paste(outdir, "seurat_luster.pdf", sep = "/"))

young.data   <- Read10X(data.dir = young_dir)

young   <- CreateSeuratObject(counts = young.data,   project = "young",   min.cells = 3, min.features = 200)

table(young$orig.ident)


young[["percent.mt"]] <- PercentageFeatureSet(young, pattern = "^mt-")
VlnPlot(young, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


save(young, file = paste(outdir, "young.raw.RData", sep = "/"))

plot1 <- FeatureScatter(young, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(young, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

young <- subset(young, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 15)

young <- NormalizeData(young, normalization.method = "LogNormalize", scale.factor = 10000)
young <- FindVariableFeatures(young, selection.method = "vst", nfeatures = 2000)

top20 <- head(VariableFeatures(young), 20)


plot1 <- VariableFeaturePlot(young)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2

all.genes <- rownames(young)
young <- ScaleData(young, features = all.genes)

young <- RunPCA(young )
print(young[["pca"]], dims = 1 : dims, nfeatures = 10)
VizDimLoadings(young, dims = 1 : dims, reduction = "pca")

DimPlot(young, reduction = "pca")
DimHeatmap(young, dims = 1 : dims, cells = 500, balanced = TRUE)

young <- FindNeighbors(young, dims = 1 : dims)
young <- FindClusters(young, resolution = res)
table(Idents(young))

table(young$seurat_clusters, young$orig.ident)
prop.table(table(Idents(young), young$orig.ident), margin = 2)
prop.table(table(young$seurat_clusters))

new.cluster.ids <- c("ETCLC", "Wnt4", "Col15a1+", "oMYB", "oMPs", "ITC", "Sars-", "IMC-1", "IMC-2")
names(new.cluster.ids) <- levels(young)

young <- RenameIdents(young, new.cluster.ids)

young <- RunTSNE(young, dims = 1 : dims, label = TRUE)
head(young@reductions$tsne@cell.embeddings)

tsnep1 <- DimPlot(young, reduction = "tsne")
tsnep1

table(young$seurat_clusters, young$orig.ident)

young <- RunUMAP(young, dims = 1 : dims, label = TRUE)
head(young@reductions$umap@cell.embeddings) 

DimPlot(young, reduction = "umap")

young.markers <- FindAllMarkers(young, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(young.markers, file = paste(outdir, "young_markers.txt", sep = "/"), sep = "\t", quote = F)

save(young, file = paste(outdir, "young.RData", sep = "/"))

young_sub <- subset(young, idents = c("Wnt4+", "ETCLC", "Col15a1+", "ITC", "oMYB", "oMPs"))
young_sub
levels(young_sub) <- c("oMYB", "Wnt4+", "Col15a1+", "ETCLC", "ITC", "oMPs")
levels(young_sub)

tSNEp_sub <- DimPlot(young_sub, reduction = "tsne", pt.size =1, cols = c("seashell4", "chocolate4", "peachpuff3", "lightskyblue", "dodgerblue", "hotpink"))
tSNEp_sub

young_sub.markers <- FindAllMarkers(young_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#View(young_sub.markers)
write.table(young_sub.markers, file = paste(outdir, "young_sub_markers.txt", sep = "/"), sep = "\t", quote = F)
save(young_sub, file = paste(outdir, "young_sub.RData", sep = "/"))

FeaturePlot(young_sub, features = c("Smim41", "Enpep", "Lamc3", "Cyp11a1", "Cyp17a1", "Hsd3b1", "Wt1", "Clec1a", "Galnt15", "Wnt4", "Myh11","Actg2", "Lmod1", "Col15a1", "Cadm4", "F2r", "Birc5", "Top2a", "Mki67", "Tcf21"), combine= F, reduction = "tsne",  pt.size = 1, min.cutoff = 0.5, max.cutoff = 1.5)

VlnPlot(young_sub, features = c("Ccnb2", "Thy1", "Ube2c", "Cks2"), combine = F, pt.size = 0, cols = c("seashell4", "chocolate4", "peachpuff3", "lightskyblue", "dodgerblue", "hotpink"))
prop.table(table(Idents(young_sub), young_sub$orig.ident), margin = 2)

a <- prop.table(table(Idents(young_sub), young_sub$orig.ident), margin = 2)

a <- as.data.frame.table(a)
write.table(a, file = paste(outdir, "young_sub_cell_proportion.txt", sep = "/"), sep = "\t", quote = F)


