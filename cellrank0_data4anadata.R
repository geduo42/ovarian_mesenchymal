ibrary(Seurat)
library(Matrix)

setwd("PATH")


args <- commandArgs()

if(length(args) < 6){    print("please input at least 1 parameter, parameter 1 for outdir  ....  ")
    q()
    }

path <- args[6]
outdir <- paste(path, "/data4anadata", sep = "")
print(outdir)

if(!file.exists(outdir)){
    dir.create(outdir)
}

# load data
load(paste(path, "/young.RData", sep = ""))
sc <- young

new.cluster.ids <- c("oMYB", "Wnt4", "Col15a1", "ETCLC", "ITC", "oMPs")
names(new.cluster.ids) <- levels(sc)
sc <- RenameIdents(sc, new.cluster.ids)
levels(sc)<- c("oMYB", "Wnt4", "Col15a1", "ETCLC", "ITC", "oMPs")
sc[['cell_type']] <- sc@active.ident
table(sc$cell_type)





# save metadata table
if(F){
sc$barcode <- colnames(sc)
sc$barcode <- gsub(sc$barcode, pattern="_1", replacement= "-0")

sc <- RenameCells(sc, new.names = sc$barcode)
rownames(sc@meta.data) <- sc$barcode
}

sc$barcode <- colnames(sc)
sc$umap_1 <- sc@reductions$umap@cell.embeddings[ , 1]
sc$umap_2 <- sc@reductions$umap@cell.embeddings[ , 2]

write.table(sc@meta.data, file = paste(outdir, "/metadata.csv", sep = ""), quote = F, row.names = F, sep = ",")

# save expression counts matrix
counts <- GetAssayData(sc, assay = "RNA", slot = "counts")
writeMM(counts, file = paste0(outdir, "/counts.mtx"))

# save dimesnionality reduction matrix
write.csv(sc@reductions$pca@cell.embeddings, file = paste0(outdir, "/pca.csv"), quote = F, row.names = F)

# save gene names
write.table(data.frame(gene = rownames(counts)), file = paste0(outdir, "/gene_names.csv"), quote = F, row.names = F, col.names = F)
