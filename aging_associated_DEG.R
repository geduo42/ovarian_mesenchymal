library(Seurat)
library(dplyr)
library(magrittr)
library(readr)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)

args <- commandArgs()

outdir 

# load data
load(paste(path, "merged_sub.RData", sep = ""))
sc <- merged_sub

database <- "/home/soft/public/GenAge/genage_models.csv"
genage <- read.csv(database)
colnames(genage)[2] <- "genes"

#get ageing associated degs
find_ageing_associated_degs <- function(cluster, p = 0.01, padj = 0.05){
    markers <- FindMarkers(sc, ident.1 = "aged", ident.2 = "young", group.by = 'orig.ident', subset.ident = cluster, logfc.threshaged = 0.25, slot = "counts", test.use= "t",  min.pct = 0.1)
    markers$cluster <- cluster
    markers_filtered <- markers[markers$p_val < p & markers$p_val_adj < padj,  ]
    markers_filtered <- cbind(genes = rownames(markers_filtered), markers_filtered)
    return(markers_filtered)
}

# get proportion of each gene expressed in each clusters
get_expression_proportion <- function(id){
    subdata <- subset(sc, idents = id)
    counts <- GetAssayData(subdata, assay = "RNA", slot = "counts")
    exp_proportion <- apply(counts, 1, function(x)(sum(x > 0))/ ncol(counts))
    return(exp_proportion)
}

#get proportion of each gene expressed in each clusters
cluster <- levels(sc)
expr_propor <- sapply(cluster, function(x)get_expression_proportion(x))
colnames(expr_propor) <- cluster

# get ageing associated degs
ageing_associated_degs <- lapply(cluster, function(x) find_ageing_associated_degs(x))
names(ageing_associated_degs) <- cluster

ageing_associated_down <- lapply(cluster, function(x) ageing_associated_degs[[x]][ageing_associated_degs[[x]]$avg_log2FC < 0, ])
names(ageing_associated_down) <- cluster

gene_down <- do.call(rbind, ageing_associated_down)
write.table(gene_down, file = paste(outdir, "gene_downregulated.txt", sep = "/"), sep = "\t", quote = F)

ageing_associated_up <- lapply(cluster, function(x) ageing_associated_degs[[x]][ageing_associated_degs[[x]]$avg_log2FC > 0, ])
names(ageing_associated_up) <- cluster
gene_up <- do.call(rbind, ageing_associated_up)
write.table(gene_up, file = paste(outdir, "gene_upregulated.txt", sep = "/"), sep = "\t", quote = F)


# get num of up or down regulated genes in each cluster
dim(gene_up)
table(gene_up$cluster)
dim(gene_down)
table(gene_down$cluster)


# get num of cluster specific up or down regulated genes
dup_up <- unique(gene_up$genes[duplicated(gene_up$genes)])
uni_up <- gene_up[!gene_up$genes %in% dup_up,]

dup_down <- unique(gene_down$genes[duplicated(gene_down$genes)])
uni_down <- gene_down[!gene_down$genes %in% dup_down,]

dim(uni_up)
table(uni_up$cluster)
dim(uni_down)
table(uni_down$cluster)



genage_dotplot <- function(ageing_degs, color, change){
#get ageing associated genes in geneage
    ageing_genage <- lapply(cluster, function(x) inner_join(ageing_degs[[x]], genage, by = "genes"))
    ageing_genage <- do.call(rbind, ageing_genage)
    ageing_genage <-  ageing_genage[ageing_genage$genes != "Gsta4", ]
    ageing_genage_long <- ageing_genage[rep(1 : dim(ageing_genage)[1], length(cluster)), ]

#save data
    fname <- paste(substitute(ageing_degs), "txt", sep = ".")
    write.table(ageing_genage, file = paste(outdir, fname, sep = "/"), sep = "\t", quote = F)

#express proportion of ageing_gnage
    expr_gene <- expr_propor[ageing_genage$genes, ]
    expr_gene_long <- melt(expr_gene, value.name = "expr_proportion", varnames= c("gene", "celltype"))

# data for plot
    df <- cbind(ageing_genage_long, expr_gene_long)
    tmp1 <- paste(ageing_genage$genes, ageing_genage$cluster, sep = "_")
    tmp2 <- paste(df$gene, df$celltype, sep = "_")
    df$change <- tmp2 %in% tmp1
    df$celltype <- factor(df$celltype, levels = levels(sc), ordered  = T)
	df$expr_proportion = df$expr_proportion * 100
    pname <- substitute(ageing_degs)
    p <- ggplot(df) + 
		geom_point(mapping = aes(x = celltype, y = gene, size = expr_proportion, color = change)) + 
		labs(title = pname, size = "Cells expressing the \n indicated genes(%)") + 
		theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
		scale_color_manual(values = c("darkgray", color), breaks = c("FALSE", "TRUE"), labels = c("Unchanged", change)) + 
		scale_size_continuous(breaks = c(20, 40, 80), label = c(20, 40, 80), range = c(2, 8)) +
		theme_bw() + 
		guides(color = guide_legend(override.aes = list(size = 6))) + 
		theme(legend.text=element_text(size = 12), legend.title=element_text(size = 12))
    return(p)
}


pdf(paste(outdir, "ageing_genage_genes.pdf", sep = "/"))
genage_dotplot(ageing_associated_degs, color = "green", change = "Changed")
genage_dotplot(ageing_associated_down, color = "deepskyblue", change = "Downregulated")
genage_dotplot(ageing_associated_up, color = "red", change = "Upregulated")

dup_gene_up <- gene_up$genes[duplicated(gene_up$genes)]
dup_gene_down <- gene_down$genes[duplicated(gene_down$genes)]

uni_up <- gene_up[!(gene_up$genes %in% dup_gene_up), ]
uni_down <- gene_down[!(gene_down$genes %in% dup_gene_down), ]

write.table(uni_up, file = paste(outdir, "ageing_associated_up_cluster_specific.txt", sep = "/"), sep = "\t", quote = F)
write.table(uni_down, file = paste(outdir, "ageing_associated_down_cluster_specific.txt", sep = "/"), sep = "\t", quote = F)

dev.off()


