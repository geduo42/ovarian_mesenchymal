library(Seurat)
library(monocle)
library(GGally)
library(colorRamps)
library(RColorBrewer)
library(reshape2)
library(clusterProfiler)
library(dplyr)

load(paste(indir, "young_sub.RData", sep = "/"))
sc <- young_sub

# prepare data for monocle
data <- GetAssayData(sc, assay = "RNA", slot = "data")
pd = new("AnnotatedDataFrame", data = sc@meta.data)
fdata = data.frame(gene_short_name = rownames(sc), row.names = rownames(sc))
fd = new("AnnotatedDataFrame", data = fdata)

# creat a monocle object
sc_monocle <- newCellDataSet(data, phenoData = pd, feature = fd, lowerDetectionLimit = 0, expressionFamily = negbinomial.size())

# get SizeFactor and Dispersions
sc_monocle <- estimateSizeFactors(sc_monocle)
sc_monocle <- estimateDispersions(sc_monocle)

# select genes which define the psudo time
variable_genes <- VariableFeatures(sc)
diff <- differentialGeneTest(sc_monocle[variable_genes, ], fullModelFormulaStr = "~seurat_clusters")
deg <- subset(diff, qval < 0.01)
deg <- deg[order(deg$qval, decreasing = F), ] 
deg_genes <- deg$gene_short_name

ordering_genes <- deg_genes
sc_monocle <- setOrderingFilter(sc_monocle, ordering_genes)
plot_ordering_genes(sc_monocle)

# dim reduction
sc_monocle <- reduceDimension(sc_monocle, max_components = 2, method = "DDRTree")

# compute psudotime value
sc_monocle <- orderCells(sc_monocle)

sc_monocle <- orderCells(sc_monocle, root_state = 3)
plot_cell_trajectory(sc_monocle, cell_size = 1, color_by = "State") +
    scale_color_manual(breaks = c(1, 2, 3), values = c("Pink4", "royalblue2", "red3"))

plot_cell_trajectory(sc_monocle, cell_size = 1, color_by = "Pseudotime")

plot_cell_trajectory(sc_monocle, cell_size = 1, color_by = "cell_type") +
    scale_color_manual(breaks = c("oMYB", "Wnt4+", "Col15a1+", "ETCLC", "ITC", "oMPs"), values = c("seashell4", "chocolate4", "peachpuff3", "lightskyblue", "dodgerblue", "hotpink"))

temp <- data.frame(x = 1, y = 1)
blank_fig <- ggplot(temp) + geom_blank()
blank_fig

BEAM_res <- BEAM(sc_monocle, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
BEAM_res <- BEAM_res[ , c("gene_short_name", "pval", "qval")]

hmcols      <- colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(65)
plot_branch <- plot_genes_branched_heatmap(sc_monocle[row.names(subset(BEAM_res, qval < 1e-4)), ],
                                         branch_point = 1,
                                         num_clusters = 3,
                                         cores = 1,
                                         hmcols = hmcols,
                                         use_gene_short_name = F,
                                         branch_colors = c("#bebebe","#009E73", "indianred2"),
                                         show_rownames = F,return_heatmap = T)

blank_fig

plot_genes_keep <- plot_genes_branched_heatmap(sc_monocle[genes_keep, ],
                                         branch_point = 1,
                                         num_clusters = 3,
                                         cores = 1,
                                         hmcols = hmcols,
                                         use_gene_short_name = F,
                                         branch_colors = c("#bebebe","#009E73", "indianred2"),
                                         show_rownames = F,return_heatmap = T)

dev.off()

save(sc_monocle, BEAM_res, plot_genes_keep, file = paste(outdir, "young_monocle.RData"))
