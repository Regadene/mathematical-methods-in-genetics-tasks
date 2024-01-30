library(DESeq2)
library(ggvenn)
library(apeglm)
library(EnhancedVolcano)
library(org.Sc.sgd.db)
library(clusterProfiler)
library(fpc)

#1
matrix <- read.table("merged_counts_matrix.tsv", head = T, sep = "\t")

Counts <- as.matrix(matrix[ , -1])
colnames(Counts) <- names(matrix)[-1]
rownames(Counts) <- matrix$gene_id

#2
Sums <- apply(Counts, 1, function(x) {
  sum(x>0)
})
table(Sums == 6)
Counts_info <- Counts[Sums == 6, ]

#построение графика нормализованных данных
Means <- rowMeans(Counts_info)
Sds <- apply(Counts_info, 1, sd)
plot(x = Means, y = Sds, col = rgb(0,0,0, alpha = 0.2) )

#построение графика нормализованных данных с использованием DESeq2
Counts.norm <- rlog(Counts_info)
Means.rlog <- rowMeans(Counts.norm)
Sds.rlog <- apply(Counts.norm, 1, sd)
plot(x = Means.rlog, y = Sds.rlog, col = rgb(0,0,0, alpha = 0.2))

#3
Dist <- dist(t(Counts.norm), method = 'euclidean')
Clust <- hclust(Dist, method = 'ward.D2')
plot(Clust)

#подбор оптимального количества кластеров
NClust <- 2:5
CH_hc_calculation <- function(k) {
  calinhara (x = t(Counts.norm), clustering = cutree(tree = Clust, k = k))
}
CH_hc = sapply(X = NClust, FUN = CH_hc_calculation)
plot(x = NClust, y = CH_hc, type = 'b')

#5
Counts_info_no_del_0 <- Counts_info[ , -1]
Means_no_del_0 <- rowMeans(Counts_info_no_del_0)
Sds_no_del_0 <- apply(Counts_info_no_del_0, 1, sd)

Counts.norm_no_del_0 <- rlog(Counts_info_no_del_0)
Means.rlog_no_del_0 <- rowMeans(Counts.norm_no_del_0)
Sds.rlog_no_del_0 <- apply(Counts.norm_no_del_0, 1, sd)
Dist_no_del_0 <- dist(t(Counts.norm_no_del_0), method = 'euclidean')
Clust_no_del_0 <- hclust(Dist_no_del_0, method = 'ward.D2')

NClust_no_del_0 <- 2:5
CH_hc_calculation_no_del_0 <- function(k) {
  calinhara (x = t(Counts.norm_no_del_0), clustering = cutree(tree = Clust_no_del_0, k = k))
}
CH_hc_no_del_0 = sapply(X = NClust_no_del_0, FUN = CH_hc_calculation_no_del_0)
plot(x = NClust_no_del_0, y = CH_hc_no_del_0, type = 'b')

#6
str(Counts_info_no_del_0, width = 40, strict.width = 'cut')
coldata <- data.frame(genotype = factor(c(rep('del',2), rep('wt', 3)), levels = c('wt', 'del')))
rownames(coldata) <- colnames(Counts_info_no_del_0)
str(coldata)
DataSet <- DESeqDataSetFromMatrix(countData = Counts_info_no_del_0, colData = coldata, design = ~genotype)
for_pca <- rlog(DataSet)
plotPCA(object = for_pca, intgroup = 'genotype', ntop = 5000, returnData = F)

#7
DataSetAnalysis <- DESeq(DataSet)
DEresults <- results(DataSetAnalysis, alpha = 0.05, lfcThreshold = 2, contrast = c('genotype', 'del', 'wt'))
DEresults_noNA <- DEresults[complete.cases(DEresults), ]
down <- rownames(DEresults_noNA[DEresults_noNA$padj < 0.05 & DEresults_noNA$log2FoldChange -1, ])
up <- rownames(DEresults_noNA[DEresults_noNA$padj < 0.05 & DEresults_noNA$log2FoldChange -1, ])
ggvenn(data = list(down = down, up = up, all = rownames(DEresults)), text_size = 3, set_name_size = 3)

#8
rownames (DEresults[DEresults$log2FoldChange < -2 & DEresults$padj < 0.05, ])
#это ген ULP2

#9
plotCounts(DataSetAnalysis, gene = 'cds-NP_012233.1', intgroup = 'genotype')
Counts_info_no_del_0[rownames(Counts_info_no_del_0) == 'cds-NP_012233.1', -1]
plotMA(DEresults)
DEresults_LFS_norm <- lfcShrink(DataSetAnalysis, coef = 'genotype_del_vs_wt', type = 'norm', res = DEresults, lfcThreshold = 1)
DEresults_LFS_apeglm <- lfcShrink(DataSetAnalysis, coef = 'genotype_del_vs_wt', type = 'apeglm', res = DEresults, lfcThreshold = 1)
plotMA(DEresults)
summary(DEresults_LFS_apeglm)
summary(DEresults_LFS_norm)

#10
columns(org.Sc.sgd.db)[c(1, 2, 4, 7, 11, 12, 19, 20, 21, 22, 24)]
p <- "cds-(NP_[0-9]*).[0-9]"
rownames(DEresults_LFS_apeglm) <- sub(pattern = p, replacement = "\\1", x = rownames(DEresults_LFS_apeglm), perl = T)
rownames(DEresults_LFS_apeglm)[1:2]
DEresults_LFS_apeglm$alias <- mapIds(org.Sc.sgd.db, key = rownames(DEresults_LFS_apeglm), column = "COMMON", keytype = "REFSEQ", multiVals = "first")
DEresults_LFS_apeglm$alias[1:2]
tmp <- DEresults_LFS_apeglm
EnhancedVolcano(tmp, lab = tmp$alias, x = "log2FoldChange", y = "svalue", pCutoff = 0.01, FCcutoff = 2, titleLabSize = 10, subtitleLabSize = 7, captionLabSize = 8, axisLabSize = 8, legendLabSize = 8, legendIconSize = 1, legendPosition = "none", pointSize = 1, labSize = 1.5)

#11
df_DERes <- as.data.frame(DEresults_LFS_apeglm)
df_DERes_down <- filter(df_DERes, svalue <= 0.05, log2FoldChange < -1)
genes_down <- rownames(df_DERes_down)

ego_down_BP = enrichGO(gene = genes_down, keyType = "REFSEQ", OrgDb = org.Sc.sgd.db, universe = rownames(as.data.frame(DEresults_LFS_apeglm)), ont = "BP")
clusterProfiler::dotplot(ego_down_BP, showCategory = 10)
#у нас с Шохрухом разные графики, так как в genes_down лежат разные данные. Не знаю, у кого ошибка(

