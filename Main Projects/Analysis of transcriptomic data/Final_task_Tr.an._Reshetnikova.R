library(DESeq2)
library(ggvenn)
library(apeglm)
library(EnhancedVolcano)
library(org.Sc.sgd.db)
library(clusterProfiler)
library(fpc)
#Скачиваю с помощью Filezilla на свой компьютер в отдельную папку для работу файлы htseq_counts
# setwd("C:/Users/Rimalon/Desktop/FinalTaskR")

#1
#Загрузка таблиц:
SRR36_htseq <- read.table("SRR8126636_htseq_counts.txt", head = F, sep = "\t")
colnames(SRR36_htseq)[colnames(SRR36_htseq) == "V2"] ="SRR36_wt"

SRR38_htseq <- read.table("SRR8126638_htseq_counts.txt", head = F, sep = "\t")
colnames(SRR38_htseq)[colnames(SRR38_htseq) == "V2"] ="SRR38_wt"

SRR44_htseq <- read.table("SRR8126644_htseq_counts.txt", head = F, sep = "\t")
colnames(SRR44_htseq)[colnames(SRR44_htseq) == "V2"] ="SRR44_wt"

SRR59_htseq <- read.table("SRR8126659_htseq_counts.txt", head = F, sep = "\t")
colnames(SRR59_htseq)[colnames(SRR59_htseq) == "V2"] ="SRR59_del"

SRR66_htseq <- read.table("SRR8126666_htseq_counts.txt", head = F, sep = "\t")
colnames(SRR66_htseq)[colnames(SRR66_htseq) == "V2"] ="SRR66_del"

SRR71_htseq <- read.table("SRR8126671_htseq_counts.txt", head = F, sep = "\t")
colnames(SRR71_htseq)[colnames(SRR71_htseq) == "V2"] ="SRR71_del"

# Объединение таблиц:
merged_36_and_38 <- merge(x = SRR36_htseq , y = SRR38_htseq, by = "V1")
merged_36_and_38_and_44 <- merge(x = merged_36_and_38 , y = SRR44_htseq, by = "V1")
merged_36_and_38_and_44_and_59 <- merge(x = merged_36_and_38_and_44 , y = SRR59_htseq , by = "V1")
merged_36_and_38_and_44_and_59_and_66 <- merge(x = merged_36_and_38_and_44_and_59 , y = SRR66_htseq , by = "V1")
merged_all <- merge(x = merged_36_and_38_and_44_and_59_and_66 , y = SRR71_htseq , by = "V1")

# Сохранение объединенной таблицы в файл
write.table(merged_all, file = "counts_Reshetnikova.csv")

# Загрузка данных из таблицы
merged_all_counts <- read.table("counts_Reshetnikova.csv", head = T,  sep = " ")
merged_all_counts <- tail(merged_all_counts, -5)

# Подготовка матрицы:
as_matrix_merged_all <- as.matrix(merged_all_counts [, -1])
row.names(as_matrix_merged_all) <- merged_all_counts$V1

# Проверка на отсутствие прочтений
sums_merged_all <- apply(as_matrix_merged_all, 1, function(x) {
  sum(x > 0)
})
as_matrix_merged_all_info <- as_matrix_merged_all[sums_merged_all == 6, ]


# построение графика нормализованных данных с использованием DESeq2
as_matrix_merged_all_info.norm <- rlog(as_matrix_merged_all_info)
Means.rlog <- rowMeans(as_matrix_merged_all_info.norm)
Sds.rlog <- apply(as_matrix_merged_all_info.norm, 1, sd)
plot(x = Means.rlog, y = Sds.rlog, col = rgb(0,0,0, alpha = 0.2))

#2 Кластерный анализ 

# построение кластерной дендрограммы
Dist <- dist(t(as_matrix_merged_all_info.norm), method = 'euclidean')
Clust <- hclust(Dist, method = 'ward.D2')
plot(Clust)

# подбор оптимального количества кластеров
NClust <- 2:5
CH_hc_calculation <- function(k) {
  calinhara (x = t(as_matrix_merged_all_info.norm), clustering = cutree(tree = Clust, k = k))
}
CH_hc = sapply(X = NClust, FUN = CH_hc_calculation)
plot(x = NClust, y = CH_hc, type = 'b')
# По графику видим, что оптимальное N кластеров = 3. Значение SRR59_del выбивается. Удаляем выбивающуюся пробу (у меня SRR59 в 4-ом столбце)
# Дальше будем работать с данными без выбивающейся пробы:

as_matrix_merged_all_info_no_SRR_59 <- as_matrix_merged_all_info[ , -4]
Means_no_SRR_59 <- rowMeans(as_matrix_merged_all_info_no_SRR_59)
Sds_no_SRR_59 <- apply(as_matrix_merged_all_info_no_SRR_59, 1, sd)

#Подбираем оптимальное количество кластеров без выбивающейся пробы:
as_matrix_merged_all.norm_no_SRR_59 <- rlog(as_matrix_merged_all_info_no_SRR_59)
Means.rlog_no_SRR_59 <- rowMeans(as_matrix_merged_all.norm_no_SRR_59)
Sds.rlog_no_SRR_59 <- apply(as_matrix_merged_all.norm_no_SRR_59, 1, sd)
Dist_no_SRR_59 <- dist(t(as_matrix_merged_all.norm_no_SRR_59), method = 'euclidean')
Clust_no_SRR_59 <- hclust(Dist_no_SRR_59, method = 'ward.D2')

NClust_no_SRR_59 <- 2:5
CH_hc_calculation_no_SRR_59 <- function(k) {
  calinhara (x = t(as_matrix_merged_all.norm_no_SRR_59), clustering = cutree(tree = Clust_no_SRR_59, k = k))
}
CH_hc_no_SRR_59 = sapply(X = NClust_no_SRR_59, FUN = CH_hc_calculation_no_SRR_59)
plot(x = NClust_no_SRR_59, y = CH_hc_no_SRR_59, type = 'b')
# Видим, что теперь оптимальное количество кластеров = 2. 

# Создаём график РСА:
str(as_matrix_merged_all_info, width = 40, strict.width = 'cut')
coldata <- data.frame(genotype = factor(c(rep('wt',3), rep('del', 3)), levels = c('wt', 'del')))
rownames(coldata) <- colnames(as_matrix_merged_all_info)
str(coldata)
DataSet <- DESeqDataSetFromMatrix(countData = as_matrix_merged_all_info, colData = coldata, design = ~genotype)
for_pca <- rlog(DataSet)
plotPCA(object = for_pca, intgroup = 'genotype', ntop = 5000, returnData = F)
#График РСА как и кластерный анализ показывает выбивающуюся пробу. 

#3 Анализ дифф. экспрессии

DataSetAnalysis <- DESeq(DataSet)
DEresults <- results(DataSetAnalysis, alpha = 0.05, lfcThreshold = 1, contrast = c('genotype', 'del', 'wt'))
summary(DEresults)
# команда summary для DEResults показывает нули для LFC>1 и LFC<1. Используем поправку apeglm для lFC:

DEresults_LFS_apeglm <- lfcShrink(DataSetAnalysis, coef = 'genotype_del_vs_wt', type = 'apeglm', res = DEresults, lfcThreshold = 1)
summary(DEresults_LFS_apeglm)
# теперь при вызове summary для оптимизированного DEresults 
# имеем 5 генов для которых LFC>1 (up), и 4 гена для которых LFC<1 (down)

# после поправки полученные данные имеют другой формат, в таблице нет столбца padj,
# так как при apeglm используется svalue 
# (который имеет повышенную точность, как я поняла из документации).
# Уровень значимости для svalue при поправке apeglm был 0.005, поэтому 
# для получения кодификаторов искомых генов использую svalue < 0.005
DEresults_LFS_apeglm_noNA <- DEresults_LFS_apeglm[complete.cases(DEresults), ]
down <- rownames(DEresults_LFS_apeglm_noNA[DEresults_LFS_apeglm_noNA$svalue < 0.005 & DEresults_LFS_apeglm_noNA$log2FoldChange < -1, ])
up <- rownames(DEresults_LFS_apeglm_noNA[DEresults_LFS_apeglm_noNA$svalue < 0.005 & DEresults_LFS_apeglm_noNA$log2FoldChange > 1, ])

# построение диаграммы Венна
ggvenn(data = list(down = down, up = up, all = rownames(DEresults)), text_size = 3, set_name_size = 3)

# составление рег. выражения для того, чтобы в дальнейшем обращаться к REFSEQ.
refseq_regexpr <- 'cds-(NP_[0-9]*).[0-9]'
DEresults_LFS_apeglm_refseq <- DEresults_LFS_apeglm
rownames(DEresults_LFS_apeglm_refseq) <- sub(pattern = refseq_regexpr, replacement = "\\1", x = rownames(DEresults_LFS_apeglm_refseq), perl = T)

# построение графика EnchancedVolcano
DEresults_LFS_apeglm_refseq$alias <- mapIds(org.Sc.sgd.db, key = rownames(DEresults_LFS_apeglm_refseq), column = "COMMON", keytype = "REFSEQ", multiVals = "first")
tmp <- DEresults_LFS_apeglm_refseq
EnhancedVolcano(tmp, lab = tmp$alias, x = "log2FoldChange", y = "svalue", pCutoff = 0.005, FCcutoff = 1, titleLabSize = 10, subtitleLabSize = 7, captionLabSize = 8, axisLabSize = 8, legendLabSize = 8, legendIconSize = 1, legendPosition = "none", pointSize = 1, labSize = 1.5)

# Анализ обогащения терминами Gene Ontology:
# экспрессия уменьшена
df_DERes_refseq <- as.data.frame(DEresults_LFS_apeglm_refseq)
df_DERes_down_refseq <- filter(df_DERes_refseq, svalue <= 0.005, log2FoldChange < -1)
genes_down_refseq <- rownames(df_DERes_down_refseq)


ego_down_BP_refseq <- enrichGO(gene = genes_down_refseq, keyType = "REFSEQ", OrgDb = org.Sc.sgd.db, universe = rownames(as.data.frame(DEresults_LFS_apeglm_refseq)), ont = "BP")
ego_down_CC_refseq <- enrichGO(gene = genes_down_refseq, keyType = "REFSEQ", OrgDb = org.Sc.sgd.db, universe = rownames(as.data.frame(DEresults_LFS_apeglm_refseq)), ont = "CC")
ego_down_MF_refseq <- enrichGO(gene = genes_down_refseq, keyType = "REFSEQ", OrgDb = org.Sc.sgd.db, universe = rownames(as.data.frame(DEresults_LFS_apeglm_refseq)), ont = "MF")
clusterProfiler::dotplot(ego_down_MF_refseq, showCategory = 10)
# для 2 из 4 генов были найдены соответствия в базе данных MF, для BP и CC соответствий нет

# экспрессия увеличена
df_DERes_up_refseq <- filter(df_DERes_refseq, svalue <= 0.005, log2FoldChange > 1)
genes_up_refseq <- rownames(df_DERes_up_refseq)

ego_up_BP_refseq <- enrichGO(gene = genes_up_refseq, keyType = "REFSEQ", OrgDb = org.Sc.sgd.db, universe = rownames(as.data.frame(DEresults_LFS_apeglm_refseq)), ont = "BP")
ego_up_CC_refseq <- enrichGO(gene = genes_up_refseq, keyType = "REFSEQ", OrgDb = org.Sc.sgd.db, universe = rownames(as.data.frame(DEresults_LFS_apeglm_refseq)), ont = "CC")
ego_up_MF_refseq <- enrichGO(gene = genes_up_refseq, keyType = "REFSEQ", OrgDb = org.Sc.sgd.db, universe = rownames(as.data.frame(DEresults_LFS_apeglm_refseq)), ont = "MF")
# не найдено соответствий ни для одного из 5 генов


#Доп задание
genes_down_sgd_ids <- mapIds(org.Sc.sgd.db, keys = genes_down_refseq, column = "SGD", keytype = "REFSEQ")


genes_down_orf_ids <- mapIds(org.Sc.sgd.db, keys = genes_down_refseq, column = "ORF", keytype = "REFSEQ")
genes_orf_ids <- as.vector(mapIds(org.Sc.sgd.db, keys = rownames(as.data.frame(DEresults_LFS_apeglm_refseq)), column = "ORF", keytype = "REFSEQ"))
sgdCHRLOC <- org.Sc.sgdCHRLOC
sgdCHRLOC_list <- as.list(sgdCHRLOC[mappedkeys(sgdCHRLOC)])
genes_down_orf_locations <- sgdCHRLOC_list[names(sgdCHRLOC_list) %in% genes_down_orf_ids]
genes_orf_locations <- sgdCHRLOC_list[names(sgdCHRLOC_list) %in% genes_orf_ids]
