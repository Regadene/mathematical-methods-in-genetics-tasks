#1
Gulob_main_array <- read.table('Gulob_main.tsv', head = T, sep = '\t')
Gulob_main_array[, -1][Gulob_main_array[, -1] < 100] <- NA
Gulob_main_array[, -1][Gulob_main_array[, -1] > 16000] <- NA
Gulob_main_array <- Gulob_main_array[complete.cases(Gulob_main_array), ]

#2
sd_Gulob <- apply(Gulob_main_array[ , -1], 1, sd)
mean_Gulob <- apply(Gulob_main_array[ , -1], 1, mean)
plot(mean_Gulob, sd_Gulob)

#3
matrix_Gulob <- as.matrix(Gulob_main_array[ , -1])
row.names(matrix_Gulob) <- Gulob_main_array$Gene.Accession.Number
colnames(matrix_Gulob) <- c(rep('ALL', 27), rep('AML', 11))

#4
library(vegan)
data.pca <- rda(X = t(matrix_Gulob), scale = T)

#5
summary(data.pca)
#0.04501

#6
# ген U22376_cds2_s_at коэффициент = -0.6365

#7 для первого
tmp <- scores(x = data.pca, choices = 1:3, tidy = T)
ggplot(tmp[tmp$score == 'sites', ], aes(x = PC1, y = PC2, color = label)) +
  geom_point() + 
  theme_bw()

ggplot(tmp[tmp$score == "species", ], aes(x = 0, y = 0, xend = PC1, yend = PC2, label = label)) + 
  geom_segment(arrow = arrow(length = unit(0.2, "cm")), alpha = 0.5, color = "red") + 
  geom_text(aes(x = PC1, y = PC2), hjust = 1, vjust = 1, size = 2) + 
  theme_bw()
