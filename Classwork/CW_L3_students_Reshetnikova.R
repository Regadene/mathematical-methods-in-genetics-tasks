#3.1
plasmodium_df <- read.table(file='Plasmodium_rna_prot.tsv',sep = "\t", header = TRUE)
plasmodium_df[plasmodium_df<=0] <- NA

#3.2
shapiro.test_p.value_fun <- function(x) { shapiro.test(x)$p.value }
apply(X = plasmodium_df[ , -1], MARGIN = 2, FUN = shapiro.test_p.value_fun)

#3.3
Plasmodium_rna_prot_log <- plasmodium_df
Plasmodium_rna_prot_log[ , -1] <- log(Plasmodium_rna_prot_log[ , -1])
apply(X = Plasmodium_rna_prot_log[ , -1],MARGIN = 2, FUN = shapiro.test_p.value_fun)


library(pwr)
#3.4
pwr.r.test(r = 0.2, sig.level = 0.05,power = 0.9, alternative = "t")

library(ggplot2)
library(tidyr)
library(dplyr)
#3.5
cor_matrix <- as.data.frame(cor(x = Plasmodium_rna_prot_log[ , 3:5], method = "pearson", use = "pairwise.complete.obs"))
cor_matrix$sample1 <- row.names(cor_matrix)
pivot_longer(data = cor_matrix, cols = -sample1, names_to = "sample2") %>%
  ggplot(mapping = aes(x = sample1, y = sample2, fill = value)) + geom_tile() + theme_bw() + labs(x = "",y = "")

#3.6
plot(formula = RNA_Troph ~ RNA_Schiz, data = Plasmodium_rna_prot_log)

plot(formula = RNA_Ring ~ RNA_Schiz, data = Plasmodium_rna_prot_log)

plot(formula = RNA_Ring ~ RNA_Troph, data = Plasmodium_rna_prot_log)

#3.7
cor.test(x = Plasmodium_rna_prot_log$RNA_Troph, y = Plasmodium_rna_prot_log$RNA_Schiz, method = "pearson")
cor.test(x = Plasmodium_rna_prot_log$RNA_Ring, y = Plasmodium_rna_prot_log$RNA_Schiz, method = "pearson")
cor.test(x = Plasmodium_rna_prot_log$RNA_Ring, y = Plasmodium_rna_prot_log$RNA_Troph, method = "pearson")

#3.8
plot(formula = RNA_Troph ~ RNA_Schiz, data = Plasmodium_rna_prot_log)
legend(x="topleft", legend = "0,765",bty = "n")

plot(formula = RNA_Ring ~ RNA_Schiz, data = Plasmodium_rna_prot_log)
legend(x="topleft", legend = "0,610",bty = "n")

plot(formula = RNA_Ring ~ RNA_Troph, data =Plasmodium_rna_prot_log)
legend(x="topleft", legend = "0,812",bty = "n")

#3.9
cor(x = Plasmodium_rna_prot_log[ , c("RNA_Ring")], y = Plasmodium_rna_prot_log[ , c("Prot_Ring")], method = "kendall", use = "pairwise.complete.obs")
cor(x = Plasmodium_rna_prot_log[ , c("RNA_Mero")], y = Plasmodium_rna_prot_log[ , c("Prot_Mero")], method = "kendall", use = "pairwise.complete.obs")
cor.test(x =Plasmodium_rna_prot_log[ , c("RNA_Ring")], y =Plasmodium_rna_prot_log[ , c("Prot_Ring")], method = "kendall")
cor.test(x =Plasmodium_rna_prot_log[ , c("RNA_Mero")], y = Plasmodium_rna_prot_log[ , c("Prot_Mero")], method = "kendall")

#3.10
#Корреляция между размером мозга и интеллектом
#Корреляция между наличием домашних животных и здоровьем сердца
#Корреляция между ростом растений и числом клеток в листьях
#Корреляция между наличием бороды у мужчин и уровнем тестостерона