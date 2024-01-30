#8.1
O <- c(5474, 1850)
E <- c(0.75, 0.25) * sum(O)
chi2 <- sum (((E-O)^2)/E)
chi2

#8.2
chisq.test(x = O, p = c(0.75, 0.25))

#8.3
Mendel <- read.table('Mendel_data.tsv', head = F, sep = '\t')
apply(X = Mendel[, 2:3], MARGIN = 1, FUN = function(x) {
  chisq.test(x = x, p = c(0.75, 0.25))$statistic
})
apply(X = Mendel[, 2:3], MARGIN = 1, FUN = function(x) {
  chisq.test(x = x, p = c(0.75, 0.25))$p.value
})

library(pwr)
#8.4
pwr.chisq.test(w = 0.3961054, df = 1, sig.level = 0.05, power = 0.9)

library(DescTools)
#8.5
Sup35_matrix <- matrix(c(417, 8548, 76, 11482), 2, 2, byrow = TRUE)
chisq.test(Sup35_matrix)
CramerV(x = Sup35_matrix, conf.level = 0.95)

library(ggplot2)
#8.6
alleles_P_S_df <- data.frame(
  Allele = c("sup35", "SUP35"),
  P = c(0.66, 4.9),
  S = c(0.08, 0.23)
)

ggplot(data = alleles_P_S_df, mapping = aes(x = Allele, y = P)) + 
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = P - S, ymax = P + S), width = 0.2, position = position_dodge(0.9)) +
  annotate(geom = "text", x = 1.5, y = 5.5, label = "***", size = 6) +
  theme_bw()

#8.7.1
Prion_induction_data <- read.table('Prion_induction.tsv', head = T, sep = '\t')
head(Prion_induction_data)
table(Prion_induction_data$allele)

#8.7.2 use Fisher
fisher.test(x = as.matrix(Prion_induction_data[Prion_induction_data$allele == "control", -1]))
#вывод: нет разницы между повторностями, значит, можем их суммировать

#8.7.3
fisher.test(x = as.matrix(Prion_induction_data[Prion_induction_data$allele == "M3", -1]), simulate.p.value = T)$p.value
fisher.test(x = as.matrix(Prion_induction_data[Prion_induction_data$allele == "M5", -1]), simulate.p.value = T)$p.value
fisher.test(x = as.matrix(Prion_induction_data[Prion_induction_data$allele == "WT", -1]), simulate.p.value = T)$p.value

#8.8
phyper(q = 28, m = 2613, n = 15310, k = 57, lower.tail = F)
#Вероятность получить такой результат - очень низкая = 1.10373e-10, поэтому можно говорить о том, что среди генов с повышенной экспрессией неслучайно много генов, ассоциированных с интересующим клеточным процессом

