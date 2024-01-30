#2.1
qrt <- read.table("sample_groups.txt", head = T)
qrt_data <- read.table("QRT_PCR_LR.txt", head = T)
qrt_grouped <- merge(qrt, qrt_data, by.x = "tissue", by.y = "Tissue")

#2.2
qrt_grouped_delta <- qrt_grouped
qrt_grouped_delta[, -c(1, 2)] <- qrt_grouped[, -c(1, 2)] - qrt_grouped$UBC

#2.3

#Смотрим на значение shapiro.test
#если меньше 0,05 то выбираем wilcox.test
#если больше то t.test
shapiro.test(qrt_grouped_delta$AKR1C1[qrt_grouped_delta$group == 1])$p.value
shapiro.test(qrt_grouped_delta$AKR1C1[qrt_grouped_delta$group == 2])$p.value
wilcox.test(formula = AKR1C1 ~ group,data=qrt_grouped_delta)$p.value

shapiro.test(qrt_grouped_delta$SDHA[qrt_grouped_delta$group == 1])$p.value
shapiro.test(qrt_grouped_delta$SDHA[qrt_grouped_delta$group == 2])$p.value
wilcox.test(formula = SDHA ~ group,data = qrt_grouped_delta)$p.value

shapiro.test(qrt_grouped_delta$CD44[qrt_grouped_delta$group == 1])$p.value
shapiro.test(qrt_grouped_delta$CD44[qrt_grouped_delta$group == 2])$p.value
var.test(formula=CD44 ~ group,data = qrt_grouped_delta)$p.value
t.test(formula=CD44 ~ group,data = qrt_grouped_delta)$p.value

#2.4
boxplot(formula = AKR1C1 ~ group, data = qrt_grouped_delta[, c("group", "AKR1C1")], outline = F)
stripchart(x = AKR1C1 ~ group, data = qrt_grouped_delta[, c("group", "AKR1C1")], add = T, vertical = T, pch = 19, method = "jitter", jitter = 0.3, cex = 0.2)
text(x = 1.5, y = 12, labels = "***", cex = 1.5)


boxplot(formula = SDHA ~ group, data = qrt_grouped_delta[, c("group", "SDHA")], outline = F)
stripchart(x = SDHA ~ group, data = qrt_grouped_delta[, c("group", "SDHA")], add = T, vertical = T, pch = 19, method = "jitter", jitter = 0.3, cex = 0.2)
text(x = 1.5, y = 12, labels = "***", cex = 1.5)


boxplot(formula = CD44 ~ group, data = qrt_grouped_delta[, c("group", "CD44")], outline = F)
stripchart(x = CD44 ~ group, data = qrt_grouped_delta[, c("group", "CD44")], add = T, vertical = T, pch = 19, method = "jitter", jitter = 0.3, cex = 0.2)
text(x = 1.5, y = 12, labels = "***", cex = 1.5)


library(pwr)
#2.5
pwr.t.test(n = NULL, d = 1/2, sig.level = 0.05, power = 0.99)

#2.6
pwr.t2n.test(n1 = 100, n2 = NULL, d = 1/2, sig.level = 0.05, power = 0.99)

#2.7
Ds <- seq(0.2, 1, 0.1)
Ns <- sapply(Ds, FUN = function(x) {
  pwr.t.test(n = NULL, d = x, sig.level = 0.05, power = 0.99)$n
})
plot(x = Ds, y = Ns, type = "p", pch = 19)

#2.8
mean1 <- mean(qrt_grouped_delta$CD44[qrt_grouped_delta$group == 1], na.rm = T)
mean2 <- mean(qrt_grouped_delta$CD44[qrt_grouped_delta$group == 2], na.rm = T)
sigma <- var(qrt_grouped_delta$CD44, na.rm = T)
pwr.t.test(sig.level = 0.05, power = 0.9, d = (mean1 - mean2)/sigma)

library(tidyr)
#2.9
qrt_long <- pivot_longer(data = qrt_grouped_delta, cols = -c("tissue", "group"))
wilcox_res <- qrt_long %>%
  group_by(name) %>%
  summarise(p.value = wilcox.test(value ~ group)$p.value)
wilcox_res$p.adj.holm <- p.adjust(p = wilcox_res$p.value, method = "holm")
wilcox_res$p.adj.BH <- p.adjust(p = wilcox_res$p.value, method = "BH")

length(wilcox_res$name[wilcox_res$p.value < 0.05])
length(wilcox_res$name[wilcox_res$p.adj.holm < 0.05])
length(wilcox_res$name[wilcox_res$p.adj.BH < 0.05])

#2.10
qrt_long <- pivot_longer(data = qrt_grouped_delta, cols = -c("tissue", "group"))
wilcox_res <- qrt_long %>%
  group_by(name) %>%
  summarise(p.value = wilcox.test(value ~ group)$p.value, FC = 2^(-1 * (median(value[group == 2]) - median(value[group == 1]))))

wilcox_res$LFC <- log2(wilcox_res$FC)
wilcox_res$p.adj.holm <- p.adjust(p = wilcox_res$p.value, method = "holm")

wilcox_res <- wilcox_res[complete.cases(wilcox_res), ]
sum(wilcox_res$LFC[wilcox_res$p.adj.holm < 0.05] > 1)

library(ggplot2)
#2.11
ggplot(data = wilcox_res, mapping = aes(x = name, y = LFC, fill = p.adj.holm < 0.05)) +
  geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 5))
