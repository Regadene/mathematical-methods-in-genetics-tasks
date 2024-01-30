library(dplyr)
library(tidyr)

#Подготовка данных
qrt_groups <- read.table("sample_groups.txt", head = T)
qrt_data <- read.table("QRT_PCR_LR.txt", head = T)
grt <- merge(qrt_groups, qrt_data, by.x = "tissue", by.y = "Tissue")

grt_delta_UBS_df <- grt
grt_delta_UBS_df[, -c(1, 2)] <- grt[, -c(1, 2)] - grt$UBC

grt_delta_SDHA_df <- grt
grt_delta_SDHA_df[, -c(1, 2)] <- grt[, -c(1, 2)] - grt$SDHA

#Проверка нормальности распределения данных
shapiro.test_p.value_fun <- function(x) { shapiro.test(x)$p.value }

shapiro_testUBS <- apply(X = grt_delta_UBS_df[, -c(1, 2, 64)], MARGIN = 2, FUN = shapiro.test_p.value_fun)

#Проверка: есть ли значения p-value больше уровня значимости
length(shapiro_testUBS[shapiro_testUBS > 0.05]) > 0
#[1] FALSE

shapiro_testSDHA <- apply(X = grt_delta_SDHA_df[, -c(1, 2, 58)], MARGIN = 2, FUN = shapiro.test_p.value_fun)

#Проверка: есть ли значения p-value больше уровня значимости
length(shapiro_testSDHA[shapiro_testSDHA > 0.05]) > 0
#[1] FALSE

#Все значения p-value, меньше уровня значимости => Н0 отвергаем, принимаем Н1. (Данные считаем не распределенными нормально)
#Используем тест Вилкоксона

#Выберем 10 генов:
qrt_ten_genes_names <- names(grt_delta_UBS_df[, c(14:23)])
#Выбранные гены:
#"CLSTN1"  "CPSG3"   "DDC"     "DPYSL3"  "ECEL1"   "ELAVL4"  "EPB41L3" "EPHA5"   "EPN2"    "FYN"

#Тесты Вилкоксона для выбранных 10 генов с нормированием по референсному гену UBS
grt_delta_UBS <- pivot_longer(data = grt_delta_UBS_df, cols = qrt_ten_genes_names)
wilcox.test_p.values_UBS <- grt_delta_UBS %>%
  group_by(name) %>%
  summarise(p.value = wilcox.test(value ~ group)$p.value)

#Поправки на множественные сравнения для тестов с нормированием по референсному гену UBS
wilcox.test_p.values_UBS$p.adj.holm <- p.adjust(p = wilcox.test_p.values_UBS$p.value, method = "holm")
wilcox.test_p.values_UBS$p.adj.BH <- p.adjust(p = wilcox.test_p.values_UBS$p.value, method = "BH")

length(wilcox.test_p.values_UBS$name[wilcox.test_p.values_UBS$p.value < 0.05])
#[1]  10
length(wilcox.test_p.values_UBS$name[wilcox.test_p.values_UBS$p.adj.holm < 0.05])
#[1]  10
length(wilcox.test_p.values_UBS$name[wilcox.test_p.values_UBS$p.adj.BH < 0.05])
#[1]  10

#Все значения p-value, меньше уровня значимости => Н0 отвергаем, принимаем Н1. (Данные считаем имеющими статистически значимые различия)
#В группах значения экспрессий для выбранных 10 генов с нормированием по референсному гену UBS имеют различия
 
#Тесты Вилкоксона для выбранных 10 генов с нормированием dCt по референсному гену SDHA
grt_delta_SDHA <- pivot_longer(data = grt_delta_SDHA_df, cols = all_of(qrt_ten_genes_names))
wilcox.test_p.values_SDHA <- grt_delta_SDHA %>%
  group_by(name) %>%
  summarise(p.value = wilcox.test(value ~ group)$p.value)

#Поправки на множественные сравнения для тестов с нормированием dCt по референсному гену SDHA
wilcox.test_p.values_SDHA$p.adj.holm <- p.adjust(p = wilcox.test_p.values_SDHA$p.value, method = "holm")
wilcox.test_p.values_SDHA$p.adj.BH <- p.adjust(p = wilcox.test_p.values_SDHA$p.value, method = "BH")

#Гены, для которых полагаем различия:
wilcox.test_p.values_SDHA$name[wilcox.test_p.values_SDHA$p.value < 0.05]
#[1] "CLSTN1"  "CPSG3"   "ECEL1"   "EPB41L3" "EPHA5"   "EPN2"
wilcox.test_p.values_SDHA$name[wilcox.test_p.values_SDHA$p.adj.holm < 0.05]
#[1] "CPSG3" "ECEL1" "EPN2"
wilcox.test_p.values_SDHA$name[wilcox.test_p.values_SDHA$p.adj.BH < 0.05]
#[1] "CLSTN1"  "CPSG3"   "ECEL1"   "EPB41L3" "EPHA5"   "EPN2"

#Гены, для которых полагаем отсутствие статистически значимых различий:
wilcox.test_p.values_SDHA$name[wilcox.test_p.values_SDHA$p.value > 0.05]
#[1] "DDC"    "DPYSL3" "ELAVL4" "FYN"
wilcox.test_p.values_SDHA$name[wilcox.test_p.values_SDHA$p.adj.holm > 0.05]
#[1] "CLSTN1"  "DDC"     "DPYSL3"  "ELAVL4"  "EPB41L3" "EPHA5"   "FYN"
wilcox.test_p.values_SDHA$name[wilcox.test_p.values_SDHA$p.adj.BH > 0.05]
#[1] "DDC"    "DPYSL3" "ELAVL4" "FYN"

#В группах для некоторых генов, наблюдаются различия для значений экспрессий при нормировании по референсному гену SDHA 
#А для некоторых генов нет

#Для генов "CPSG3" "ECEL1" "EPN2" p-value ниже уровня значимости для каждого метода поправок => более вероятно, полагаем, что есть статистически значимые различия между группами
#Для генов "CLSTN1" "EPB41L3" "EPHA5" p-value ниже уровня значимости для стандартного теста и для поправки BH, но не для поправки Holm =>
#=> Предполагаем статистически значимые различия и для них.
#Для генов "DDC" "DPYSL3" "ELAVL4" "FYN" p-value больше уровня значимости =>
#предполагаем отсутствие статистически значимых различий между ними.


#Выбор референсного гена влияет на итоговый результат, поскльку для тех же выбранных генов
#при нормировании по UBS значения для всех геннов в группах имеются различия
#а при нормировании по SDHA различия и их отсутствие зависят от генов. (различия описаны выше)

library(ggplot2)

#Графики значений экспрессий, нормированных по гену UBS

grt_delta_UBS_by_group_name <- group_by(grt_delta_UBS, group, name) %>% select(value, group, name)

ggplot(data = grt_delta_UBS_by_group_name, mapping = aes(x = name, y = value, color = name)) +
  facet_grid(. ~ group) +
  theme_classic() +
  geom_violin() +
  geom_boxplot(width = 0.3) +
  labs(title = "dCt values for the UBS reference for the genes in the 1st and 2nd groups", x = "Gene name", y = "dCt value") +
  scale_color_manual(name = "Gene: Difference source", labels = c("CLSTN1: every adj.", "CPSG3: every adj.",  "DDC: every adj.",    "DPYSL3: every adj.", "ECEL1: every adj.",  "ELAVL4: every adj.", "EPB41L3: every adj.","EPHA5: every adj.",  "EPN2: every adj.","FYN: every adj."), values = rep("red", 10))


#Графики значений экспрессий, нормированных по гену SDHA

grt_delta_SDHA_by_group_name <- group_by(grt_delta_SDHA, group, name) %>% select(value, group, name)

ggplot(data = grt_delta_SDHA_by_group_name, mapping = aes(x = name, y = value, color = name)) +
  facet_grid(. ~ group) +
  theme_classic() +
  geom_violin() +
  geom_boxplot(width = 0.3) +
  labs(title = "dCt values for the SDHA reference for the genes in the 1st and 2nd groups", x = "Gene name", y = "dCt value") +
  scale_color_manual(name = "Gene: Difference source", labels = c("CLSTN1: BH adj.", "CPSG3: every adj.", "ECEL1: every adj.", "EPB41L3: BH adj.", "EPHA5: BH adj.", "EPN2: every adj."), values = c("CLSTN1" = "blue", "CPSG3" = "red", "ECEL1" = "red", "EPB41L3" = "blue", "EPHA5" = "blue", "EPN2" = "red"))
