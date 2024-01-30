library(ggplot2)

# Загружаем таблицу в R:
qrt_pcr_rye <- read.table("data.tsv", head = T)

# Рассчитываю ∆CT. Используем reference_gene в качестве референса.
qrt_pcr_rye_delta <- qrt_pcr_rye
qrt_pcr_rye_delta[, -c(5, 6)] <- qrt_pcr_rye_delta[, -c(5, 6)] - qrt_pcr_rye_delta$reference_gene
# Создаём факторы для линий и цветов:
qrt_pcr_rye_delta$line <- as.factor(qrt_pcr_rye_delta$line)
qrt_pcr_rye_delta$color <- as.factor(qrt_pcr_rye_delta$color)

# Сравнение экспрессии генов в зависимости от линии:
# Проведем дисперсионный анализ. Проверка на нормальность данных для генов F3H, DFR, ANS. 
shapiro.test_p.value_fun <- function(x) { shapiro.test(x)$p.value }

shapiro_test_F3H <- aggregate(F3H ~ line, data = qrt_pcr_rye_delta, FUN = shapiro.test_p.value_fun)
# Для гена F3H 
#    line         F3H
# 1 line1 0.599464724 - норм. распределение
# 2 line2 0.258073980 - норм. распределение
# 3 line3 0.817154285 - норм. распределение
# 4 line4 0.002298943 - не норм. распределение
# 5 line5 0.026194917 - не норм. распределение
# 6 line6 0.069972115 - норм. распределение
# 7 line7 0.494993564 - норм. распределение

shapiro_test_DFR <- aggregate(DFR ~ line, data = qrt_pcr_rye_delta, FUN = shapiro.test_p.value_fun)
# Для гена DFR
#    line        DFR
# 1 line1 0.55567002 - норм. распределение
# 2 line2 0.49814106 - норм. распределение
# 3 line3 0.45967333 - норм. распределение
# 4 line4 0.02415555 - ненорм. распределение
# 5 line5 0.44348075 - норм. распределение
# 6 line6 0.42001400 - норм. распределение
# 7 line7 0.14412558 - норм. распределение

shapiro_test_ANS <- aggregate(ANS ~ line, data = qrt_pcr_rye_delta, FUN = shapiro.test_p.value_fun)
# Для гена ANS
#    line         ANS
# 1 line1 0.008772804 - ненорм. распределение
# 2 line2 0.819011083 - норм. распределение
# 3 line3 0.018959032 - ненорм. распределение
# 4 line4 0.003010515 - ненорм. распределение
# 5 line5 0.994697754 - норм. распределение
# 6 line6 0.070845590 - норм. распределение
# 7 line7 0.512757904 - норм. распределение

# Хоть мы и получили ненормальное распределение и будем использовать непараметрические тесты,
# дополнительно можем провести barlett.test для оценки равенства дисперсий. 

bartlett.test(x = qrt_pcr_rye_delta$F3H, g = qrt_pcr_rye_delta$line)
# p-value меньше уровня значимости (p-value = 1.577e-05) => дисперсии различаются

bartlett.test(x = qrt_pcr_rye_delta$DFR, g = qrt_pcr_rye_delta$line)
# p-value больше уровня значимости (p-value = 0.3482) => дисперсии равны

bartlett.test(x = qrt_pcr_rye_delta$ANS, g = qrt_pcr_rye_delta$line)
# p-value больше уровня значимости (p-value = 0.1609) => дисперсии равны 

#Данные не распределены нормально. Используем непараметрический тест Краскала-Уоллиса

kruskal.test(F3H ~ line , data = qrt_pcr_rye_delta)
# p-value больше уровня значимости (p-value = 0.2984) => нет значимых различий по ΔCq между линиями


kruskal.test(DFR ~ line , data = qrt_pcr_rye_delta)
# p-value больше уровня значимости (p-value = 0.6211) => нет значимых различий по ΔCq между линиями


kruskal.test(ANS ~ line , data = qrt_pcr_rye_delta)
# p-value меньше уровня значимости (p-value = 0.8145) => нет значимых различий по ΔCq между линиями


# Помимо теста Краскала-Уоллиса проведем тест Mann-Whitney U-тест (pairwise.wilcox.test). По умолчанию используется поправка Холма, но 
# различий в результатах попарных сравнений при её использовании почти не было. 
# Для DFR, например, все результаты были равны единице. Поэтому использовала другую поправку Бенджамини-Хохмана.
# Используя её, получаем более точные результаты:

pairwise.wilcox.test(x = qrt_pcr_rye_delta$F3H, g = qrt_pcr_rye_delta$line, p.adjust.method = "BH")

#         line1 line2 line3 line4 line5 line6
#   line2 0.827 -     -     -     -     -    
#   line3 0.755 0.904 -     -     -     -    
#   line4 0.904 0.904 0.904 -     -     -    
#   line5 0.755 0.755 0.091 0.812 -     -    
#   line6 0.904 0.904 0.812 0.937 0.432 -    
#   line7 0.937 0.827 0.812 0.904 0.904 0.904
# Статистически значимых различий -- нет. Но между 3 и 5 линией можем полагать существование хоть каких-то
# (но статистические незначимых) различий (p-value = 0.091)


pairwise.wilcox.test(x = qrt_pcr_rye_delta$DFR, g = qrt_pcr_rye_delta$line, p.adjust.method = "BH")

#         line1 line2 line3 line4 line5 line6
#   line2 0.94  -     -     -     -     -    
#   line3 0.92  0.92  -     -     -     -    
#   line4 0.94  0.94  0.94  -     -     -    
#   line5 0.92  0.84  0.94  0.84  -     -    
#   line6 0.94  0.94  0.84  0.94  0.84  -    
#   line7 0.94  0.94  0.84  0.94  0.84  0.94 
# Статистически значимых различий -- нет.

pairwise.wilcox.test(x = qrt_pcr_rye_delta$ANS, g = qrt_pcr_rye_delta$line, p.adjust.method = "BH")
#         line1 line2 line3 line4 line5 line6
#   line2 0.98  -     -     -     -     -    
#   line3 0.98  0.98  -     -     -     -    
#   line4 0.98  0.98  0.98  -     -     -    
#   line5 0.98  0.98  0.98  0.98  -     -    
#   line6 0.98  0.98  0.98  0.98  0.98  -    
#   line7 0.98  0.98  0.98  0.98  1.00  0.98 
# Статистически значимых различий -- нет.


# Графики значений ΔCt
error_max_fun <- function(x) mean(x) + sd(x)/sqrt(length(x))
error_min_fun <- function(x) mean(x) - sd(x)/sqrt(length(x))

# Графики значений ΔCt F3H
ggplot(qrt_pcr_rye_delta, aes(x = line, y = F3H)) +
  xlab("Line number") +
  ylab("delCt F3H") +
  geom_violin() +
  geom_jitter(width = 0.1) +
  stat_summary(fun = mean, geom = "errorbar", fun.max = error_max_fun, fun.min = error_min_fun, color = "red", width = 0.2) +
  stat_summary(fun = median, geom = "point", color = "red") +
  theme_bw()

# Графики значений ΔCt DFR
ggplot(qrt_pcr_rye_delta, aes(x = line, y = DFR)) +
  xlab("Line number") +
  ylab("delCt DFR") +
  geom_violin() +
  geom_jitter(width = 0.1) +
  stat_summary(fun = mean, geom = "errorbar", fun.max = error_max_fun, fun.min = error_min_fun, color = "red", width = 0.2) +
  stat_summary(fun = median, geom = "point", color = "red") +
  theme_bw()

# Графики значений ΔCt ANS
ggplot(qrt_pcr_rye_delta, aes(x = line, y = ANS)) +
  xlab("Line number") +
  ylab("delCt ANS") +
  geom_violin() +
  geom_jitter(width = 0.1) +
  stat_summary(fun = mean, geom = "errorbar", fun.max = error_max_fun, fun.min = error_min_fun, color = "red", width = 0.2) +
  stat_summary(fun = median, geom = "point", color = "red") +
  theme_bw()


# Сравнение экспрессии генов в зависимости от цвета
# Проверка на нормальность данных для генов F3H, DFR, ANS. 
shapiro_test_F3H_color <- aggregate(F3H ~ color, data = qrt_pcr_rye_delta, FUN = shapiro.test_p.value_fun)

#   color           F3H
# 1 brown  0.2751092130 - норм. распределение
# 2 yellow 0.0003272108 - не норм. распределение

shapiro_test_DFR_color <- aggregate(DFR ~ color, data = qrt_pcr_rye_delta, FUN = shapiro.test_p.value_fun) 

#   color         DFR
# 1 brown  0.87358215 - норм. распределение
# 2 yellow 0.03208932 - не норм. распределение

shapiro_test_ANS_color <- aggregate(ANS ~ color, data = qrt_pcr_rye_delta, FUN = shapiro.test_p.value_fun)

#   color          ANS
# 1 brown  0.001162996 - не норм. распределение
# 2 yellow 0.061563677 - норм. распределение

# Данные распределены ненормально, будем использовать непараметрический тест Вилкоксона

qrt_genes_names <- c("F3H", "DFR", "ANS")
qrt_pcr_rye_delta_tidy <- pivot_longer(data = qrt_pcr_rye_delta, cols = qrt_genes_names)
wilcox.test_p.values <- qrt_pcr_rye_delta_tidy %>%
  group_by(name) %>%
  summarise(p.value = wilcox.test(value ~ color)$p.value)

# Поправки на множественные сравнения
wilcox.test_p.values$p.adj.holm <- p.adjust(p = wilcox.test_p.values$p.value, method = "holm")
wilcox.test_p.values$p.adj.BH <- p.adjust(p = wilcox.test_p.values$p.value, method = "BH")
#   name  p.value p.adj.holm p.adj.BH
# 1 ANS     0.606          1    0.909
# 2 DFR     0.930          1    0.930
# 3 F3H     0.357          1    0.909
# Статистически значимых различий нет

#Графики значений ΔCt, нормированных по референсному гену:
qrt_delta_by_name_color <- group_by(qrt_pcr_rye_delta_tidy, name, color) %>% select(value, name, color)

ggplot(data = qrt_delta_by_name_color, mapping = aes(x = color, y = value)) +
  facet_grid(. ~ name) +
  theme_classic() +
  geom_violin() +
  geom_boxplot(width = 0.3) +
  labs(title = "dCt values for the anthocyan genes by colors", x = "Color", y = "dCt value")

