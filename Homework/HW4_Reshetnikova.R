library(ggplot2)

# Подготовка данных
qrt_groups <- read.table("sample_groups_4.tsv", head = T)
qrt_data <- read.table("QRT_PCR_LR.txt", head = T)
qrt <- merge(qrt_groups, qrt_data, by.x = "sample", by.y = "tissue")
qrt$groups_4 <- as.factor(qrt$groups_4)

# ΔCq используется для определения разницы в экспрессии генов между разными образцами.
# Он вычисляется вычитанием Cq значения эталонного гена из Cq значения целевого гена
# Используем UBC в качестве референса, вычитаем значения UBC

qrt_delta_UBS_df <- qrt
qrt_delta_UBS_df[, -c(1, 2)] <- qrt[, -c(1, 2)] - qrt$UBC

#4.1
# Выберем гены MAPT и PLAT

# Проверка нормальности распределения данных
shapiro.test_p.value_fun <- function(x) { shapiro.test(x)$p.value }

shapiro_test_MAPT <- aggregate(MAPT ~ groups_4, data = qrt_delta_UBS_df, FUN = shapiro.test_p.value_fun)
length(shapiro_test_MAPT[shapiro_test_MAPT$MAPT > 0.05]) == length(shapiro_test_MAPT$MAPT)
# Данные в группах не распределены нормально

shapiro_test_PLAT <- aggregate(PLAT ~ groups_4, data = qrt_delta_UBS_df, FUN = shapiro.test_p.value_fun)
length(shapiro_test_PLAT[shapiro_test_PLAT$PLAT > 0.05]) == length(shapiro_test_PLAT$PLAT)
# Данные в группах не распределены нормально

bartlett.test(x = qrt_delta_UBS_df$MAPT, g = qrt_delta_UBS_df$groups_4)
# p-value меньше уровня значимости (< 2.2e-16) => дисперсии различаются

bartlett.test(x = qrt_delta_UBS_df$PLAT, g = qrt_delta_UBS_df$groups_4)
# p-value меньше уровня значимости (< 2.2e-16) => дисперсии различаются

# Для выполнения дисперсионного анализа нам необходимо, чтобы данные были
# распределены нормально и дисперсии были равны
# поэтому будем использовать тест Краскала-Уоллиса

kruskal.test(MAPT ~ groups_4 , data = qrt_delta_UBS_df)
# p-value меньше уровня значимости (1.442e-12) => существуют значимые различия по ΔCq между группами

kruskal.test(PLAT ~ groups_4 , data = qrt_delta_UBS_df)
# p-value меньше уровня значимости (0.0008746) => существуют значимые различия по ΔCq между группами

#4.2
# С помощью post-hoc тестов определим, какие именно группы отличаются и в какую сторону.

anova_MAPT <- aov(formula = MAPT ~ groups_4, data = qrt_delta_UBS_df)
post_hoc_MAPT <- TukeyHSD(anova_MAPT)
plot(post_hoc_MAPT, las = 1)
# 2-1 - ΔCq больше у второй группы
# 3-1 - нет значимых различий
# 4-1 - нет значимых различий
# 3-2 - ΔCq больше у второй группы
# 4-2 - нет значимых различий
# 4-3 - нет значимых различий

anova_PLAT <- aov(formula = PLAT ~ groups_4, data = qrt_delta_UBS_df)
post_hoc_PLAT <- TukeyHSD(anova_PLAT)
plot(post_hoc_PLAT, las = 1)
# 2-1 - нет значимых различий
# 3-1 - нет значимых различий
# 4-1 - ΔCq больше у первой группы
# 3-2 - нет значимых различий
# 4-2 - ΔCq больше у второй группы
# 4-3 - ΔCq больше у третьей группы

#4.3
# Графики значений ΔCq
error_max_fun <- function(x) mean(x) + sd(x)/sqrt(length(x))
error_min_fun <- function(x) mean(x) - sd(x)/sqrt(length(x))

# Графики значений ΔCq MAPT и UBS
ggplot(qrt_delta_UBS_df, aes(x = groups_4, y = MAPT)) +
  xlab("Group number") +
  ylab("delCq MAPT_UBS") +
  geom_violin() +
  geom_jitter(width = 0.1) +
  stat_summary(fun = mean, geom = "errorbar", fun.max = error_max_fun, fun.min = error_min_fun, color = "red", width = 0.2) +
  stat_summary(fun = mean, geom = "point", color = "red") +
  theme_bw()
# На графике видно, что медианы для 1 и 3 группы примерно равны
# медиана 2 группы заметно больше медиан 1 и 3 групп
# в 4 группе довольно мало данных и их разброс большой,
# медиана между значениями для 2 группы и 1, 3 групп,
# но разброс ошибки сильно больше => однозначно дать ответ сложно,
# считаем, что значимых различий между 4 группой и остальными нет

# Графики значений ΔCq PLAT и UBS
ggplot(qrt_delta_UBS_df, aes(x = groups_4, y = PLAT)) +
  xlab("Group number") +
  ylab("delCq PLAT_UBS") +
  geom_violin() +
  geom_jitter(width = 0.1) +
  stat_summary(fun = mean, geom = "errorbar", fun.max = error_max_fun, fun.min = error_min_fun, color = "red", width = 0.2) +
  stat_summary(fun = mean, geom = "point", color = "red") +
  theme_bw()
# На графике видно, что медианы для 1 и 3 группы примерно равны
# а значимых различий между медианой с ошибкой для 2 группы и для 1 и 3 нет
# отчетливо видно, что значения медиан и ошибки для 4 группы меньше,
# чем для всех остальных групп

# Имеем, различия в ΔCq между группами, полученные по анализу post-hoc тестов и
# на графиках соответствуют друг другу

