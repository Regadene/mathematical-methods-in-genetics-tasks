library(tidyr)
library(ggplot2)

#3.1
diabetes_df <- read.table(file = 'diabetes.txt', sep = "\t", skip = 10, head = TRUE)
diabetes_df_complete.cases <- diabetes_df[complete.cases(diabetes_df[, c("id", "gender", "chol", "stab.glu", "hdl", "glyhb")]), c("id", "gender", "chol", "stab.glu", "hdl", "glyhb")]

shapiro.test_p.value_fun <- function(x) { shapiro.test(x)$p.value > 0.05 }

diabetes_df_complete.cases_shapiro_reses <- apply(diabetes_df_complete.cases[ , c("chol", "stab.glu", "hdl", "glyhb")], MARGIN = 2, FUN = shapiro.test_p.value_fun)
# Для каджого из "chol", "stab.glu", "hdl", "glyhb", p.value ниже уровня значимости
# => отклоняем нулевую гипотезу и пологаем, что данные не распределены нормально.

#Данные не соответствуют нормальному распределению, следует использовать непараметрические тесты - коэффициенты корреляции Кендалла или Спирмена
kendall_cor <- cor(x = diabetes_df_complete.cases$chol, y = diabetes_df_complete.cases[ , c( "stab.glu", "glyhb", "hdl")], method = "kendall", use = "pairwise.complete.obs")
spearman_cor <- cor(x = diabetes_df_complete.cases$chol, y = diabetes_df_complete.cases[ , c( "stab.glu", "glyhb", "hdl")], method = "spearman", use = "pairwise.complete.obs")

cor.test(x = diabetes_df_complete.cases$chol, y = diabetes_df_complete.cases[ , c( "stab.glu")], method = "kendall", use = "pairwise.complete.obs")$p.value > 0.05
cor.test(x = diabetes_df_complete.cases$chol, y = diabetes_df_complete.cases[ , c( "glyhb")], method = "kendall", use = "pairwise.complete.obs")$p.value > 0.05
cor.test(x = diabetes_df_complete.cases$chol, y = diabetes_df_complete.cases[ , c( "hdl")], method = "kendall", use = "pairwise.complete.obs")$p.value > 0.05

#p.value ниже уровня значимости для "stab.glu", "hdl", "glyhb"
# => отклоняем нулевую гипотезу и пологаем существование корреляции


#Тепловая карта корреляционной матрицы по методу Кенделла
kendall_cor_df <- as.data.frame(kendall_cor)
kendall_cor_df <- pivot_longer(data = kendall_cor_df, cols = c("stab.glu", "glyhb", "hdl"))
ggplot(data = kendall_cor_df, mapping = aes(x = rep(1, 3), y = name, fill = value)) + xlab("chol") + geom_tile() + theme_bw()

#Тепловая карта корреляционной матрицы по методу Спирмана
spearman_cor_df <- as.data.frame(spearman_cor)
spearman_cor_df <- pivot_longer(data = spearman_cor_df, cols = c("stab.glu", "glyhb", "hdl"))
ggplot(data = spearman_cor_df, mapping = aes(x = rep(1, 3), y = name, fill = value)) + xlab("chol") + geom_tile() + theme_bw()

#Диаграммы рассеяния
plot(chol ~ stab.glu, data = diabetes_df_complete.cases)
plot(chol ~ hdl, data = diabetes_df_complete.cases)
plot(chol ~ glyhb, data = diabetes_df_complete.cases)

#3.2
#Регрессионные модели 
chol_model <- lm(formula = chol ~ glyhb + hdl + stab.glu, data = diabetes_df_complete.cases)
summary(chol_model)$coefficients
summary(chol_model)$r.squared

#Анализ регрессионных остатков
shapiro.test(x = chol_model$residuals)$p.value > 0.05
# p.value в shapiro.test для остатков ниже уровня значимости
# => отклоняем нулевую гипотезу и пологаем, что данные не распределены нормально.

#Изначально - 389 строк без NA в "chol", "stab.glu", "hdl", "glyhb"
#Ищем выбросы на графиках и исключаем их, добавляя новые условия в фильтрацию
#(я перешла на исключение интервалами, так как долго не могла прийти к значениям близким к уровню значимости) 
#& (diabetes_df_complete.cases$id < 20720 | diabetes_df_complete.cases$id > 20820)

#Исключаем выбросы
diabetes_df_complete.cases_normal <- diabetes_df_complete.cases[ !(diabetes_df_complete.cases$id %in% c(60, 141, 285, 2780, 2785, 2787, 2778, 2791, 2793, 2794, 2795, 3250, 3750, 3751, 3752, 4000, 12778, 12769, 13254, 13501, 13505, 14758, 15010, 15012, 15013, 15278, 15519, 15800, 15802, 15805, 20318, 20332, 20340, 20350, 20361, 20754, 41003, 41023, 41034)) & (diabetes_df_complete.cases$id < 15740 | diabetes_df_complete.cases$id > 15773) & (diabetes_df_complete.cases$id < 15500 | diabetes_df_complete.cases$id > 15550) & (diabetes_df_complete.cases$id < 20720 | diabetes_df_complete.cases$id > 20820) & (diabetes_df_complete.cases$id < 21250 | diabetes_df_complete.cases$id > 21500) & (diabetes_df_complete.cases$id < 20250 | diabetes_df_complete.cases$id > 20315) & (diabetes_df_complete.cases$id < 40400 | diabetes_df_complete.cases$id > 40600), ]
plot(chol_model, which = 1, labels.id = diabetes_df_complete.cases_normal$id)
plot(chol_model, which = 2, labels.id = diabetes_df_complete.cases_normal$id)
plot(chol_model, which = 4, labels.id = diabetes_df_complete.cases_normal$id)

chol_model_normal <- lm(formula = chol ~ glyhb + hdl + stab.glu, data = diabetes_df_complete.cases_normal)
summary(chol_model_normal)$coefficients
summary(chol_model_normal)$r.squared

#Для 267 строк в diabetes_df_complete.cases_normal мы уже имеем нормальное распределение остатков
shapiro.test(x = chol_model_normal$residuals)$p.value > 0.05
#p.value выше уровня значимости => полагаем, что распределение нормальное
t.test(chol_model_normal$residuals)$p.value
#p.value == 1 => полагаем, что среднее значение этого распределения равно 0
cor.test(chol_model_normal$residuals, diabetes_df_complete.cases_normal$glyhb, method = 'kendall')$p.value > 0.05
cor.test(chol_model_normal$residuals, diabetes_df_complete.cases_normal$hdl, method = 'kendall')$p.value > 0.05
cor.test(chol_model_normal$residuals, diabetes_df_complete.cases_normal$stab.glu, method = 'kendall')$p.value > 0.05
#p.value для всех корреляционных тестов выше уровня значимости 
#=> полагаем, что нет взаимосвязи между значением остатка и независимыми переменными

#3.3
summary(chol_model_normal)$coefficients[2, 4] > 0.05
summary(chol_model_normal)$coefficients[4, 4] > 0.05
#p.value для glyhb и stab.glu больше уровня значимости => полагаем, что glyhb и stab.glu не влияют на chol значительно

#Исключаем параметр stab.glu из модели
chol_model_by_glyhb_and_hdl <- lm(formula = chol ~ glyhb + hdl, data = diabetes_df_complete.cases_normal)
summary(chol_model_by_glyhb_and_hdl)$coefficients
summary(chol_model_by_glyhb_and_hdl)$r.squared

shapiro.test(x = chol_model_by_glyhb_and_hdl$residuals)$p.value
t.test(chol_model_by_glyhb_and_hdl$residuals)$p.value
cor.test(chol_model_by_glyhb_and_hdl$residuals, diabetes_df_complete.cases_normal$glyhb, method = 'kendall')$p.value > 0.05
cor.test(chol_model_by_glyhb_and_hdl$residuals, diabetes_df_complete.cases_normal$hdl, method = 'kendall')$p.value > 0.05
cor.test(chol_model_by_glyhb_and_hdl$residuals, diabetes_df_complete.cases_normal$stab.glu, method = 'kendall')$p.value > 0.05
#При исключении stab.glu из модели значение r.squared почти не изменилось
#(0.1086591 -> 0.1075446)
#=>параметр stab.glu не вносит значительного вклада в объяснение изменчивости chol
anova(chol_model_normal, chol_model_by_glyhb_and_hdl)
#p.value сильно выше уровня значимости (0.5668) => полагаем, что влияние stab.glu незначительно, можно его отбросить

#Рассмотрим исключение hdl из модели
chol_model_by_glyhb_and_stab.glu <- lm(formula = chol ~ glyhb + stab.glu, data = diabetes_df_complete.cases_normal)
summary(chol_model_by_glyhb_and_stab.glu)$r.squared
#Исключение hdl сильно влияет на значение r.squared
#(0.1086591 -> 0.03012276)
anova(chol_model_normal, chol_model_by_glyhb_and_stab.glu)
#p.value ниже уровня значимости => полагаем, что влияние hdl значительно

#Рассмотрим исключение glyhb из модели
chol_model_by_stab.glu_and_hdl <- lm(formula = chol ~ stab.glu + hdl, data = diabetes_df_complete.cases_normal)
summary(chol_model_by_stab.glu_and_hdl)$r.squared
#Исключение glyhb влияет на значение r.squared, сильно меньше, чем исключение hdl, но сильнее, чем исключение stab.glu
#В целом, влияние можно считать малозначительным
#(0.1086591 -> 0.096773)
anova(chol_model_normal, chol_model_by_stab.glu_and_hdl)
#p.value выше уровня значимости (0.06221) => полагаем, что влияние glyhb малозначительно, можно отбросить и этот параметр

#Рассмотрим исключение stab.glu и glyhb из модели
chol_model_by_hdl <- lm(formula = chol ~ hdl, data = diabetes_df_complete.cases_normal)
summary(chol_model_by_hdl)$r.squared
#Исключение glyhb вместе с stab.glu влияет на значение r.squared уже значительно
#(0.1086591 -> 0.06164308)
anova(chol_model_normal, chol_model_by_hdl)
#p.value меньше уровня значимости (0.00116) => полагаем, что влияние пары glyhb и stab.glu значительно
#не стоит отбрасывать оба этих параметра

#Считаем, что оптимальная модель: chol ~ glyhb + hdl

#3.4 
#При проведении корреляционного анализа мы заключили,
#что корреляция существует для каждого из параметра "stab.glu", "hdl", "glyhb"
#При проведении регрессионного анализа мы определили достоверное влияние для параметра "hdl",
#параметр "glyhb" не объясняет изменчивость параметра "chol" достаточно достоверно,
#но в отличие "stab.glu" значения p.value довольно близки к уровню значимости, так что однозначный вывод дать сложно.
#параметр "stab.glu" мы считаем не объясняющим изменчивость параметра "chol" достоверно

#3.5
#Строим для мужчин
diabetes_df_complete.cases_normal_male <- diabetes_df_complete.cases_normal[diabetes_df_complete.cases_normal$gender=='male', ]

chol_model_normal_male <- lm(formula = chol ~ glyhb + hdl, data = diabetes_df_complete.cases_normal_male)
summary(chol_model_normal_male)$coefficients
summary(chol_model_normal_male)$r.squared

shapiro.test(x = chol_model_normal_male$residuals)$p.value > 0.05
#p.value выше уровня значимости => полагаем, что распределение нормальное
t.test(chol_model_normal_male$residuals)$p.value
#p.value == 1 => полагаем, что среднее значение этого распределения равно 0
cor.test(chol_model_normal_male$residuals, diabetes_df_complete.cases_normal_male$glyhb, method = 'kendall')$p.value > 0.05
cor.test(chol_model_normal_male$residuals, diabetes_df_complete.cases_normal_male$hdl, method = 'kendall')$p.value > 0.05
#p.value для всех корреляционных тестов выше уровня значимости 
#=> полагаем, что нет взаимосвязи между значением остатка и независимыми переменными

#Оптимальная модель: chol ~ glyhb + hdl
#Строим для женщин
diabetes_df_complete.cases_female <- diabetes_df_complete.cases_normal[diabetes_df_complete.cases_normal$gender=='female', ]

chol_model_female <- lm(formula = chol ~ glyhb + hdl, data = diabetes_df_complete.cases_female)
summary(chol_model_female)$coefficients
summary(chol_model_female)$r.squared

shapiro.test(x = chol_model_female$residuals)$p.value > 0.05
#p.value ниже уровня значимости => полагаем, что данные не распределены нормально
#=>полученная на этих данных регрессионная модель недостаточно информативна для заключения выводов (p.value 0.03314781)

#Исключим выбросы из данных, для проведения более информативного сравнения моделей для женщин и мужчин
diabetes_df_complete.cases_normal_female <- diabetes_df_complete.cases_female[ !(diabetes_df_complete.cases_female$id %in% c(10001, 10012, 10014, 10016, 12002, 12004, 12005, 12006, 12501, 12502, 12507)), ]
plot(chol_model_female, which = 1, labels.id = diabetes_df_complete.cases_normal_female$id)
plot(chol_model_female, which = 2, labels.id = diabetes_df_complete.cases_normal_female$id)
plot(chol_model_female, which = 4, labels.id = diabetes_df_complete.cases_normal_female$id)

#Построим модель с исключенными выбросами
chol_model_normal_female <- lm(formula = chol ~ glyhb + hdl, data = diabetes_df_complete.cases_normal_female)
summary(chol_model_normal_female)$coefficients
summary(chol_model_normal_female)$r.squared

shapiro.test(x = chol_model_normal_female$residuals)$p.value > 0.05
#p.value выше уровня значимости => полагаем, что распределение нормальное
t.test(chol_model_normal_female$residuals)$p.value
#p.value == 1 => полагаем, что среднее значение этого распределения равно 0
cor.test(chol_model_normal_female$residuals, diabetes_df_complete.cases_normal_female$glyhb, method = 'kendall')$p.value > 0.05
cor.test(chol_model_normal_female$residuals, diabetes_df_complete.cases_normal_female$hdl, method = 'kendall')$p.value > 0.05
#p.value для всех корреляционных тестов выше уровня значимости 
#=> полагаем, что нет взаимосвязи между значением остатка и независимыми переменными

#Сравним две модели
b1 <- summary(chol_model_normal_male)$coefficients[2, 1]
b1_e <- summary(chol_model_normal_male)$coefficients[2, 2]
b2 <- summary(chol_model_normal_female)$coefficients[2, 1]
b2_e <- summary(chol_model_normal_female)$coefficients[2, 2]
t = (b1 - b2)/sqrt(b1_e^2 + b2_e^2)
t
pt(q = abs(t), df = length(diabetes_df_complete.cases_normal_male$chol) + length(diabetes_df_complete.cases_normal_female$chol) - 2, lower.tail = F) * 2 > 0.05
#p.value больше уровня значимости (0.8486665) => полагаем, что модели не отличаются статистически значимо

#Графики для моделей
plot(chol_model_normal_male)
plot(chol_model_normal_female)