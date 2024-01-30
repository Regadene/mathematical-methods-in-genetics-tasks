diabetes_df <- read.table(file='diabetes.txt',sep = "\t", skip = 10, header = TRUE)

#9.1
age_by_location_df <- aggregate(x = age ~ location, data = diabetes_df, FUN = c)
age_location_1_vector <- age_by_location_df$age[[1]]
age_location_2_vector <- age_by_location_df$age[[2]]

shapiro.test(age_location_1_vector) #p-value = 0.0005226
shapiro.test(age_location_2_vector) #p-value = 0.001161
#Значения p-value низкие, распределения не нормальные
#значит не выполнено главное условие для использования теста Стьюдента

#Поэтому будем использовать тест Вилкоксона
stat_test_results <- wilcox.test(age_location_1_vector, age_location_2_vector)
#высокая вероятность, что наша гипотеза верна: p-value = 0.7908. А значит, отличия незначительные
age_differences <- FALSE

#9.2
gender_amount_by_location_df <- aggregate(x = id ~ location + gender, data = diabetes_df, FUN = length)
gender_amount_by_location_matrix <- matrix(data = gender_amount_by_location_df$id, nrow = 2, ncol = 2)
MF_ratio_stat_test_results <- chisq.test(gender_amount_by_location_matrix)
#высокая вероятность, что наша гипотеза верна: p-value = 0.7422. А значит, отличия незначительные
MF_ratio_differences <- FALSE
