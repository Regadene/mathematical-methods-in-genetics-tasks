diabetes_df <- read.table(file='diabetes.txt',sep = "\t", skip = 10, header = TRUE)

#7.1
diabetes_M_vector <- diabetes_df$gender == "male"
diabetes_F_vector <- diabetes_df$gender == "female"
diabetes_Small_vector <- diabetes_df$frame == "small"
diabetes_Medium_vector <- diabetes_df$frame == "medium"
diabetes_Large_vector <- diabetes_df$frame == "large"

SBP_M_small <- mean(diabetes_df$bp.1s[diabetes_M_vector & diabetes_Small_vector], na.rm = TRUE)
SBP_M_medium <- mean(diabetes_df$bp.1s[diabetes_M_vector & diabetes_Medium_vector], na.rm = TRUE)
SBP_M_large <- mean(diabetes_df$bp.1s[diabetes_M_vector & diabetes_Large_vector], na.rm = TRUE)

SBP_F_small <- mean(diabetes_df$bp.1s[diabetes_F_vector & diabetes_Small_vector], na.rm = TRUE)
SBP_F_medium <- mean(diabetes_df$bp.1s[diabetes_F_vector & diabetes_Medium_vector], na.rm = TRUE)
SBP_F_large <- mean(diabetes_df$bp.1s[diabetes_F_vector & diabetes_Large_vector], na.rm = TRUE)

DPB_M_small <- mean(diabetes_df$bp.1d[diabetes_M_vector & diabetes_Small_vector], na.rm = TRUE)
DPB_M_medium <- mean(diabetes_df$bp.1d[diabetes_M_vector & diabetes_Medium_vector], na.rm = TRUE)
DPB_M_large <- mean(diabetes_df$bp.1d[diabetes_M_vector & diabetes_Large_vector], na.rm = TRUE)

DPB_F_small <- mean(diabetes_df$bp.1d[diabetes_F_vector & diabetes_Small_vector], na.rm = TRUE)
DPB_F_medium <- mean(diabetes_df$bp.1d[diabetes_F_vector & diabetes_Medium_vector], na.rm = TRUE)
DPB_F_large <- mean(diabetes_df$bp.1d[diabetes_F_vector & diabetes_Large_vector], na.rm = TRUE)

#7.2
diabetes_age_quantile <- quantile(diabetes_df$age)

first_age_quart_diabetes_vector <- diabetes_df$age >= diabetes_age_quantile[[1]] & diabetes_df$age <= diabetes_age_quantile[[2]]
second_age_quart_diabetes_vector <- diabetes_df$age >= diabetes_age_quantile[[2]] & diabetes_df$age <= diabetes_age_quantile[[3]]
third_age_quart_diabetes_vector <- diabetes_df$age >= diabetes_age_quantile[[3]] & diabetes_df$age <= diabetes_age_quantile[[4]]
fourth_age_quart_diabetes_vector <- diabetes_df$age >= diabetes_age_quantile[[4]] & diabetes_df$age <= diabetes_age_quantile[[5]]

stab.glu_1th_quart <- mean(diabetes_df$stab.glu[first_age_quart_diabetes_vector], na.rm = TRUE)
stab.glu_2nd_quart <- mean(diabetes_df$stab.glu[second_age_quart_diabetes_vector], na.rm = TRUE)
stab.glu_3rd_quart <- mean(diabetes_df$stab.glu[third_age_quart_diabetes_vector], na.rm = TRUE)
stab.glu_4th_quart <- mean(diabetes_df$stab.glu[fourth_age_quart_diabetes_vector], na.rm = TRUE)

#7.3
chol_min <- min(diabetes_df$chol, na.rm = TRUE)
chol_max <- max(diabetes_df$chol, na.rm = TRUE)

id_chol_min_with_NA <- diabetes_df$id[diabetes_df$chol == chol_min]
id_chol_max_with_NA <- diabetes_df$id[diabetes_df$chol == chol_max]

id_chol_min <- id_chol_min_with_NA[complete.cases(id_chol_min_with_NA)]
id_chol_max <- id_chol_max_with_NA[complete.cases(id_chol_max_with_NA)]

#7.4
aggregate(x = diabetes_df$age, by = list("gender" = diabetes_df$gender, "frame" = diabetes_df$frame), FUN = mean)
# Не совсем ясно, что значит "оценить средний возраст в полученных группах"
# Оценка: средний возраст испытуемых в диапазоне от 42 до 54 лет.
# Средний возраст мужчин от 42 и до 54 лет, а женщин от 42 и до 52 лет.

#7.5
diabetes_is_col_numeric <- sapply(X = diabetes_df, FUN = is.numeric)
means_df <- aggregate(x = diabetes_df[diabetes_is_col_numeric], by = list("location" = diabetes_df$location), FUN = mean, na.rm = TRUE)

#7.6
character_column <- sapply(X = diabetes_df, FUN = is.character)

#7.7
medians <- apply(X = diabetes_df[diabetes_is_col_numeric], MARGIN = 2, FUN = median, na.rm = TRUE)