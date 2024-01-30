#1
#Arabidopsis thaliana
#Интервальный - высота растения
#Шкальный - короткостебельные, нормальные, длинностебельные
#Категориальный - Окраска лепестков
#Drosophila melanogaster
#Интервальный - длина тела
#Шкальный - короткощетинистые, нормальные  
#Категориальный - закрученные крылья/незакрученные 
#Mus musculus
#Интервальный - масса
#Шкальный - короткохвостые, нормальные, длиннохвостые
#Категориальный - пол: мужской женский
#Saccharomyces cerevisiae
#Интервальный - количество колоний в чашке Петри
#Шкальный - размеры колоний
#Категориальный - формы колоний

diabetes_df <- read.table(file="diabetes.txt",sep = "\t", skip = 10, header = TRUE)
#2
males_df <- diabetes_df[diabetes_df$gender == "male", ]
males_average_time.ppn <- mean(males_df$time.ppn, na.rm = T)
males_median_time.ppn <- median(males_df$time.ppn, na.rm = T)
#Я считаю, для описания выборки стоит использовать медиану. 
#На гистограмме ниже видно, что в районе значения 750 большой выброс данных.
#А среди значений от 0 до 250-300 собрано довольно много значений.
#Поэтому среднее может быть более плохим выбором, поскольку выброс в районе 750 может влиять на среднее довольно сильно.
#time.ppn - интервальный тип данных

#3
hist(x = males_df$time.ppn, breaks = 30, ylim = c(0, 22))
abline(v = males_average_time.ppn, lwd = 4, col = 2)
text(x=males_average_time.ppn, y=21, "average")
abline(v = males_median_time.ppn, lwd = 4, col = 3)
text(x=males_median_time.ppn, y=22, "median")

#4
males_average_bp.1d <- mean(males_df$bp.1d, na.rm = T)
males_median_bp.1d <- median(males_df$bp.1d, na.rm = T)

hist(x = males_df$bp.1d, breaks = 30, ylim = c(0, 22))
abline(v = mean(males_average_bp.1d, na.rm = T), lwd = 4, col = 2)
text(x=males_average_bp.1d, y=21, "average")
abline(v = median(males_median_bp.1d, na.rm = T), lwd = 4, col = 3)
text(x=males_median_bp.1d, y=22, "median")

#5
males_sd_bp.1d <- sd(males_df$bp.1d, na.rm = T)
males_se_bp.1d <- males_sd_bp.1d/sqrt(length(males_df$bp.1d[!is.na(males_df$bp.1d)]))

#6
buckingham_patients <- diabetes_df[diabetes_df$location == "Buckingham", ]
buckingham_males <- buckingham_patients[buckingham_patients$gender == "male", ]
buckingham_males_percentage <- 100 * nrow(buckingham_males) / nrow(buckingham_patients)
buckingham_males_percentage_error <- sqrt(buckingham_males_percentage * (100 - buckingham_males_percentage)/nrow(buckingham_patients))
