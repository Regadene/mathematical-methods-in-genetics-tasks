diabetes_df <- read.table(file='diabetes.txt',sep = "\t", skip = 10, header = TRUE)

#8.1
Age <- diabetes_df$age
hist(x = Age, breaks = seq(0, 100, 5), labels = TRUE, col="white")

#8.2
gender_frame_pressure_df <- aggregate(x = list(bp.1s = diabetes_df$bp.1s, bp.1d = diabetes_df$bp.1d), by = list("gender" = diabetes_df$gender, "frame" = diabetes_df$frame), FUN = mean, na.rm=TRUE)
gender_frame_pressure_df$genderAndFrame <- paste(gender_frame_pressure_df$gender, '\n', gender_frame_pressure_df$frame)
gender_frame_pressures_matrix <- matrix(data = gender_frame_pressure_df$bp.1s, nrow = 1)
gender_frame_pressures_matrix <- rbind(gender_frame_pressures_matrix, gender_frame_pressure_df$bp.1d)
barplot(gender_frame_pressures_matrix, ylab = "BP", names = gender_frame_pressure_df$genderAndFrame, ylim = c(0, 160), beside = TRUE, legend.text = c("SBP", "DBP"), args.legend = list(x = "topleft", bty = "n", horiz = TRUE))

#8.3
b_pressures_by_gender_df <- aggregate(x = list(bp.1s = diabetes_df$bp.1s, bp.2s = diabetes_df$bp.2s), by = list("gender" = diabetes_df$gender), FUN = c)

boxplot_data_list <- list("female_first_b_pressures" <- b_pressures_by_gender_df$bp.1s[[1]],
"female_second_b_pressures" <- b_pressures_by_gender_df$bp.2s[[1]],
"male_first_b_pressures" <- b_pressures_by_gender_df$bp.1s[[2]],
"male_second_b_pressures" <- b_pressures_by_gender_df$bp.2s[[2]])

boxplot(x = boxplot_data_list, notch = TRUE, col = c("white", "turquoise"), names = c("female", "female", "male", "male"), ylab = "SBP", lwd = 4)
stripchart(x = boxplot_data_list, add = TRUE, pch=15, vertical = T, method = "jitter", jitter = 0.2, cex=0.5)
legend(x = 2.6, y = 250, bty="n", c("1st meas.", "2nd meas."), fill = c("black", "turquoise"))