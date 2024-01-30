library(dplyr)

diabetes_df <- read.table(file='diabetes.txt',sep = "\t", skip = 10, header = TRUE)

#10.1
group_by(diabetes_df, gender, frame) %>%
  select(gender, frame, age) %>%
  summarise(mean(age))

#10.2
medians <- diabetes_df %>%
  select(-location, -gender, -frame) %>%
  summarise_each(funs(mean(., na.rm = TRUE))) %>%
  unlist()
#В задаче было написано сохранить в вектор, а как я понимаю на последнем шаге мы имеем таблицу с одной строкой.
#Нашла команду unlist(), чтобы преобразовать в вектор

library(ggplot2)

#10.3
age_by_location_df <- group_by(diabetes_df, location) %>% select(age, location)
age_by_location_plot <- ggplot(data = age_by_location_df, mapping = aes(x = location, y = age))
age_by_location_plot + geom_violin() + geom_boxplot(width = 0.5)

#10.4
bp.1s_by_location_gender_df <- group_by(diabetes_df, location, gender) %>% select(bp.1s, location, gender)

#на сайте с темами есть пример данных и функция facet_grid для создания двух половин. В примере Automatic и manual. А у нас то же самое для location
bp.1s_plot <- ggplot(data = bp.1s_by_location_gender_df, mapping = aes(x = gender, y = bp.1s)) + facet_grid(. ~ location)
bp.1s_plot + theme_classic() + geom_violin() + geom_boxplot(width = 0.3, notch = TRUE)

#10.5
chol_age_gender_and_frame_df <- group_by(diabetes_df, gender, frame) %>% select(chol, age, gender, frame)
chol_age_gender_and_frame_plot <- ggplot(data = chol_age_gender_and_frame_df, mapping = aes(x = age, y = chol, color = frame, shape = gender))
chol_age_gender_and_frame_plot <- chol_age_gender_and_frame_plot + geom_point() + theme_classic()
chol_age_gender_and_frame_plot <- chol_age_gender_and_frame_plot + labs(title = "Age related changes in the cholesterol level", x = "Age", y = "Cholesterol level")
chol_age_gender_and_frame_plot <- chol_age_gender_and_frame_plot + scale_color_manual(name="", values = c("red", "magenta", "blue", "gray"), labels = c("large", "medium", "small", "NA"))
chol_age_gender_and_frame_plot <- chol_age_gender_and_frame_plot + scale_shape(name="")
chol_age_gender_and_frame_plot