#5.1
chol_df <- read.table("SISG-Data-cholesterol.csv", head = T, skip=14, sep = ",")

#5.2
cor(x = chol_df$chol,y = chol_df[, c("BMI","TG")], method = "s")
cor.test(x = chol_df$chol, y = chol_df$BMI, method = "s")$p.value
cor.test(x = chol_df$chol, y = chol_df$TG, method = "s")$p.value

#5.3
plot(chol ~ BMI, data = chol_df)
plot(chol ~ TG, data = chol_df)

#5.4
chol_by_TG_model <- lm(chol ~ TG, data = chol_df)
summary(chol_by_TG_model)$coefficients

shapiro.test(chol_by_TG_model$residuals)$p.value
t.test(chol_by_TG_model$residuals)$p.value
cor.test(chol_by_TG_model$residuals, chol_df$TG, method = "s")$p.value


chol_by_BMI_model <- lm(chol ~ BMI, data = chol_df)
summary(chol_by_BMI_model)$coefficients

shapiro.test(chol_by_BMI_model$residuals)$p.value
t.test(chol_by_BMI_model$residuals)$p.value
cor.test(chol_by_BMI_model$residuals,chol_df$BMI,method ="s")$p.value

#5.5
chol_by_BMI_and_TG_model <- lm(chol ~ BMI + TG, data = chol_df)
summary(chol_by_BMI_and_TG_model)$coefficients

shapiro.test(chol_by_BMI_and_TG_model$residuals)$p.value
t.test(chol_by_BMI_and_TG_model$residuals)$p.value
cor.test(chol_by_BMI_and_TG_model$residuals,chol_df$BMI,method ="s")$p.value

summary(chol_by_BMI_and_TG_model)$r.squared
summary(chol_by_BMI_model)$r.squared
summary(chol_by_TG_model)$r.squared

#5.6
chol_by_BMI_TG_and_sex_model <- lm(chol ~ BMI + TG + sex, data = chol_df)
summary(chol_by_BMI_TG_and_sex_model)$coefficients
summary(chol_by_BMI_TG_and_sex_model)$r.squared

#5.7
shapiro.test(x = chol_df$chol[chol_df$sex == 0])$p.value
shapiro.test(x = chol_df$chol[chol_df$sex == 1])$p.value

var.test(chol ~ sex, data = chol_df)$p.value
t.test(chol ~ sex, data = chol_df, var.equl = T)
#p-value 7.128e-07 => имеются статистически значимые различия между мужчинами и женщинами

summary(lm(chol ~ sex, data = chol_df))$coefficients
