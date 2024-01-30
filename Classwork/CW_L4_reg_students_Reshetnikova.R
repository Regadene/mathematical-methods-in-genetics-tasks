library(tidyr)
library(dplyr)

#4.1
qrt_df <- read.table(file = "QRT-PCR_LR_standarts.tsv", sep = "\t", head = T)

#4.2
plot(qrt_df$AHCY ~ qrt_df$Conc)

#4.3 
plot(qrt_df$AHCY ~ log10(qrt_df$Conc))

#4.4
AHCY_model <- lm(qrt_df$AHCY ~ log10(qrt_df$Conc))
coefficients(AHCY_model)

#4.5
AHCY_model_high_сoncentration <- lm(data = qrt_df[1:12, ], formula = AHCY ~ log10(Conc))
coefficients(AHCY_model_high_сoncentration)

#4.6
elastic_slope_function = function(x) { (10^(-1/x) - 1) * 100 }

qrt_slopes <- pivot_longer(data = qrt_df, cols = -c("Conc")) %>%
  group_by(name) %>%
  summarise(slope = lm(formula = value ~ log10(Conc))$coefficients[2]) %>%
  mutate(elastic = elastic_slope_function(slope))

hist(x = qrt_slopes$slope)
hist(x = qrt_slopes$elastic)
#Наибольший к-т для гена ALUsq slope: -1.210085, elastic: 570.48%

#4.7
ALUsq_model <- lm(ALUsq ~ log10(Conc), data = qrt_df)

shapiro.test(x = ALUsq_model$residuals)
t.test(x = ALUsq_model$residuals)
cor.test(x = ALUsq_model$residuals, y = log10(qrt_df$Conc))

new_data <- data.frame(Conc = seq(1, 2e+05, 100))
CI <- predict(object = ALUsq_model, new = data.frame(Conc = seq(1, 2e+05, 100)), interval = "confidence", level = 0.95)

plot(ALUsq ~ log10(Conc), data = qrt_df)
abline(reg = ALUsq_model, col = 2)
lines(x = log10(new_data$Conc), y = CI[, 2], col = 2, lty = 2)
lines(x=log10(new_data$Conc), y = CI[, 3], col = 2, lty = 2)

#4.8
SDHA_model <- lm(SDHA ~ log10(Conc), data = qrt_df)
b1 <- summary(SDHA_model)$coefficients[2, 1]
b1_e <- summary(SDHA_model)$coefficients[2, 2]
b2 <- summary(AHCY_model)$coefficients[2, 1]
b2_e <- summary(AHCY_model)$coefficients[2, 2]
t <- (b1 - b2)/sqrt(b1_e^2 + b2_e^2)
t

pt(q = abs(t), df = 30 - 2, lower.tail = F) * 2
