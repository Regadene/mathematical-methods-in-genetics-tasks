#7.1
qrt_groups <- read.table("sample_groups.txt", head = T)
qrt_data <- read.table("QRT_PCR_LR.txt", head = T)
grt <- merge(qrt_groups, qrt_data, by.x = "tissue", by.y = "Tissue")

grt_delta_UBS_df <- grt
grt_delta_UBS_df[, -c(1, 2)] <- grt[, -c(1, 2)] - grt$UBC
grt_delta_UBS_df$group <- grt_delta_UBS_df$group - 1

shapiro.test(grt_delta_UBS_df$CD44[grt_delta_UBS_df$group == 0])$p.value > 0.05
shapiro.test(grt_delta_UBS_df$CD44[grt_delta_UBS_df$group == 1])$p.value > 0.05

var.test(CD44 ~ group, data = grt_delta_UBS_df[ , c("group","CD44")])$p.value > 0.05
t.test(CD44 ~ group, data = grt_delta_UBS_df)$p.value > 0.05
#Экспрессия отличается в группах

#7.2
set.seed(1)
rows_test <- sample(x = 1 : nrow(grt_delta_UBS_df), size = 72)
qrt_test <- grt_delta_UBS_df[rows_test, ]
qrt_train <- grt_delta_UBS_df[-rows_test, ]

fit_AKR1C1 <- glm(group ~ AKR1C1, data = qrt_train, family = binomial())
summary(fit_AKR1C1)

library(pROC) 
#7.3
qrt_test$prob <- predict(fit_AKR1C1, newdata = qrt_test, type = "response")
g <- roc(group ~ prob, data = qrt_test)
g$auc
plot(g) 

#7.4
set.seed(1)
rows_test2 <- sample(x = 1 : nrow(grt_delta_UBS_df), size = 72)
qrt_test2 <- grt_delta_UBS_df[rows_test2, ]
qrt_train2 <- grt_delta_UBS_df[-rows_test2, ]

fit_AKR1C1_CPSG3<- glm(group ~ AKR1C1 + CPSG3, data = qrt_train2, family = binomial())
summary(fit_AKR1C1_CPSG3)

qrt_test2$prob <- predict(fit_AKR1C1_CPSG3, newdata = qrt_test2, type = "response")
g2 <- roc(group ~ prob, data = qrt_test2)
g2$auc
plot(g2) 
