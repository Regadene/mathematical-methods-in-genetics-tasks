#5.1
default_logical <- vector(length = 1)
default_integer <- vector(mode = 'integer', length = 1) 
default_numeric <- vector(mode = 'numeric', length = 1) 
default_character <- vector(mode = 'character', length = 1)

#5.2
x <- 1:6
y <- 10:4
integer_vector <- c(x,y)

numeric_vector_1<-seq(from = 1, to = 20, by = 0.5)

numeric_vector_2<-numeric_vector_1*c(1,1,0)

#5.3
integer_vector !=4 & integer_vector !=5

#5.4
load(file = "variables_3rd_week")
logi_vector_measur <- measurements >=10 & measurements <=40
table(logi_vector_measur) #616 TRUE

#5.5
table(gender) 
#     F   Female   M     Male     ND 
#    537      9    341    118    340
genders <- names(table(gender)) 
most_frequent_gender <- 'F'

#5.6
table(day_temperature)
# Cold Freezy    Hot Normal   Warm 
# 214     19    581    724   1107 
good_days <- 1107+581

#5.7
factor_temp_ordered <- factor(day_temperature, levels = c('Freezy','Cold','Normal','Warm','Hot'), ordered = TRUE)

#5.8
first_matrix=matrix(data = rep(c(0,2,0,4,0,6,0,8,0,10,0),10),ncol=10,nrow=10)

#5.9
Lengths<- list(integer_vector = length(integer_vector),numeric_vector_1 = length(numeric_vector_1),good_days=length(good_days),day_temperature=length(day_temperature), first_matrix = length(first_matrix))

#5.10
separator <- "\t"
NLines <- 10
diabetes_df <- read.table(file='diabetes.txt',sep = separator, skip = NLines, header = TRUE)
readLines('diabetes.txt', n = 12)

#5.11
Number_of_patients <- nrow(diabetes_df)

#5.12
Pound_to_kg <- function(Pound){
      Pound/2.205
   }