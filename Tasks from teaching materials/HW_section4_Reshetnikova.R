#4.1
test_numeric <- 1  #Создание numeric переменной
test_integer <- 2L #Создание integer переменной
test_character <- 'Yki'#Создание character переменной
test_logical <- F #Создание logical переменной

#4.2
defaultWD <- getwd() #присвоило такой путь "C:/Users/Rimalon/Desktop"
setwd("C:/Users/Rimalon/Desktop/R_Project_Resh") 
currentWD <- getwd()

#4.3
load(file = 'variables_2nd_week')
N_variables <- 3L

#4.4
str(current_year)
str(desired_grade_for_the_cource)
str(true_string)
N_int <- 1
N_char <- 1
N_logi <- 0
N_num <- 1

#4.5
save(currentWD, file = 'currentWD')
rm(currentWD)

#4.6
result_of_division <- current_year/10