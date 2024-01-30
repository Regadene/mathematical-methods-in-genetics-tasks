#6.1
noNA <- complete.cases(diabetes_df)
diabetes_df[noNA, ]

#6.2
MaleAndFemale <- diabetes_df['gender']
F_amount<-length(MaleAndFemale[MaleAndFemale == 'female'])
M_amount<-length(MaleAndFemale[MaleAndFemale == 'male'])
MF_ratio <- M_amount/F_amount
MF_ratio_char <- paste(as.character(M_amount),':', as.character(F_amount)) #чтобы было строчное отображение можно использовать конкатенацию строк с использованием функции "paste".

#6.3
Height_And_Gender_df <- diabetes_df[c('height','gender')]
noNA_H_a_G <- complete.cases(Height_And_Gender_df)
noNA_H_a_G_df <- Height_And_Gender_df[noNA_H_a_G, ]

male_70plus <- (noNA_H_a_G_df$gender == 'male' & noNA_H_a_G_df$height > 70)
length(male_70plus[male_70plus == TRUE])

female_70plus <- (noNA_H_a_G_df$gender == 'female' & noNA_H_a_G_df$height > 70)
length(female_70plus[female_70plus == TRUE])

#6.4
noNA_frames <- complete.cases(diabetes_df['frame'])
noNA_frames_df<-diabetes_df[noNA_frames, ]
new_table <- noNA_frames_df[c('id', 'gender', 'age', 'height', 'weight', 'location')]
str(new_table)

#6.5
new_table$height <- round(new_table$height * 2.54)
new_table$weight <- round(new_table$weight / 2.205)
names(new_table)[4] <- "height_cm"
names(new_table)[5] <- "weight_kg"

#6.6
new_table$gender[new_table$gender == "female"] <- "F"
new_table$gender[new_table$gender == "male"] <- "M"

#6.7
Bucki_vector_ids <- new_table$id[new_table$location == 'Buckingham']
Bucki_frame <- data.frame(id=Bucki_vector_ids)
Bucki_threatment_placebo_vector <- rep(c('treatment', 'placebo'), length(Bucki_vector_ids)/2+1)
length(Bucki_threatment_placebo_vector) <- length(Bucki_vector_ids)
Bucki_frame$group <- Bucki_threatment_placebo_vector

Louisa_vector_ids <- new_table$id[new_table$location == 'Louisa']
Louisa_frame <- data.frame(id=Louisa_vector_ids)
Louisa_threatment_placebo_vector <- rep(c('treatment', 'placebo'), length(Louisa_vector_ids)/2+1)
length(Louisa_threatment_placebo_vector) <- length(Louisa_vector_ids)
Louisa_frame$group <- Louisa_threatment_placebo_vector

patient_groups <- rbind(Bucki_frame, Louisa_frame)

#6.8
new_table_w_groups <- merge(x=patient_groups, y=new_table, by="id")




