###Script for Processing NAM Phenotype Data to get Mean Colony Counts, Heritability, Binary and Normaliation Transformations##
#Author: Patrick Woods#
#Date: January 25th, 2024#

library(lme4)
library(emmeans)
library(dplyr)
library(bestNormalize)

setwd("../Projects/NAM_project/phenotype_data/")
#####NAM POPULATION 1 (W2021006)#####
#Read in data and convert the variables to their correct types (such as factor for Plant ID)
nam1 <- read.csv("W2021006/W2021006_NAM_pheno_replaced_nodot.csv", header=T)
head(nam1)
str(nam1)
nam1$Plant_ID <- as.factor(nam1$Plant_ID)
nam1$Bench <- as.factor(nam1$Bench)
nam1$NAM_Family <- as.factor(nam1$NAM_Family)

#Split data by Bench to estimate means of Replicate plants
nam1_split <- split(nam1, nam1$Bench)

#Bench A#
nam1_A <- nam1_split$A
nam1_A_bin1_means <- aggregate(Bine.1.Colony.counts ~ Plant_ID, data = nam1_A, mean)
nam1_A_bin1_means$means <- nam1_A_bin1_means$Bine.1.Colony.counts 
nam1_A_bin1_means$Bine.1.Colony.counts <- NULL
nam1_A_bin2_means <- aggregate(Bine.2.Colony.counts ~ Plant_ID, data = nam1_A, mean)
nam1_A_bin2_means$means <- nam1_A_bin2_means$Bine.2.Colony.counts
nam1_A_bin2_means$Bine.2.Colony.counts <- NULL
nam1_A_bin1_2_comb <- rbind(nam1_A_bin1_means, nam1_A_bin2_means) #make a combined data frame that has a single column of means for Bine1 and Bine 2 on Bench A
nam1_A_overall_means <- aggregate(means ~ Plant_ID, data = nam1_A_bin1_2_comb, mean) #calculate overall mean for bench A replicates

#Bench B#
nam1_B <- nam1_split$B
nam1_B_bin1_means <- aggregate(Bine.1.Colony.counts ~ Plant_ID, data = nam1_B, mean)
nam1_B_bin1_means$means <- nam1_B_bin1_means$Bine.1.Colony.counts 
nam1_B_bin1_means$Bine.1.Colony.counts <- NULL
nam1_B_bin2_means <- aggregate(Bine.2.Colony.counts ~ Plant_ID, data = nam1_B, mean)
nam1_B_bin2_means$means <- nam1_B_bin2_means$Bine.2.Colony.counts
nam1_B_bin2_means$Bine.2.Colony.counts <- NULL
nam1_B_bin1_2_comb <- rbind(nam1_B_bin1_means, nam1_B_bin2_means) #make a combined data frame that has a single column of means for Bine1 and Bine 2 on Bench B
nam1_B_overall_means <- aggregate(means ~ Plant_ID, data = nam1_B_bin1_2_comb, mean) #calculate overall mean for bench B replicates

#Combine overall Bench means into a Single Dataframe to calculate broad sense heritability
nam1_A_B_means_comb <- rbind(nam1_A_overall_means, nam1_B_overall_means)
nam1_heritability_model <- lmer(means ~ (1|Plant_ID), data = nam1_A_B_means_comb)
summary(nam1_heritability_model)
nam1_total_variance <- 80.00 + 29.79
nam1_heritability <- 80.00/nam1_total_variance * 100
nam1_heritability #72.87% heritable

#Calculate final mean between Bench A and B
nam1_final_means <- aggregate(means ~ Plant_ID, data = nam1_A_B_means_comb, mean)

#Convert Mean Values less than 1 to "0", and Mean Values Greater than or Equal to 1 to "1".
nam1_final_means_bin <- nam1_final_means %>% 
                        mutate(bin_means = case_when(means >=1 ~ 1,
                                                    means < 1 ~ 0))


nam1_table <- table(nam1_final_means_bin$bin_means)
prop.table(nam1_table) #Nam Family 1 (W2021006) shows a 3:1 (Susceptible:Resistant) segregation which suggests a single locus segregating at 1:2:1 where the resistant allele is recessive.
prop.test(nam1_table) #significant p-value indicates that the 3:1 segregation we observe is real
hist(nam1_final_means_bin$bin_means)

#Now Need to Quantile Normalize Phenotypes of 1 or Greater to get a normal distribution for Susceptibility
#Add a new column with Simple Means previously Caclculated and Set Means Less Than 1 to "NA" so they are not included in the QN Transformation
nam1_final_means_bin$QN_means <- nam1_final_means_bin$means
nam1_final_means_bin$QN_means[nam1_final_means_bin$bin_means < 1] <- NA

#Quantile Normalize The QN_means column
nam1_means_qn <- orderNorm(nam1_final_means_bin$QN_means)
hist(nam1_means_qn$x.t)
shapiro.test(nam1_means_qn$x.t) #QN Transformation Successful!
View(data.frame(nam1_means_qn$x.t))

#Add QN Means Back to the Original Data File
nam1_final_means_bin$quantile_norm <- nam1_means_qn$x.t
cor.test(nam1_final_means_bin$QN_means, nam1_final_means_bin$quantile_norm, method = "spearman")

#Make Histogram Plots for Each Phenotype
hist(nam1_final_means_bin$means, main = "Histogram of W2021006 Mean Colony Counts (raw)", xlab = "Colony Count Score")
hist(nam1_final_means_bin$bin_means, main = "Histogram of W2021006 Mean Colony Counts (binary)", xlab = "Colony Count Score")
hist(nam1_final_means_bin$quantile_norm, main = "Histogram of W2021006 Mean Colony Counts (normalized)", xlab = "Colony Count Score")

#Write the final phenotype object out to a .csv file
write.csv(nam1_final_means_bin, file = "W2021006_phenotype_means_raw_binary_qn.csv", quote = F, row.names = F)

#####NAM POPULATION 2 (W2021009)#####
#Read in data and convert the variables to their correct types (such as factor for Plant ID)
nam2 <- read.csv("W2021009/W2021009_NAM_pheno_replaced_nodot.csv", header=T)
head(nam2)
str(nam2)
nam2$Plant_ID <- as.factor(nam2$Plant_ID)
nam2$Bench <- as.factor(nam2$Bench)
nam2$NAM_Family <- as.factor(nam2$NAM_Family)

#Split data by Bench to estimate means of Replicate plants
nam2_split <- split(nam2, nam2$Bench)

#Bench A#
nam2_A <- nam2_split$A
nam2_A_bin1_means <- aggregate(Bine.1.Colony.counts ~ Plant_ID, data = nam2_A, mean)
nam2_A_bin1_means$means <- nam2_A_bin1_means$Bine.1.Colony.counts 
nam2_A_bin1_means$Bine.1.Colony.counts <- NULL
nam2_A_bin2_means <- aggregate(Bine.2.Colony.counts ~ Plant_ID, data = nam2_A, mean)
nam2_A_bin2_means$means <- nam2_A_bin2_means$Bine.2.Colony.counts
nam2_A_bin2_means$Bine.2.Colony.counts <- NULL
nam2_A_bin1_2_comb <- rbind(nam2_A_bin1_means, nam2_A_bin2_means) #make a combined data frame that has a single column of means for Bine1 and Bine 2 on Bench A
nam2_A_overall_means <- aggregate(means ~ Plant_ID, data = nam2_A_bin1_2_comb, mean) #calculate overall mean for bench A replicates

#Bench B#
nam2_B <- nam2_split$B
nam2_B_bin1_means <- aggregate(Bine.1.Colony.counts ~ Plant_ID, data = nam2_B, mean)
nam2_B_bin1_means$means <- nam2_B_bin1_means$Bine.1.Colony.counts 
nam2_B_bin1_means$Bine.1.Colony.counts <- NULL
nam2_B_bin2_means <- aggregate(Bine.2.Colony.counts ~ Plant_ID, data = nam2_B, mean)
nam2_B_bin2_means$means <- nam2_B_bin2_means$Bine.2.Colony.counts
nam2_B_bin2_means$Bine.2.Colony.counts <- NULL
nam2_B_bin1_2_comb <- rbind(nam2_B_bin1_means, nam2_B_bin2_means) #make a combined data frame that has a single column of means for Bine1 and Bine 2 on Bench B
nam2_B_overall_means <- aggregate(means ~ Plant_ID, data = nam2_B_bin1_2_comb, mean) #calculate overall mean for bench B replicates

#Combine overall Bench means into a Single Dataframe to calculate broad sense heritability
nam2_A_B_means_comb <- rbind(nam2_A_overall_means, nam2_B_overall_means)
nam2_heritability_model <- lmer(means ~ (1|Plant_ID), data = nam2_A_B_means_comb)
summary(nam2_heritability_model)
nam2_total_variance <- 51.00 + 54.07
nam2_heritability <- 51.00/nam2_total_variance * 100
nam2_heritability #48.54% heritable

#Calculate final mean between Bench A and B
nam2_final_means <- aggregate(means ~ Plant_ID, data = nam2_A_B_means_comb, mean)

#Convert Mean Values less than 1 to "0", and Mean Values Greater than or Equal to 1 to "1".
nam2_final_means_bin <- nam2_final_means %>% 
  mutate(bin_means = case_when(means >=1 ~ 1,
                               means < 1 ~ 0))


nam2_table <- table(nam2_final_means_bin$bin_means)
prop.table(nam2_table)
prop.test(nam2_table) #close to a 1:1 segregation of Resistant:Susceptible, 44%:56%, non-significant p-value
hist(nam2_final_means_bin$bin_means)

#Now Need to Quantile Normalize Phenotypes of 1 or Greater to get a normal distribution for Susceptibility
#Add a new column with Simple Means previously Caclculated and Set Means Less Than 1 to "NA" so they are not included in the QN Transformation
nam2_final_means_bin$QN_means <- nam2_final_means_bin$means
nam2_final_means_bin$QN_means[nam2_final_means_bin$bin_means < 1] <- NA

#Quantile Normalize The QN_means column
nam2_means_qn <- orderNorm(nam2_final_means_bin$QN_means)
hist(nam2_means_qn$x.t)
shapiro.test(nam2_means_qn$x.t) #QN Transformation Successful!
View(data.frame(nam2_means_qn$x.t))

#Add QN Means Back to the Original Data File
nam2_final_means_bin$quantile_norm <- nam2_means_qn$x.t
cor.test(nam2_final_means_bin$QN_means, nam2_final_means_bin$quantile_norm, method = "spearman")

#Make Histogram Plots for Each Phenotype
hist(nam2_final_means_bin$means, main = "Histogram of W2021009 Mean Colony Counts (raw)", xlab = "Colony Count Score")
hist(nam2_final_means_bin$bin_means, main = "Histogram of W2021009 Mean Colony Counts (binary)", xlab = "Colony Count Score")
hist(nam2_final_means_bin$quantile_norm, main = "Histogram of W2021009 Mean Colony Counts (normalized)", xlab = "Colony Count Score")

#Write the final phenotype object out to a .csv file
write.csv(nam2_final_means_bin, file = "W2021009_phenotype_means_raw_binary_qn.csv", quote = F, row.names = F)


###NAM Population 3 (W2021010)###
#Read in data and convert the variables to their correct types (such as factor for Plant ID)
nam3 <- read.csv("W2021010/W2021010_NAM_pheno_replaced_nodot.csv", header=T)
head(nam3)
str(nam3)
nam3$Plant_ID <- as.factor(nam3$Plant_ID)
nam3$Bench <- as.factor(nam3$Bench)
nam3$NAM_Family <- as.factor(nam3$NAM_Family)

#Split data by Bench to estimate means of Replicate plants
nam3_split <- split(nam3, nam3$Bench)

#Bench A#
nam3_A <- nam3_split$A
nam3_A_bin1_means <- aggregate(Bine.1.Colony.counts ~ Plant_ID, data = nam3_A, mean)
nam3_A_bin1_means$means <- nam3_A_bin1_means$Bine.1.Colony.counts 
nam3_A_bin1_means$Bine.1.Colony.counts <- NULL
nam3_A_bin2_means <- aggregate(Bine.2.Colony.counts ~ Plant_ID, data = nam3_A, mean)
nam3_A_bin2_means$means <- nam3_A_bin2_means$Bine.2.Colony.counts
nam3_A_bin2_means$Bine.2.Colony.counts <- NULL
nam3_A_bin1_2_comb <- rbind(nam3_A_bin1_means, nam3_A_bin2_means) #make a combined data frame that has a single column of means for Bine1 and Bine 2 on Bench A
nam3_A_overall_means <- aggregate(means ~ Plant_ID, data = nam3_A_bin1_2_comb, mean) #calculate overall mean for bench A replicates

#Bench B#
nam3_B <- nam3_split$B
nam3_B_bin1_means <- aggregate(Bine.1.Colony.counts ~ Plant_ID, data = nam3_B, mean)
nam3_B_bin1_means$means <- nam3_B_bin1_means$Bine.1.Colony.counts 
nam3_B_bin1_means$Bine.1.Colony.counts <- NULL
nam3_B_bin2_means <- aggregate(Bine.2.Colony.counts ~ Plant_ID, data = nam3_B, mean)
nam3_B_bin2_means$means <- nam3_B_bin2_means$Bine.2.Colony.counts
nam3_B_bin2_means$Bine.2.Colony.counts <- NULL
nam3_B_bin1_2_comb <- rbind(nam3_B_bin1_means, nam3_B_bin2_means) #make a combined data frame that has a single column of means for Bine1 and Bine 2 on Bench B
nam3_B_overall_means <- aggregate(means ~ Plant_ID, data = nam3_B_bin1_2_comb, mean) #calculate overall mean for bench B replicates

#Combine overall Bench means into a Single Dataframe to calculate broad sense heritability
nam3_A_B_means_comb <- rbind(nam3_A_overall_means, nam3_B_overall_means)
nam3_heritability_model <- lmer(means ~ (1|Plant_ID), data = nam3_A_B_means_comb)
summary(nam3_heritability_model)
nam3_total_variance <- 38.825 + 9.843
nam3_heritability <- 38.825/nam3_total_variance * 100
nam3_heritability #79.77% heritable

#Calculate final mean between Bench A and B
nam3_final_means <- aggregate(means ~ Plant_ID, data = nam3_A_B_means_comb, mean)

#Convert Mean Values less than 1 to "0", and Mean Values Greater than or Equal to 1 to "1".
nam3_final_means_bin <- nam3_final_means %>% 
  mutate(bin_means = case_when(means >=1 ~ 1,
                               means < 1 ~ 0))


nam3_table <- table(nam3_final_means_bin$bin_means)
prop.table(nam3_table) #Close to a 1:1 segregation, 59% resistant: 41% susceptible
prop.test(nam3_table) #p-value = 0.02544 which indicates a slight departure from the 1:1 segregation but this could be due to sample size
hist(nam3_final_means_bin$bin_means)

#Now Need to Quantile Normalize Phenotypes of 1 or Greater to get a normal distribution for Susceptibility
#Add a new column with Simple Means previously Caclculated and Set Means Less Than 1 to "NA" so they are not included in the QN Transformation
nam3_final_means_bin$QN_means <- nam3_final_means_bin$means
nam3_final_means_bin$QN_means[nam3_final_means_bin$bin_means < 1] <- NA

#Quantile Normalize The QN_means column
nam3_means_qn <- orderNorm(nam3_final_means_bin$QN_means)
hist(nam3_means_qn$x.t)
shapiro.test(nam3_means_qn$x.t) #QN Transformation Successful!
View(data.frame(nam3_means_qn$x.t))

#Add QN Means Back to the Original Data File
nam3_final_means_bin$quantile_norm <- nam3_means_qn$x.t
cor.test(nam3_final_means_bin$QN_means, nam3_final_means_bin$quantile_norm, method = "spearman")

#Make Histogram Plots for Each Phenotype
hist(nam3_final_means_bin$means, main = "Histogram of W2021010 Mean Colony Counts (raw)", xlab = "Colony Count Score")
hist(nam3_final_means_bin$bin_means, main = "Histogram of W2021010 Mean Colony Counts (binary)", xlab = "Colony Count Score")
hist(nam3_final_means_bin$quantile_norm, main = "Histogram of W2021010 Mean Colony Counts (normalized)", xlab = "Colony Count Score")

#Write the final phenotype object out to a .csv file
write.csv(nam3_final_means_bin, file = "W2021010_phenotype_means_raw_binary_qn.csv", quote = F, row.names = F)

###NAM Population 4 (W2021011)###
#Read in data and convert the variables to their correct types (such as factor for Plant ID)
nam4 <- read.csv("W2021011/W2021011_NAM_pheno_replaced_nodot.csv", header=T)
head(nam4)
str(nam4)
nam4$Plant_ID <- as.factor(nam4$Plant_ID)
nam4$Bench <- as.factor(nam4$Bench)
nam4$NAM_Family <- as.factor(nam4$NAM_Family)

#Split data by Bench to estimate means of Replicate plants
nam4_split <- split(nam4, nam4$Bench)

#Bench A#
nam4_A <- nam4_split$A
nam4_A_bin1_means <- aggregate(Bine.1.Colony.counts ~ Plant_ID, data = nam4_A, mean)
nam4_A_bin1_means$means <- nam4_A_bin1_means$Bine.1.Colony.counts 
nam4_A_bin1_means$Bine.1.Colony.counts <- NULL
nam4_A_bin2_means <- aggregate(Bine.2.Colony.counts ~ Plant_ID, data = nam4_A, mean)
nam4_A_bin2_means$means <- nam4_A_bin2_means$Bine.2.Colony.counts
nam4_A_bin2_means$Bine.2.Colony.counts <- NULL
nam4_A_bin1_2_comb <- rbind(nam4_A_bin1_means, nam4_A_bin2_means) #make a combined data frame that has a single column of means for Bine1 and Bine 2 on Bench A
nam4_A_overall_means <- aggregate(means ~ Plant_ID, data = nam4_A_bin1_2_comb, mean) #calculate overall mean for bench A replicates

#Bench B#
nam4_B <- nam4_split$B
nam4_B_bin1_means <- aggregate(Bine.1.Colony.counts ~ Plant_ID, data = nam4_B, mean)
nam4_B_bin1_means$means <- nam4_B_bin1_means$Bine.1.Colony.counts 
nam4_B_bin1_means$Bine.1.Colony.counts <- NULL
nam4_B_bin2_means <- aggregate(Bine.2.Colony.counts ~ Plant_ID, data = nam4_B, mean)
nam4_B_bin2_means$means <- nam4_B_bin2_means$Bine.2.Colony.counts
nam4_B_bin2_means$Bine.2.Colony.counts <- NULL
nam4_B_bin1_2_comb <- rbind(nam4_B_bin1_means, nam4_B_bin2_means) #make a combined data frame that has a single column of means for Bine1 and Bine 2 on Bench B
nam4_B_overall_means <- aggregate(means ~ Plant_ID, data = nam4_B_bin1_2_comb, mean) #calculate overall mean for bench B replicates

#Combine overall Bench means into a Single Dataframe to calculate broad sense heritability
nam4_A_B_means_comb <- rbind(nam4_A_overall_means, nam4_B_overall_means)
nam4_heritability_model <- lmer(means ~ (1|Plant_ID), data = nam4_A_B_means_comb)
summary(nam4_heritability_model)
nam4_total_variance <- 43.18 + 13.46
nam4_heritability <- 43.18/nam4_total_variance * 100
nam4_heritability #76.24% heritable

#Calculate final mean between Bench A and B
nam4_final_means <- aggregate(means ~ Plant_ID, data = nam4_A_B_means_comb, mean)

#Convert Mean Values less than 1 to "0", and Mean Values Greater than or Equal to 1 to "1".
nam4_final_means_bin <- nam4_final_means %>% 
  mutate(bin_means = case_when(means >=1 ~ 1,
                               means < 1 ~ 0))


nam4_table <- table(nam4_final_means_bin$bin_means)
prop.table(nam4_table) #1:3 segregation of 29% resistant: 71% susceptible
prop.test(nam4_table) #very significant p-value which indicates strong departure from 1:1 segregation
hist(nam4_final_means_bin$bin_means)

#Now Need to Quantile Normalize Phenotypes of 1 or Greater to get a normal distribution for Susceptibility
#Add a new column with Simple Means previously Caclculated and Set Means Less Than 1 to "NA" so they are not included in the QN Transformation
nam4_final_means_bin$QN_means <- nam4_final_means_bin$means
nam4_final_means_bin$QN_means[nam4_final_means_bin$bin_means < 1] <- NA

#Quantile Normalize The QN_means column
nam4_means_qn <- orderNorm(nam4_final_means_bin$QN_means)
hist(nam4_means_qn$x.t)
shapiro.test(nam4_means_qn$x.t) #QN Transformation Successful!
View(data.frame(nam4_means_qn$x.t))

#Add QN Means Back to the Original Data File
nam4_final_means_bin$quantile_norm <- nam4_means_qn$x.t
cor.test(nam4_final_means_bin$QN_means, nam4_final_means_bin$quantile_norm, method = "spearman")

#Make Histogram Plots for Each Phenotype
hist(nam4_final_means_bin$means, main = "Histogram of W2021011 Mean Colony Counts (raw)", xlab = "Colony Count Score")
hist(nam4_final_means_bin$bin_means, main = "Histogram of W2021011 Mean Colony Counts (binary)", xlab = "Colony Count Score")
hist(nam4_final_means_bin$quantile_norm, main = "Histogram of W2021011 Mean Colony Counts (normalized)", xlab = "Colony Count Score")

#Write the final phenotype object out to a .csv file
write.csv(nam4_final_means_bin, file = "W2021011_phenotype_means_raw_binary_qn.csv", quote = F, row.names = F)

###Now Do a Combined Analysis that performs All above calculations but Treating All NAMs as one big family

nama <- rbind(nam1,nam2,nam3,nam4)
str(nama)

nama$Plant_ID <- as.factor(nama$Plant_ID)
nama$Bench <- as.factor(nama$Bench)
nama$NAM_Family <- as.factor(nama$NAM_Family)

#Split data by Bench to estimate means of Replicate plants
nama_split <- split(nama, nama$Bench)

#Bench A#
nama_A <- nama_split$A
nama_A_bin1_means <- aggregate(Bine.1.Colony.counts ~ Plant_ID, data = nama_A, mean)
nama_A_bin1_means$means <- nama_A_bin1_means$Bine.1.Colony.counts 
nama_A_bin1_means$Bine.1.Colony.counts <- NULL
nama_A_bin2_means <- aggregate(Bine.2.Colony.counts ~ Plant_ID, data = nama_A, mean)
nama_A_bin2_means$means <- nama_A_bin2_means$Bine.2.Colony.counts
nama_A_bin2_means$Bine.2.Colony.counts <- NULL
nama_A_bin1_2_comb <- rbind(nama_A_bin1_means, nama_A_bin2_means) #make a combined data frame that has a single column of means for Bine1 and Bine 2 on Bench A
nama_A_overall_means <- aggregate(means ~ Plant_ID, data = nama_A_bin1_2_comb, mean) #calculate overall mean for bench A replicates

#Bench B#
nama_B <- nama_split$B
nama_B_bin1_means <- aggregate(Bine.1.Colony.counts ~ Plant_ID, data = nama_B, mean)
nama_B_bin1_means$means <- nama_B_bin1_means$Bine.1.Colony.counts 
nama_B_bin1_means$Bine.1.Colony.counts <- NULL
nama_B_bin2_means <- aggregate(Bine.2.Colony.counts ~ Plant_ID, data = nama_B, mean)
nama_B_bin2_means$means <- nama_B_bin2_means$Bine.2.Colony.counts
nama_B_bin2_means$Bine.2.Colony.counts <- NULL
nama_B_bin1_2_comb <- rbind(nama_B_bin1_means, nama_B_bin2_means) #make a combined data frame that has a single column of means for Bine1 and Bine 2 on Bench B
nama_B_overall_means <- aggregate(means ~ Plant_ID, data = nama_B_bin1_2_comb, mean) #calculate overall mean for bench B replicates

#Combine overall Bench means into a Single Dataframe to calculate broad sense heritability
nama_A_B_means_comb <- rbind(nama_A_overall_means, nama_B_overall_means)
nama_heritability_model <- lmer(means ~ (1|Plant_ID), data = nama_A_B_means_comb)
summary(nama_heritability_model)
nama_total_variance <- 53.99 + 28.08
nama_heritability <- 53.99/nama_total_variance * 100
nama_heritability #65.78% heritable

#Calculate final mean between Bench A and B
nama_final_means <- aggregate(means ~ Plant_ID, data = nama_A_B_means_comb, mean)

#Convert Mean Values less than 1 to "0", and Mean Values Greater than or Equal to 1 to "1".
nama_final_means_bin <- nama_final_means %>% 
  mutate(bin_means = case_when(means >=1 ~ 1,
                               means < 1 ~ 0))


nama_table <- table(nama_final_means_bin$bin_means)
prop.table(nama_table) # close to a 1:1 segregation, 41% resistant, 59% susceptible, might indicate segregation distortion at R locus?
prop.test(nama_table) # highly significant p-value indicates a small yet real departure from a 1:1 segregation
hist(nama_final_means_bin$bin_means)

#Now Need to Quantile Normalize Phenotypes of 1 or Greater to get a normal distribution for Susceptibility
#Add a new column with Simple Means previously Caclculated and Set Means Less Than 1 to "NA" so they are not included in the QN Transformation
nama_final_means_bin$QN_means <- nama_final_means_bin$means
nama_final_means_bin$QN_means[nama_final_means_bin$bin_means < 1] <- NA

#Quantile Normalize The QN_means column
nama_means_qn <- orderNorm(nama_final_means_bin$QN_means)
hist(nama_means_qn$x.t)
shapiro.test(nama_means_qn$x.t) #QN Transformation Successful!
View(data.frame(nama_means_qn$x.t))

#Add QN Means Back to the Original Data File
nama_final_means_bin$quantile_norm <- nama_means_qn$x.t
cor.test(nama_final_means_bin$QN_means, nama_final_means_bin$quantile_norm, method = "spearman")

#Make Histogram Plots for Each Phenotype
hist(nama_final_means_bin$means, main = "Histogram of ALL NAM Mean Colony Counts (raw)", xlab = "Colony Count Score")
hist(nama_final_means_bin$bin_means, main = "Histogram of All NAM Mean Colony Counts (binary)", xlab = "Colony Count Score")
hist(nama_final_means_bin$quantile_norm, main = "Histogram of All NAM Mean Colony Counts (normalized)", xlab = "Colony Count Score")

#Write the final phenotype object out to a .csv file
write.csv(nama_final_means_bin, file = "All_NAM_phenotype_means_raw_binary_qn.csv", quote = F, row.names = F)

#Extract lines with Binary values = 1 to calculate heritability of susceptibility trait
nama_sub <- subset(nama_final_means_bin, nama_final_means_bin$bin_means == "1")
nama_left_join <- left_join(nama_sub, nama_A_B_means_comb, by = 'Plant_ID')
nama_susceptible_heritability_model <- lmer(means.y ~ (1|Plant_ID), data = nama_left_join)
summary(susceptible_heritability_model)
total_nama_susceptibility_variance <- 61.07 + 47.34
susceptible_heritability <- 61.07/total_nama_susceptibility_variance * 100
susceptible_heritability #56.33% heritable

#Write out a .csv file that contains the raw means for susceptible lines only
write.csv(nama_left_join, file = "NAM_all_suscpetible_line_means.csv")

#Modify in Excel to have a Family Column and then read data file back in
susher <- read.csv("All_NAM/NAM_all_suscpetible_line_means.csv", header=T)
str(susher)
susher$Family <- as.factor(susher$Family)
susher$Plant_ID <- as.factor(susher$Plant_ID)
susher_split <- split(susher, susher$Family)
nam1_sus <- susher_split$W2021006
nam2_sus <- susher_split$W2021009
nam3_sus <- susher_split$W2021010
nam4_sus <- susher_split$W2021011

#Now make 4 LMER models to estimate within family heritability of susceptibility

nam1_susmod <- lmer(raw_means.y ~ (1|Plant_ID), data = nam1_sus)
summary(nam1_susmod)
nam1_susmod_totalvar <- 90.96 + 41.34
nam1_sus_heritability <- 90.96/nam1_susmod_totalvar * 100
nam1_sus_heritability #68.75% heritable

nam2_susmod <- lmer(raw_means.y ~ (1|Plant_ID), data = nam2_sus)
summary(nam2_susmod)
nam2_susmod_totalvar <- 50.32 + 96.05
nam2_sus_heritability <- 50.32/nam2_susmod_totalvar * 100
nam2_sus_heritability #34.38% heritable

nam3_susmod <- lmer(raw_means.y ~ (1|Plant_ID), data = nam3_sus)
summary(nam3_susmod)
nam3_susmod_totalvar <- 53.31 + 24.21
nam3_sus_heritability <- 53.31/nam3_susmod_totalvar * 100
nam3_sus_heritability #68.77% heritable


nam4_susmod <- lmer(raw_means.y ~ (1|Plant_ID), data = nam4_sus)
summary(nam4_susmod)
nam4_susmod_totalvar <- 44.46 + 18.85
nam4_sus_heritability <- 44.46/nam4_susmod_totalvar * 100
nam4_sus_heritability #70.23% heritable

#Galena Mean

#Split data by Bench to estimate means of Replicate plants
gal <- read.csv(file.choose(),header=T)
gal_split <- split(gal, gal$Bench)

#Bench A#
gal_A <- gal_split$A
gal_A_bin1_means <- aggregate(Bine.1.Colony.counts ~ Plant_ID, data = gal_A, mean)
gal_A_bin1_means$means <- gal_A_bin1_means$Bine.1.Colony.counts 
gal_A_bin1_means$Bine.1.Colony.counts <- NULL
gal_A_bin2_means <- aggregate(Bine.2.Colony.counts ~ Plant_ID, data = gal_A, mean)
gal_A_bin2_means$means <- gal_A_bin2_means$Bine.2.Colony.counts
gal_A_bin2_means$Bine.2.Colony.counts <- NULL
gal_A_bin1_2_comb <- rbind(gal_A_bin1_means, gal_A_bin2_means) #make a combined data frame that has a single column of means for Bine1 and Bine 2 on Bench A
gal_A_overall_means <- aggregate(means ~ Plant_ID, data = gal_A_bin1_2_comb, mean) #calculate overall mean for bench A replicates

#Bench B#
gal_B <- gal_split$B
gal_B_bin1_means <- aggregate(Bine.1.Colony.counts ~ Plant_ID, data = gal_B, mean)
gal_B_bin1_means$means <- gal_B_bin1_means$Bine.1.Colony.counts 
gal_B_bin1_means$Bine.1.Colony.counts <- NULL
gal_B_bin2_means <- aggregate(Bine.2.Colony.counts ~ Plant_ID, data = gal_B, mean)
gal_B_bin2_means$means <- gal_B_bin2_means$Bine.2.Colony.counts
gal_B_bin2_means$Bine.2.Colony.counts <- NULL
gal_B_bin1_2_comb <- rbind(gal_B_bin1_means, gal_B_bin2_means) #make a combined data frame that has a single column of means for Bine1 and Bine 2 on Bench B
gal_B_overall_means <- aggregate(means ~ Plant_ID, data = gal_B_bin1_2_comb, mean) #calculate overall mean for bench B replicates

#Combine overall Bench means into a Single Dataframe
gal_A_B_means_comb <- rbind(gal_A_overall_means, gal_B_overall_means)


#Calculate final mean between Bench A and B
gal_final_means <- aggregate(means ~ Plant_ID, data = gal_A_B_means_comb, mean)

