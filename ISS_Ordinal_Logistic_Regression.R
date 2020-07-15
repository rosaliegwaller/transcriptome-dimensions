#Install relevant packages
install.packages("MASS")
            
#Load package
library("MASS")

#Read data file of pc scores and diagnosis ISS  
setwd("/Users/rosal/OneDrive - University of Utah/2020/analyze/data/transcriptome-dimensions/")
load(file='rdata/process_20200422.RData') #pcs and clinical data in dt
data <- clin.pc.dt %>% dplyr::select("D_PT_iss",ends_with("_sd")) %>% na.omit()
head(data)

#Ordering the dependent variable
data$D_PT_iss = factor(data$D_PT_iss, levels = c(1,2,3), ordered = TRUE) 

#Summarizing the data
summary(data)

#Build ordinal logistic regression model on all data
model = polr(D_PT_iss ~ ., data = data, Hess = TRUE)
summary(model)
Anova(model)

summary_table <- coef(summary(model))
pval <- pnorm(abs(summary_table[, "t value"]),lower.tail = FALSE)*2
summary_table <- cbind(summary_table, "p value" = pval)
summary_table

#Plotting the effects 
Effect(focal.predictors = "PC3_sd", model)
plot(Effect(focal.predictors = "PC3_sd", model))
plot(Effect(focal.predictors = c("PC1_sd", "PC2_sd"),model))


#Build ordinal logistic regression model on training data
#Random sampling 
samplesize = 0.60*nrow(data)
set.seed(100)
index = sample(seq_len(nrow(data)), size = samplesize)

#Creating training and test set 
datatrain = data[index,]
datatest = data[-index,]

#Build ordinal logistic regression model
modeltrain = polr(D_PT_iss ~ ., data = datatrain, Hess = TRUE)
summary(modeltrain)
Anova(modeltrain)

summary_table <- coef(summary(modeltrain))
pval <- pnorm(abs(summary_table[, "t value"]),lower.tail = FALSE)*2
summary_table <- cbind(summary_table, "p value" = pval)
summary_table

#Compute confusion table and misclassification error
predict_iss = predict(modeltrain,datatest)
table(datatest$D_PT_iss, predict_iss)
mean(as.character(datatest$D_PT_iss) != as.character(predict_iss))

#Plotting the effects 
Effect(focal.predictors = "PC3_sd", modeltrain)
plot(Effect(focal.predictors = "PC3_sd", modeltrain))
plot(Effect(focal.predictors = c("PC1_sd", "PC2_sd"),modeltrain))




