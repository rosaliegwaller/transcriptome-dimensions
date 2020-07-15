#Predict survival or time to treatment failure
#Using the Cox Model and elastic-net penalty

#Load packages
library("glmnet")
library("survival")
library("selectiveInference")

#Load data
setwd("/Users/rosal/OneDrive - University of Utah/2020/analyze/data/transcriptome-dimensions/")
load(file='rdata/process_20200428.RData') #pcs and clinical data in dt
data <- merge(samples[,c("sample_id","patient_id","ttcos","censos","ttcpfs","censpfs","ttctf1","censtf1")],
              score_sd)
rm(cbat.dt,pc,samples,score_sd)

#Overall survival
Y <- Surv(time = data$ttcos, event = data$censos)
X <- data.matrix(dplyr::select(data,starts_with('PC')))

fit <- glmnet(X,Y,family = "cox")
cvfit <- cv.glmnet(X,Y,family = "cox",nfolds = 12)
beta_hat = as.numeric(coef(fit,x=X,y=Y,s=cvfit$lambda.min,extract=TRUE))
out = fixedLassoInf(x=X,y=data$ttcos,status=data$censos,family="cox",beta=beta_hat,lambda=cvfit$lambda.min)
coxph(Y~X)

plot(fit,label=TRUE)
print(fit)

cvfit <- cv.glmnet(X,Y,family = "cox",nfolds = 12)
coef(cvfit, s = "lambda.min")


plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se

coef.min = coef(cvfit, s = "lambda.min")
coef.min
active.min = which(coef.min != 0)
active.min
index.min = coef.min[active.min]
index.min
