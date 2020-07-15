#Interpert regression results

# packages
library("data.table")
library("ggplot2")
library("jtools")
library("MASS")
library("sjPlot")
library("dplyr")

# data
setwd("/Users/rosal/OneDrive - University of Utah/2020/analyze/data/transcriptome-dimensions/")
load(file = "rdata/regression_models-20200429.RData") #regression results from "regression_overview_plot.R"
load(file='rdata/process_20200428.RData') #pcs and clinical data in dt
dt <- merge(samples,score_sd,by='sample_id')
rm(samples,score_sd,pc)
baseline <- dt[collection_reason=='Baseline']$sample_id

# UAMC 70 gene risk panel ----
base.cbat <- cbat.dt[sample_id%in%baseline]
anno_up <- base.cbat %>% dplyr::select("AIM2","LARS2","OPN3","CCT2","UBE2I","RFC4","CKS1B","PFN1","ILF3","IFI16","TBRG4","ENO1","EXOSC4",
                                     "TAGLN2","RUVBL1","ALDOA","CPSF3","LGALS1","RAD18","SNX5","PSMD4","RAN","CBX3","TMPO","TCOF1","YWHAZ")
anno_dw <- base.cbat %>% dplyr::select("GNG10","AHCYL1","EVI5","CTBS","UBE2R2","FUCA1","LTBP1","TRIM33")
# Down - 7 not in normalized data: 
# "PNPLA4","KIAA1754","MCLC","AD-020","PARG1","RFP2""FLJ20489"
# Up - 24 not in normalized data: 
# "FABP5","PDHA1","TRIP13","SELI","SLCI19A1","ASPM","STK6",
# "FLJ13052","LAS1L","BIRC5","CKAP1","MGC57827","DKFZp779O175","PAPD1",
# "EIF2C2","MGC4308","DSG2","C6orf173","MGC15606","KIF14","DKFZP586L0724",
# "WEE1","ROBO1","MPHOSPH1"

# compute geometric means
x <- base.cbat[,"sample_id"]
x$up <- anno_up %>% rowMeans()
x$dw <- anno_dw %>% rowMeans()

# compute proporition of average low / high, note, already in log2 scale
x$score <- x$up - x$dw
hist(x$score,breaks = 200,main = "UAMC 70 Gene Risk Score",xlab = "log2(mean up reg) - log2(mean down reg)")

unique(x$score == fit$model$score) #score is the same as stored in model

# model statistics
fit = models$`UAMC-70`
summary(fit)
f = summary(fit)$fstatistic
pf(f[1],f[2],f[3],lower.tail=F) #calcualte overall model p-value

actual_preds <- data.frame(cbind(actual=fit$model$score,predicted=fit$fitted.values))
cor(actual_preds)
mean(abs((actual_preds$predicted - actual_preds$actual))/actual_preds$actual) #mean absolute percentage error

summ(fit)

# plot fitted values x acutal values
theme_set(theme_minimal() + theme(legend.position = "right"))
ggplot(actual_preds,aes(y = actual, x = predicted)) +
  geom_point()

# clean up workspace
rm(actual_preds,anno_dw,anno_up,base.cbat,cbat.dt,dt,fit,mdt,x,baseline,f)

# cytogenetic features ----
summ(models$`d(17p)`)
summary(dt[,"D_TRI_CF_ABNORMALITYPR11"])

summ(models$`t(14;16)`)
pchisq(124.933,df = 40,lower.tail = F)
summary(dt[,"D_TRI_CF_ABNORMALITYPR8"])
summary(models$`t(14;16)`)

summ(models$`a(1q)`)
pchisq(256.52,df = 40,lower.tail = F)
summary(dt[,"D_TRI_CF_ABNORMALITYPR13"])
summary(models$`a(1q)`)

summ(models$`t(4;14)`)
pchisq(210.67,df = 40,lower.tail = F)
summary(dt[,"D_TRI_CF_ABNORMALITYPR3"])
summary(models$`t(4;14)`)

summ(models$`t(11;14)`)
pchisq(225.78,df = 40,lower.tail = F)
summary(dt[,"D_TRI_CF_ABNORMALITYPR6"])
summary(models$`t(11;14)`)

# ISS ----
models$ISS
summary(models$ISS,digits=3)
confint(models$ISS)
plot(models$ISS)

# Disease course ----
models$OS$summary

models$PFS$summary
models$PFS$summary$logtest

models$TFTF$summary
models$TFTF$summary$logtest

# Disparties ----
models$race
summ(models$race)
pchisq(130.31,df = 40,lower.tail = F)
summary(dt[,"D_PT_race"])
summary(models$race)
exp(coef(models$race))
exp(confint(models$race))

models$ethnic
summ(models$ethnic)
pchisq(69.76,df = 40,lower.tail = F)
summary(dt[,"D_PT_ethnic"])
summary(models$ethnic)
exp(coef(models$ethnic))[17]
exp(confint(models$ethnic))[17,]

models$gender
summ(models$gender)
pchisq(111.51,df = 40,lower.tail = F)
summary(dt[,"D_PT_gender"])
summary(models$gender)
exp(coef(models$gender))[40]
exp(confint(models$gender))[40,]

models$age
summ(models$age)
f = summary(models$age)$fstatistic
pf(f[1],f[2],f[3],lower.tail=F)
summary(dt[,"D_PT_age"])
summary(models$age)
exp(coef(models$age))[c(2,22)]
exp(confint(models$age))[c(2,22),]
