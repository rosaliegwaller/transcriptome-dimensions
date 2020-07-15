#Regression and plot

# packages
library(dplyr)
library(MASS)
library(survivalAnalysis)
library(data.table)
library(ggplot2)

# data
setwd("/Users/rosal/OneDrive - University of Utah/2020/analyze/data/transcriptome-dimensions/")
load(file='rdata/process_20200428.RData') #pcs and clinical data in dt
dt <- merge(samples,score_sd,by='sample_id')
rm(samples,score_sd,pc)

# Shaughnessy 70-gene risk score
up <- c("FABP5","PDHA1","TRIP13","AIM2","SELI","SLCI19A1","LARS2","OPN3","ASPM","CCT2","UBE2I","STK6","FLJ13052",
        "LAS1L","BIRC5","RFC4","CKS1B","CKAP1","MGC57827","DKFZp779O175","PFN1","ILF3","IFI16","TBRG4","PAPD1",
        "EIF2C2","MGC4308","ENO1","DSG2","C6orf173","EXOSC4","TAGLN2","RUVBL1","ALDOA","CPSF3","MGC15606","LGALS1",
        "RAD18","SNX5","PSMD4","RAN","KIF14","CBX3","TMPO","DKFZP586L0724","WEE1","ROBO1","TCOF1","YWHAZ","MPHOSPH1")

down <- c("GNG10","PNPLA4","KIAA1754","AHCYL1","MCLC","EVI5","AD-020","PARG1","CTBS","UBE2R2","FUCA1","RFP2","FLJ20489","LTBP1","TRIM33")

anno_up <- cbat.dt %>% dplyr::select("AIM2","LARS2","OPN3","CCT2","UBE2I",
                                     "RFC4","CKS1B","PFN1","ILF3","IFI16","TBRG4",
                                     "ENO1","EXOSC4","TAGLN2","RUVBL1","ALDOA","CPSF3","LGALS1",
                                     "RAD18","SNX5","PSMD4","RAN","CBX3","TMPO","TCOF1","YWHAZ")

anno_dw <- cbat.dt %>% dplyr::select("GNG10","AHCYL1","EVI5","CTBS","UBE2R2","FUCA1","LTBP1","TRIM33")

# compute geometric means
x <- cbat.dt[,"sample_id"]
x$up <- anno_up %>% rowMeans()
x$dw <- anno_dw %>% rowMeans()

# compute proporition of average low / high, note, already in log2 scale
x$score <- x$up - x$dw

# Not in data down: "PNPLA4","KIAA1754","MCLC","AD-020","PARG1","RFP2""FLJ20489"
# Not in data up: "FABP5","PDHA1","TRIP13","SELI","SLCI19A1","ASPM","STK6","FLJ13052","LAS1L","BIRC5","CKAP1","MGC57827","DKFZp779O175","PAPD1","EIF2C2","MGC4308","DSG2","C6orf173","MGC15606","KIF14","DKFZP586L0724","WEE1","ROBO1","MPHOSPH1"

##Run linear regression
mdt <- merge(x,dt,by="sample_id") %>% dplyr::select("score",ends_with("_sd")) %>% na.omit()
fit <- lm(data = mdt, formula = score ~ .)

beta <- coef(fit)
RR <- exp(coef(fit))
CI <- exp(confint(fit))
colnames(CI) <- c("Lower", "Upper")
P <- coef(summary(fit))[,4]
tb <- as.data.table(cbind(beta,RR, CI, P)[-1,])
tb$Y <- "UAMC-70"
tb$PC <- 1:40

sum.dt <- tb
models <- list("UAMC-70"=fit)
rm(tb,anno_dw,anno_up,cbat.dt,x,down,up,mdt,fit,CI,P,RR,beta)

# FISH ----
##D_TRI_CF_ABNORMALITYPR6: t(11;14) abnormality present standard risk
mdt <- dt %>% dplyr::select("D_TRI_CF_ABNORMALITYPR6",ends_with("_sd")) %>% na.omit()
fit <- glm(data = mdt, formula = D_TRI_CF_ABNORMALITYPR6 ~ ., family = "binomial")

beta <- coef(fit)
RR <- exp(coef(fit))
CI <- exp(confint(fit))
colnames(CI) <- c("Lower", "Upper")
P <- coef(summary(fit))[,4]
tb <- as.data.table(cbind(beta,RR, CI, P)[-1,])
tb$Y <- "t(11;14)"
tb$PC <- 1:40

sum.dt <- rbind(tb,sum.dt)
models <- c(list("t(11;14)"=fit),models)
rm(mdt,fit,RR,CI,P,tb,beta)

##D_TRI_CF_ABNORMALITYPR4: t(6;14) abnormality present standard risk
mdt <- dt %>% dplyr::select("D_TRI_CF_ABNORMALITYPR4",ends_with("_sd")) %>% na.omit()
fit <- glm(data = mdt, formula = D_TRI_CF_ABNORMALITYPR4 ~ ., family = "binomial") #did not converge
#beta <- coef(fit)
#RR <- exp(coef(fit))
#CI <- exp(confint(fit))
#colnames(CI) <- c("Lower", "Upper")
#P <- coef(summary(fit))[,4]
#tb <- as.data.table(cbind(beta,RR, CI, P)[-1,])
#tb$Y <- "t(6;14)"
#tb$PC <- 1:40

#sum.dt <- rbind(tb,sum.dt)
#models <- c(list("t(6;14)"=fit),models))
rm(mdt,fit,RR,CI,P,tb,beta)

##D_TRI_CF_ABNORMALITYPR3: t(4;14) abnormality present intermediate risk
mdt <- dt %>% dplyr::select("D_TRI_CF_ABNORMALITYPR3",ends_with("_sd")) %>% na.omit()
fit <- glm(data = mdt, formula = D_TRI_CF_ABNORMALITYPR3 ~ ., family = "binomial")
beta <- coef(fit)
RR <- exp(coef(fit))
CI <- exp(confint(fit))
colnames(CI) <- c("Lower", "Upper")
P <- coef(summary(fit))[,4]
tb <- as.data.table(cbind(beta,RR, CI, P)[-1,])
tb$Y <- "t(4;14)"
tb$PC <- 1:40

sum.dt <- rbind(tb,sum.dt)
models <- c(list("t(4;14)"=fit),models)
rm(mdt,fit,RR,CI,P,tb,beta)

##D_TRI_CF_ABNORMALITYPR13: 1q amplification abnormality present intermediate risk
mdt <- dt %>% dplyr::select("D_TRI_CF_ABNORMALITYPR13",ends_with("_sd")) %>% na.omit()
fit <- glm(data = mdt, formula = D_TRI_CF_ABNORMALITYPR13 ~ ., family = "binomial")
beta <- coef(fit)
RR <- exp(coef(fit))
CI <- exp(confint(fit))
colnames(CI) <- c("Lower", "Upper")
P <- coef(summary(fit))[,4]
tb <- as.data.table(cbind(beta,RR, CI, P)[-1,])
tb$Y <- "a(1q)"
tb$PC <- 1:40

sum.dt <- rbind(tb,sum.dt)
models <- c(list("a(1q)"=fit),models)
rm(mdt,fit,RR,CI,P,tb,beta)

##D_TRI_CF_ABNORMALITYPR8: t(14;16) abnormality present high risk
mdt <- dt %>% dplyr::select("D_TRI_CF_ABNORMALITYPR8",ends_with("_sd")) %>% na.omit()
fit <- glm(data = mdt, formula = D_TRI_CF_ABNORMALITYPR8 ~ ., family = "binomial")
beta <- coef(fit)
RR <- exp(coef(fit))
CI <- exp(confint(fit))
colnames(CI) <- c("Lower", "Upper")
P <- coef(summary(fit))[,4]
tb <- as.data.table(cbind(beta,RR, CI, P)[-1,])
tb$Y <- "t(14;16)"
tb$PC <- 1:40

sum.dt <- rbind(tb,sum.dt)
models <- c(list("t(14;16)"=fit),models)
rm(mdt,fit,RR,CI,P,tb,beta)

##D_TRI_CF_ABNORMALITYPR9: t(14;20) abnormality present high risk
mdt <- dt %>% dplyr::select("D_TRI_CF_ABNORMALITYPR9",ends_with("_sd")) %>% na.omit()
fit <- glm(data = mdt, formula = D_TRI_CF_ABNORMALITYPR9 ~ ., family = "binomial") #fitted probabilities numerically 0 or 1 occurred 
#beta <- coef(fit)
#RR <- exp(coef(fit))
#CI <- exp(confint(fit))
#colnames(CI) <- c("Lower", "Upper")
#P <- coef(summary(fit))[,4]
#tb <- as.data.table(cbind(beta,RR, CI, P)[-1,])
#tb$Y <- "t(14;20)"
#tb$PC <- 1:40

#sum.dt <- rbind(tb,sum.dt)
#models <- list("t(14;20)"=fit,models)
rm(mdt,fit,RR,CI,P,tb,beta)

##D_TRI_CF_ABNORMALITYPR11: del 17p abnormality present high risk
mdt <- dt %>% dplyr::select("D_TRI_CF_ABNORMALITYPR11",ends_with("_sd")) %>% na.omit()
fit <- glm(data = mdt, formula = D_TRI_CF_ABNORMALITYPR11 ~ ., family = "binomial")
beta <- coef(fit)
RR <- exp(coef(fit))
CI <- exp(confint(fit))
colnames(CI) <- c("Lower", "Upper")
P <- coef(summary(fit))[,4]
tb <- as.data.table(cbind(beta,RR, CI, P)[-1,])
tb$Y <- "d(17p)"
tb$PC <- 1:40

sum.dt <- rbind(tb,sum.dt)
models <- c(list("d(17p)"=fit),models)
rm(mdt,fit,RR,CI,P,tb,beta)

# International Staging System Score ----
#Ordering the dependent variable
dt$D_PT_iss = factor(dt$D_PT_iss, levels = c(1,2,3), ordered = TRUE)
mdt <- dt %>% dplyr::select("D_PT_iss",ends_with("_sd")) %>% na.omit()
#Build ordinal logistic regression model on all data
fit = polr(D_PT_iss ~ ., data = mdt, Hess = TRUE)
beta <- coef(fit)
RR <- exp(coef(fit))
CI <- exp(confint(fit))
colnames(CI) <- c("Lower", "Upper")
P <- pnorm(abs(coef(summary(fit))[-c(41,42),"t value"]),lower.tail = FALSE)*2
tb <- as.data.table(cbind(beta,RR, CI, P))
tb$Y <- "ISS"
tb$PC <- 1:40

sum.dt <- rbind(tb,sum.dt)
models <- c(list("ISS"=fit),models)
rm(mdt,fit,RR,CI,P,tb,beta)

# Survival/Time to treatment failure ----
##Overall survival
dt$censos <- as.integer(dt$censos)
os <- analyse_multivariate(data = dt,
                           time_status = vars(ttcos,censos),
                           covariates = colnames(dplyr::select(dt,ends_with("_sd"))))

fit <- os$coxph
beta <- coef(fit)
RR <- exp(coef(fit))
CI <- exp(confint(fit))
colnames(CI) <- c("Lower", "Upper")
P <- coef(summary(fit))[,5]
tb <- as.data.table(cbind(beta,RR, CI, P))
tb$Y <- "OS"
tb$PC <- 1:40

sum.dt <- rbind(tb,sum.dt)
models <- c(list("OS"=os),models)
rm(fit,RR,CI,P,tb,os,beta)

##Progression free survival
dt$censpfs <- as.integer(dt$censpfs)
pfs <- analyse_multivariate(data = dt,
                            time_status = vars(ttcpfs,censpfs),
                            covariates = colnames(dplyr::select(dt,ends_with("_sd"))))

fit <- pfs$coxph
beta <- coef(fit)
RR <- exp(coef(fit))
CI <- exp(confint(fit))
colnames(CI) <- c("Lower", "Upper")
P <- coef(summary(fit))[,5]
tb <- as.data.table(cbind(beta,RR, CI, P))
tb$Y <- "PFS"
tb$PC <- 1:40

sum.dt <- rbind(tb,sum.dt)
models <- c(list("PFS"=pfs),models)
rm(fit,RR,CI,P,tb,pfs,beta)

##Time to treatment failure on line 1
dt$censtf1 <- as.integer(dt$censtf1)
tf1 <- analyse_multivariate(data = dt,
                            time_status = vars(ttctf1,censtf1),
                            covariates = colnames(dplyr::select(dt,ends_with("_sd"))))

fit <- tf1$coxph
beta <- coef(fit)
RR <- exp(coef(fit))
CI <- exp(confint(fit))
colnames(CI) <- c("Lower", "Upper")
P <- coef(summary(fit))[,5]
tb <- as.data.table(cbind(beta,RR, CI, P))
tb$Y <- "TFTF"
tb$PC <- 1:40

sum.dt <- rbind(tb,sum.dt)
models <- c(list("TFTF"=tf1),models)
rm(fit,RR,CI,P,tb,tf1,beta)


# Disparities ----
##Race
mdt <- dt[D_PT_race%in%c(1,2)] %>% dplyr::select("D_PT_race",ends_with("_sd")) %>% na.omit()
fit <- glm(data = mdt, formula = D_PT_race ~ ., family = "binomial")
beta <- coef(fit)
RR <- exp(coef(fit))
CI <- exp(confint(fit))
colnames(CI) <- c("Lower", "Upper")
P <- coef(summary(fit))[,4]
tb <- as.data.table(cbind(beta,RR, CI, P)[-1,])
tb$Y <- "race"
tb$PC <- 1:40

sum.dt <- rbind(tb,sum.dt)
models <- c(list("race"=fit),models)
rm(mdt,fit,RR,CI,P,tb,beta)

#Ethnic
mdt <- dt[D_PT_ethnic%in%c(1,2)] %>% dplyr::select("D_PT_ethnic",ends_with("_sd")) %>% na.omit()
fit <- glm(data = mdt, formula = D_PT_ethnic ~ ., family = "binomial")
beta <- coef(fit)
RR <- exp(coef(fit))
CI <- exp(confint(fit))
colnames(CI) <- c("Lower", "Upper")
P <- coef(summary(fit))[,4]
tb <- as.data.table(cbind(beta,RR, CI, P)[-1,])
tb$Y <- "ethnic"
tb$PC <- 1:40

sum.dt <- rbind(tb,sum.dt)
models <- c(list("ethnic"=fit),models)
rm(mdt,fit,RR,CI,P,tb,beta)

#Gender
mdt <- dt[D_PT_gender%in%c(1,2)] %>% dplyr::select("D_PT_gender",ends_with("_sd")) %>% na.omit()
fit <- glm(data = mdt, formula = D_PT_gender ~ ., family = "binomial")
beta <- coef(fit)
RR <- exp(coef(fit))
CI <- exp(confint(fit))
colnames(CI) <- c("Lower", "Upper")
P <- coef(summary(fit))[,4]
tb <- as.data.table(cbind(beta,RR, CI, P)[-1,])
tb$Y <- "gender"
tb$PC <- 1:40

sum.dt <- rbind(tb,sum.dt)
models <- c(list("gender"=fit),models)
rm(mdt,fit,RR,CI,P,tb,beta)

#Age
mdt <- dt %>% dplyr::select("D_PT_age",ends_with("_sd")) %>% na.omit()
fit <- lm(data = mdt, formula = D_PT_age ~ .)
beta <- coef(fit)
RR <- exp(coef(fit))
CI <- exp(confint(fit))
colnames(CI) <- c("Lower", "Upper")
P <- coef(summary(fit))[,4]
tb <- as.data.table(cbind(beta,RR, CI, P)[-1,])
tb$Y <- "age"
tb$PC <- 1:40

sum.dt <- rbind(tb,sum.dt)
models <- c(list("age"=fit),models)
rm(mdt,fit,RR,CI,P,tb,beta)

# Tidy summary table ----
df <- as.data.frame(sum.dt[!sum.dt$Y%in%c("t(6;14)","t(14;20)"),] %>% na.omit()) #remove problematic regressions
df$sig = as.factor(if_else(df$P<0.05,1,0))
df$Direction = as.factor(if_else(df$sig==1 & df$RR<1,"-",if_else(df$sig==1 & df$RR>1,"+","ns")))
df$Effect = abs(df$beta)

ordr <- c("UAMC-70","d(17p)","t(14;16)","a(1q)","t(4;14)","t(11;14)","ISS","OS","PFS","TFTF","race","ethnic","gender","age")
df$Y <- factor(df$Y,levels = ordr)
df[df$sig==0,] <- NA 
#summary(df)
pldt <- droplevels(na.omit(df))
summary(pldt)

# Plot ----
theme_set(theme_minimal() + theme(legend.position = "right"))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#size based on risk ratio
ggplot(pldt, aes(y = PC, x = Y)) +
  geom_point(aes(colour = Direction,size=Effect),na.rm = T) +
  scale_color_manual(values = c("#E69F00","#0072B2")) +
  xlab("") + ylab("dimension") + scale_y_continuous(limits = c(1,40), breaks=seq(1,40,1)) +
  theme(panel.grid.minor = element_blank()) + theme(legend.position = "top")

#color based on direction of association
ggplot(pldt, aes(y = PC, x = Y)) +
  geom_point(aes(colour = dir),size=5,na.rm = T) +
  scale_color_manual(values = c("#E69F00","#0072B2")) +
  xlab("") + ylab("dimension") + scale_y_continuous(limits = c(1,40), breaks=seq(1,40,1)) +
  theme(panel.grid.minor = element_blank()) #+ theme(legend.position = "none")

#color based on risk ratio
ggplot(pldt, aes(y = PC, x = Y)) +
  geom_point(aes(colour = RR),size=5,na.rm = T) +
  scale_colour_gradient2(limits=c(min(pldt$RR),2),low="#E69F00",mid="white",high="#0072B2",midpoint=1,na.value="#0072B2") +
  xlab("") + ylab("dimension") + scale_y_continuous(limits = c(1,40), breaks=seq(1,40,1)) +
  theme(panel.grid.minor = element_blank())# + theme(legend.position = "none")

ggplot(pldt, aes(y = PC, x = Y)) +
  geom_point(shape = 21, colour = "gray", aes(fill = RR), size = 5, stroke = 1)+
  scale_fill_gradient2(limits=c(min(pldt$RR),2),low="#E69F00",mid="white",high="#0072B2",midpoint=1,na.value="#0072B2") +
  xlab("") + ylab("dimension") + scale_y_continuous(limits = c(1,40), breaks=seq(1,40,1)) +
  theme(panel.grid.minor = element_blank())# + theme(legend.position = "none")

#tile, color based on risk ratio
ggplot(pldt, aes(y = PC, x = Y, fill = RR)) +
  geom_tile() +
  scale_fill_gradient2(limits=c(min(pldt$RR),2),low="#E69F00",mid="white",high="#0072B2",midpoint=1,na.value="#0072B2") +
  xlab("") + ylab("dimension") + scale_y_continuous(limits = c(1,40), breaks=seq(1,40,1)) +
  theme(panel.grid.minor = element_blank())# + theme(legend.position = "none")

#save models
save(models,file = "rdata/regression_models-20200429.RData")
