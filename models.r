# regression and plot

# packages
library(dplyr)
library(MASS)
library(survivalAnalysis)
library(data.table)
library(ggplot2)
library(glmnet)
library(gridExtra)

# data
setwd("/Users/rosal/OneDrive - University of Utah/2020/analyze/data/transcriptome-dimensions/")
load(file='rdata/process_20200428.RData') #pcs and clinical data in dt
dt <- merge(samples,score_sd,by='sample_id')
rm(samples,score_sd,pc)

# Shaughnessy 70-gene risk score ----
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
hist(x$score,breaks = 200,main = "UQ-TPM Risk Score",xlab = "log2(mean up reg) - log2(mean down reg)")

# Not in data down: "PNPLA4","KIAA1754","MCLC","AD-020","PARG1","RFP2""FLJ20489"
# Not in data up: "FABP5","PDHA1","TRIP13","SELI","SLCI19A1","ASPM","STK6","FLJ13052","LAS1L","BIRC5","CKAP1","MGC57827","DKFZp779O175","PAPD1","EIF2C2","MGC4308","DSG2","C6orf173","MGC15606","KIF14","DKFZP586L0724","WEE1","ROBO1","MPHOSPH1"

# run linear regression
mdt <- merge(x,dt,by="sample_id") %>% dplyr::select("score",ends_with("_sd")) %>% na.omit()
score.lm <- lm(data = mdt, formula = score ~ .)
summary(score.lm)
plot(score.lm$residuals, pch = 16, col = "red")

# penalized - stepAIC
#score.aic <- stepAIC(score.lm,direction = "both",trace = F)

# penalized - glmnet
#score.net <- glmnet(x)

tb <- cbind(
  round(exp(coef(score.lm)),2),
  round(exp(confint(score.lm)),2),
  formatC(digits = 2, format = "e",x = coef(summary(score.lm))[,4]))
tb <- as.data.table(tb)#[-1,])
colnames(tb) <- c("RR","Lower","Upper","P")
tb$Y <- "Shaughnessy"
tb$PC <- c(0:45)
tb$RR <- as.numeric(tb$RR)
tb$Lower <- as.numeric(tb$Lower)
tb$Upper <- as.numeric(tb$Upper)
tb$P <- as.numeric(tb$P)

# plot
df.sha <- as.data.frame(tb[Y=="Shaughnessy",])
plot.sha <- ggplot(df.sha, aes(y = PC, x = Y, fill = RR)) +
  geom_tile() +
  geom_text(size = 3, aes(label = RR)) + # aes(label = formatC(x = P,digits = 1,format = "e"))) +
  scale_fill_gradient2(limits = c(0.7,1.3), low = "purple", mid = "white", high = "orange", midpoint = 1) +
  xlab("") + ylab("PC") + #scale_y_continuous(limits = c(0, 45), breaks=seq(0,45,1)) + 
  theme_bw() + theme(panel.grid.major = element_blank(),axis.line = element_line(colour = "black"))

grid.arrange(plot.sha,nrow=1)

# ISS score ----
##ordinal logistic regression model
mdt <- dt %>% dplyr::select("D_PT_iss",ends_with("_sd")) %>% na.omit()
model.iss = polr(data = mdt, formula = D_PT_iss ~ ., Hess = T)
summary(model.iss)

glm.iss <- glm(data = mdt, formula = D_PT_iss ~ ., family = "binomial")
summary(glm.iss)
exp(coef(glm.iss))
exp(confint(glm.iss))

# overall survival ----
dt$censos <- as.integer(dt$censos)
os <- analyse_multivariate(data = dt,
                           time_status = vars(ttcos,censos),
                           covariates = colnames(dplyr::select(dt,ends_with("_sd"))))
summary(os$coxph)
summary(os$coxph)$logtest

#plot_model(os$coxph) + ylim(0.95,1.05)
forest_plot(os,labels_displayed = "factor",label_headers = c("","PC"),title = "OS") #values_displayed = "p",

# progression free survival ----
dt$censpfs <- as.integer(dt$censpfs)
pfs <- analyse_multivariate(data = dt,
                            time_status = vars(ttcpfs,censpfs),
                            covariates = colnames(dplyr::select(dt,ends_with("_sd")))
)
summary(pfs$coxph)$logtest

forest_plot(pfs,labels_displayed = "factor",
            #values_displayed = "p",
            label_headers = c("","PC"),title = "PFS")

# time to treatment failure on line 1 ----
dt$censtf1 <- as.integer(dt$censtf1)
tf1 <- analyse_multivariate(data = dt,
                            time_status = vars(ttctf1,censtf1),
                            covariates = colnames(dplyr::select(dt,ends_with("_sd")))
)
summary(tf1$coxph)$logtest

forest_plot(tf1,labels_displayed = "factor",
            #values_displayed = "p",
            label_headers = c("","PC"),title = "Treatment failure")

fit <- tf1$coxph
HR <- round(exp(coef(fit)),2)
CI <- round(exp(confint(fit)),2)
colnames(CI) <- c("Lower", "Higher")
P <- round(coef(summary(fit))[,5],5)
tb2 <- as.data.table(cbind(HR, CI, P))
tb2$PC <- rownames(tb2)
tb2$Var <- "ttctf1"
tb <- rbind(tb,tb2)



# FISH ----
# "D_TRI_CF_ABNORMALITYPR6",	#t(11;14) abnormality present standard risk ----
mdt <- dt %>% dplyr::select("D_TRI_CF_ABNORMALITYPR6",ends_with("_sd")) %>% na.omit()
t11_14 <- glm(data = mdt, formula = D_TRI_CF_ABNORMALITYPR6 ~ ., family = "binomial")

summary(t11_14)

fit <- t11_14
RR <- round(exp(coef(fit)),2)
CI <- round(exp(confint(fit)),2)
colnames(CI) <- c("Lower", "Higher")
P <- round(coef(summary(fit))[,4],5)
tb <- as.data.table(cbind(RR, CI, P)[-1,])
tb$PC <- 1:45
tb$Var <- "t11_14"


tb[P>=0.05] <- NA
df <- as.data.frame(na.omit(tb))

ggplot(df, aes(y = PC, x = Var)) +
  geom_point(aes(colour = RR,size=5),na.rm = T) +
  scale_colour_gradient2(low = "purple", mid = "white", high = "orange", midpoint = 1) +
  xlab("") + ylab("PC") + scale_y_continuous(limits = c(1, 45), breaks=seq(1,45,1)) +
  theme_bw() + theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank())#,panel.grid.minor = element_blank())


# "D_TRI_CF_ABNORMALITYPR4",	#t(6;14) abnormality present standard risk
mdt <- dt %>% dplyr::select("D_TRI_CF_ABNORMALITYPR4",ends_with("_sd")) %>% na.omit()
t_6_14 <- glm(data = mdt, formula = D_TRI_CF_ABNORMALITYPR4 ~ ., family = "binomial")

summary(t_6_14)

fit <- t_6_14
RR <- round(exp(coef(fit)),2)
CI <- round(exp(confint(fit)),2)
colnames(CI) <- c("Lower", "Higher")
P <- round(coef(summary(fit))[,4],5)
tb <- as.data.table(cbind(RR, CI, P)[-1,])
tb$PC <- 1:45
tb$Var <- "t_6_14"

tb[P>=0.05] <- NA
df <- as.data.frame(na.omit(tb))

ggplot(df, aes(y = PC, x = Var)) +
  geom_point(aes(colour = RR,size=5),na.rm = T) +
  scale_colour_gradient2(low = "purple", mid = "white", high = "orange", midpoint = 1) +
  xlab("") + ylab("PC") + scale_y_continuous(limits = c(1, 45), breaks=seq(1,45,1)) +
  theme_bw() + theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank())#,panel.grid.minor = element_blank())

# "D_TRI_CF_ABNORMALITYPR3",	#t(4;14) abnormality present intermediate risk
#t4_14 <- run_glm('D_TRI_CF_ABNORMALITYPR3',dt)
mdt <- dt %>% dplyr::select("D_TRI_CF_ABNORMALITYPR3",ends_with("_sd")) %>% na.omit()
t4_14 <- glm(data = mdt, formula = D_TRI_CF_ABNORMALITYPR3 ~ ., family = "binomial")

summary(t4_14)

fit <- t4_14
RR <- round(exp(coef(fit)),2)
CI <- round(exp(confint(fit)),2)
colnames(CI) <- c("Lower", "Higher")
P <- round(coef(summary(fit))[,4],5)
tb <- as.data.table(cbind(RR, CI, P)[-1,])
tb$PC <- 1:45
tb$Var <- "t4_14"

tb[P>=0.05] <- NA
df <- as.data.frame(na.omit(tb))

ggplot(df, aes(y = PC, x = Var)) +
  geom_point(aes(colour = RR,size=5),na.rm = T) +
  scale_colour_gradient2(low = "purple", mid = "white", high = "orange", midpoint = 1) +
  xlab("") + ylab("PC") + scale_y_continuous(limits = c(1, 45), breaks=seq(1,45,1)) +
  theme_bw() + theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank())#,panel.grid.minor = element_blank())

# "D_TRI_CF_ABNORMALITYPR13",	#1q amplification abnormality present intermediate risk
#amp1q <- run_glm('D_TRI_CF_ABNORMALITYPR13',dt)
mdt <- dt %>% dplyr::select("D_TRI_CF_ABNORMALITYPR13",ends_with("_sd")) %>% na.omit()
amp1q <- glm(data = mdt, formula = D_TRI_CF_ABNORMALITYPR13 ~ ., family = "binomial")

summary(amp1q)

fit <- amp1q
RR <- round(exp(coef(fit)),2)
CI <- round(exp(confint(fit)),2)
colnames(CI) <- c("Lower", "Higher")
P <- round(coef(summary(fit))[,4],5)
tb <- as.data.table(cbind(RR, CI, P)[-1,])
tb$PC <- 1:45
tb$Var <- "amp1q"

tb[P>=0.05] <- NA
df <- as.data.frame(na.omit(tb))

ggplot(df, aes(y = PC, x = Var)) +
  geom_point(aes(colour = RR,size=5),na.rm = T) +
  scale_colour_gradient2(low = "purple", mid = "white", high = "orange", midpoint = 1) +
  xlab("") + ylab("PC") + scale_y_continuous(limits = c(1, 45), breaks=seq(1,45,1)) +
  theme_bw() + theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank())#,panel.grid.minor = element_blank())


# "D_TRI_CF_ABNORMALITYPR8",	#t(14;16) abnormality present high risk
#t14_16 <- run_glm('D_TRI_CF_ABNORMALITYPR8',dt)
mdt <- dt %>% dplyr::select("D_TRI_CF_ABNORMALITYPR8",ends_with("_sd")) %>% na.omit()
t14_16 <- glm(data = mdt, formula = D_TRI_CF_ABNORMALITYPR8 ~ ., family = "binomial")

summary(t14_16)

fit <- t14_16
RR <- round(exp(coef(fit)),2)
CI <- round(exp(confint(fit)),2)
colnames(CI) <- c("Lower", "Higher")
P <- round(coef(summary(fit))[,4],5)
tb <- as.data.table(cbind(RR, CI, P)[-1,])
tb$PC <- 1:45
tb$Var <- "t14_16"

tb[P>=0.05] <- NA
df <- as.data.frame(na.omit(tb))

ggplot(df, aes(y = PC, x = Var)) +
  geom_point(aes(colour = RR,size=5),na.rm = T) +
  scale_colour_gradient2(low = "purple", mid = "white", high = "orange", midpoint = 1) +
  xlab("") + ylab("PC") + scale_y_continuous(limits = c(1, 45), breaks=seq(1,45,1)) +
  theme_bw() + theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank())#,panel.grid.minor = element_blank())

# "D_TRI_CF_ABNORMALITYPR9",	#t(14;20) abnormality present high risk
#t14_20 <- run_glm('D_TRI_CF_ABNORMALITYPR9',dt)
mdt <- dt %>% dplyr::select("D_TRI_CF_ABNORMALITYPR9",ends_with("_sd")) %>% na.omit()
t14_20 <- glm(data = mdt, formula = D_TRI_CF_ABNORMALITYPR9 ~ ., family = "binomial")

summary(t14_20)

fit <- t14_20
RR <- round(exp(coef(fit)),2)
CI <- round(exp(confint(fit)),2)
colnames(CI) <- c("Lower", "Higher")
P <- round(coef(summary(fit))[,4],5)
tb <- as.data.table(cbind(RR, CI, P)[-1,])
tb$PC <- 1:45
tb$Var <- "t14_20"

tb[P>=0.05] <- NA
df <- as.data.frame(na.omit(tb))

ggplot(df, aes(y = PC, x = Var)) +
  geom_point(aes(colour = RR,size=5),na.rm = T) +
  scale_colour_gradient2(low = "purple", mid = "white", high = "orange", midpoint = 1) +
  xlab("") + ylab("PC") + scale_y_continuous(limits = c(1, 45), breaks=seq(1,45,1)) +
  theme_bw() + theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank())#,panel.grid.minor = element_blank())


# "D_TRI_CF_ABNORMALITYPR11"	#del 17p abnormality present high risk
#del17p <- run_glm('D_TRI_CF_ABNORMALITYPR11',dt)
mdt <- dt %>% dplyr::select("D_TRI_CF_ABNORMALITYPR11",ends_with("_sd")) %>% na.omit()
del17p <- glm(data = mdt, formula = D_TRI_CF_ABNORMALITYPR11 ~ ., family = "binomial")

summary(del17p)

fit <- del17p
RR <- round(exp(coef(fit)),2)
CI <- round(exp(confint(fit)),2)
colnames(CI) <- c("Lower", "Higher")
P <- round(coef(summary(fit))[,4],5)
tb <- as.data.table(cbind(RR, CI, P)[-1,])
tb$PC <- 1:45
tb$Var <- "del17p"

tb[P>=0.05] <- NA
df <- as.data.frame(na.omit(tb))

ggplot(df, aes(y = PC, x = Var)) +
  geom_point(aes(colour = RR,size=5),na.rm = T) +
  scale_colour_gradient2(low = "purple", mid = "white", high = "orange", midpoint = 1) +
  xlab("") + ylab("PC") + scale_y_continuous(limits = c(1, 45), breaks=seq(1,45,1)) +
  theme_bw() + theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank())#,panel.grid.minor = element_blank())


# plot
#tb[P>=0.05] <- NA
df <- as.data.frame(na.omit(tb))

ggplot(df, aes(y = PC, x = Var)) +
  geom_point(aes(colour = RR,size=5),na.rm = T) +
  scale_colour_gradient2(low = "purple", mid = "white", high = "orange", midpoint = 1) +
  xlab("FISH Abnormality") + ylab("PC") + scale_y_continuous(limits = c(1, 45), breaks=seq(1,45,1)) +
  theme_bw() + theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank())#,panel.grid.minor = element_blank())


# gender ----
dt <- mutate(dt,female = ifelse(D_PT_gender==2,1,0))
dt$female <- as.factor(dt$female)

mdt <- dt %>% dplyr::select("female",ends_with("_sd")) %>% na.omit()
mod <- glm(data = mdt, formula = female ~ ., family = "binomial")
summary(mod)
tb <- as.data.table(cbind(
  round(exp(coef(mod)),2),
  round(exp(confint(mod)),2),
  round(coef(summary(mod))[,4],5)))
colnames(tb) <- c("OR","Lower","Upper","P")
tb$var <- "female"
tb$PC <- c("Intercept",1:45)
tb

# plot

## heat map
tb.na <- tb
tb.na[P>=0.05]$OR <- 1
df <- as.data.frame(tb.na)
ggplot(df, aes(y = PC, x = var, fill = OR)) +
  geom_tile() +
  scale_fill_gradient2(low = "purple", mid = "white", high = "orange", midpoint = 1) +
  xlab("") + ylab("PC") + #scale_y_continuous(limits = c(1, 45), breaks=seq(1,45,1)) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


# forest plot
df <- as.data.frame(tb)
ggplot(df, aes(x=PC, y=OR, ymin=Lower, ymax=Upper))+
  geom_pointrange()+
  geom_hline(yintercept = 1, linetype=2)+
  #scale_x_continuous(limits = c(1, 45), breaks=seq(1,45,1)) +
  xlab('PC') + ylab('female') +
  coord_flip()
