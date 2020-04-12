#!/usr/bin/env Rscript

## Libraries ------------------------------
library(dplyr)
library(data.table)
library(survival)
library(survminer)
library(jtools)
library(MASS)
library(ggcorrplot)
library(reghelper)
library(survivalAnalysis)

## Functions ------------------------------
run_lm <- function(Y,dt){
  fom <- as.formula(paste(
    paste(Y," ~ ",sep = ""),
    paste(colnames(dplyr::select(dt,starts_with("PC"))),collapse = "+")
  ))
  mod <- lm(
    formula = fom,
    data = dt,
    na.action = na.omit
  )
  stp <- stepAIC(mod, direction = "both", trace = FALSE)
  sum <- summ(stp,confint = TRUE, digits = 5)
  out <- list("formula"=fom,"lm"=mod,"stp"=stp,"sum"=sum)
  return(out)
}

run_glm <- function(Y,dt){
  mod <- glm(
    formula = as.formula(paste(
      paste(Y," ~ ",sep = ""),
      paste(colnames(dplyr::select(dt,starts_with("PC"))),collapse = "+"))),
    data = dt,
    family = "binomial",
    na.action = na.omit
  )
  out <- list("glm"=mod)
  return(out)
}


## Data ------------------------------
setwd("/Users/rosal/OneDrive - University of Utah/2020/career/analyze/data/transcriptome-dimensions/")
load(file = "rdata/combat_pcs_20200411.rdata") #pcs and clinical data in dt

## TO DO: STANDARDIZE --------------------


## SURVIVAL ------------------------------
# overall survival
dt$censos <- as.integer(dt$censos)
os <- analyse_multivariate(data = dt,
                           time_status = vars(ttcos,censos),
                           covariates = colnames(dplyr::select(dt,starts_with("PC")))
                           )
forest_plot(os)
summary(os$coxph)

# progression free survival
dt$censpfs <- as.integer(dt$censpfs)
pfs <- analyse_multivariate(data = dt,
                           time_status = vars(ttcpfs,censpfs),
                           covariates = colnames(dplyr::select(dt,starts_with("PC")))[1:10]
                           )
forest_plot(pfs,orderer = ~order(HR))
summary(pfs$coxph)

# time to treatment failure on line 1
dt$censtf1 <- as.integer(dt$censtf1)
tf1 <- analyse_multivariate(data = dt,
                            time_status = vars(ttctf1,censtf1),
                            covariates = colnames(dplyr::select(dt,starts_with("PC")))[1:10]
                            )
forest_plot(tf1,orderer = ~order(HR))
summary(tf1$coxph)

# D_PT_iss
iss <- run_glm('D_PT_iss',dt)

## ABNORMALITIES ------------------------------
# IDEA: could create smart high, intermediate, and standard variables here
# "D_TRI_CF_ABNORMALITYPR6",	#t(11;14) abnormality present standard risk
t11_14 <- run_glm('D_TRI_CF_ABNORMALITYPR6',dt)

# "D_TRI_CF_ABNORMALITYPR4",	#t(6;14) abnormality present standard risk
t6_14 <- run_glm('D_TRI_CF_ABNORMALITYPR4',dt)

# "D_TRI_CF_ABNORMALITYPR3",	#t(4;14) abnormality present intermediate risk
t4_14 <- run_glm('D_TRI_CF_ABNORMALITYPR3',dt)

# "D_TRI_CF_ABNORMALITYPR13",	#1q amplification abnormality present intermediate risk
amp1q <- run_glm('D_TRI_CF_ABNORMALITYPR13',dt)

# "D_TRI_CF_ABNORMALITYPR8",	#t(14;16) abnormality present high risk
t14_16 <- run_glm('D_TRI_CF_ABNORMALITYPR8',dt)

# "D_TRI_CF_ABNORMALITYPR9",	#t(14;20) abnormality present high risk
t14_20 <- run_glm('D_TRI_CF_ABNORMALITYPR9',dt)

# "D_TRI_CF_ABNORMALITYPR11"	#del 17p abnormality present high risk
del17p <- run_glm('D_TRI_CF_ABNORMALITYPR11',dt)


## LAB ------------------------------
#NOTE: not sure the usefulness of these variables...
# D_LAB_serum_m_protein
mpr <- run_lm('D_LAB_serum_m_protein',dt)

# D_LAB_serum_kappa
kap <- run_lm('D_LAB_serum_kappa',dt)

# D_LAB_serum_lambda
lam <- run_lm('D_LAB_serum_lambda',dt)

# D_LAB_serum_beta2_microglobulin
b2mg <- run_lm('D_LAB_serum_beta2_microglobulin',dt)


## DEMOGRAPHICS ------------------------------
#D_PT_age
age <- run_lm('D_PT_age',dt)

#D_PT_gender
gender <- run_glm('D_PT_gender',dt)

#D_PT_race
race <- run_glm('D_PT_race',dt)

#D_PT_ethnic
ethnic <- run_glm('D_PT_ethnic',dt)


# plot another way
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
library(gridExtra)

#data(efc)
#theme_set(theme_sjplot())
fun_plot_model <- function(mod){
  p <- plot_model(
    model = mod,
    #type = "std2",
    show.values = TRUE,
    show.p = TRUE,
    value.offset = 0.5,
    value.size = 2,
    line.size = 0.5,
    dot.size = 1.5,
    digits = 3
  )
  return(p)
}

p <- fun_plot_model(gender$glm)
#p1 <- plot_model(age$lm,show.values=T,value.offset = 0.5,show.p=T,value.size = 2, line.size = 0.5, dot.size = 1.5, digits = 3) + ylim(-0.3,0.25)
p2 <- plot_model(gender$glm,show.values=T,value.offset = 0.5,show.p=T,value.size = 2, line.size = 0.5, dot.size = 1.5, digits = 3) + ylim(0.95,1.1)
p3 <- plot_model(race$glm,show.values=T,value.offset = 0.5,show.p=T,value.size = 2, line.size = 0.5, dot.size = 1.5, digits = 3) + ylim(0.95,1.1)
p4 <- plot_model(ethnic$glm,show.values=T,value.offset = 0.5,show.p=T,value.size = 2, line.size = 0.5, dot.size = 1.5, digits = 3) + ylim(0.95,1.1)
#p10 <- plot_model(os$coxph,show.values=T,value.offset = 0.5,show.p=T,value.size = 2, line.size = 0.5, dot.size = 1.5, digits = 3) + ylim(0.95,1.05)
grid.arrange(p2, p3, p4, nrow = 1)




## combine beta and p-values from each model
r <- data.table(cbind(exp(cbind(OR = coef(race$glm), confint(race$glm))),P = summary(race$glm)$coef[,"Pr(>|z|)"]))[-1,]
r$PC <- 1:45
r$var <- "race"

e <- data.table(cbind(exp(cbind(OR = coef(ethnic$glm), confint(ethnic$glm))),P = summary(ethnic$glm)$coef[,"Pr(>|z|)"]))[-1,]
e$PC <- 1:45
e$var <- "ethnic"

g <- data.table(cbind(exp(cbind(OR = coef(gender$glm), confint(gender$glm))),P = summary(gender$glm)$coef[,"Pr(>|z|)"]))[-1,]
g$PC <- 1:45
g$var <- "gender"

d <- rbind(r,e,g)
d[P>0.05][,1:3] <- NA

df <- as.data.frame(na.omit(d))
colnames(df) <- c("RiskRatio","LowerLimit","UpperLimit","P","Condition","Group")
## plot
ggplot(data=df,
           aes(x = Group,y = RiskRatio, ymin = LowerLimit, ymax = UpperLimit ))+
  geom_pointrange(aes(col=Group))+
  geom_hline(aes(fill=Group),yintercept =1, linetype=2)+
  xlab('Group')+ ylab("Risk Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=LowerLimit, ymax=UpperLimit,col=Group),width=0.5,cex=1)+ 
  facet_wrap(~Condition,strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()

## plot another way --------
dt.age <- data.table()
dt.age$PC <- 1:45 #colnames(dplyr::select(dt,starts_with("PC")))
dt.age$var <- "age"
dt.age$beta <- summary(age$lm)$coefficients[-1,"Estimate"]
dt.age$pval <- summary(age$lm)$coefficients[-1,"Pr(>|t|)"]

dt.gen <- data.table()
dt.gen$PC <- 1:45 #colnames(dplyr::select(dt,starts_with("PC")))
dt.gen$var <- "gender"
dt.gen$beta <- summary(gender$glm)$coefficients[-1,"Estimate"]
dt.gen$pval <- summary(gender$glm)$coefficients[-1,"Pr(>|z|)"]

dt.rac <- data.table()
dt.rac$PC <- 1:45 #colnames(dplyr::select(dt,starts_with("PC")))
dt.rac$var <- "race"
dt.rac$beta <- summary(race$glm)$coefficients[-1,"Estimate"]
dt.rac$pval <- summary(race$glm)$coefficients[-1,"Pr(>|z|)"]

dt.eth <- data.table()
dt.eth$PC <- 1:45 #colnames(dplyr::select(dt,starts_with("PC")))
dt.eth$var <- "ethnic"
dt.eth$beta <- summary(ethnic$glm)$coefficients[-1,"Estimate"]
dt.eth$pval <- summary(ethnic$glm)$coefficients[-1,"Pr(>|z|)"]

#TO DO: add surival and risk variables...

d <- rbind(dt.age,dt.gen,dt.rac,dt.eth)
d[pval>0.05] <- NA

df <- as.data.frame(na.omit(d))

ggplot(df, aes(y = PC, x = var)) +
  geom_point(aes(colour = beta, size = -pval),na.rm = T) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size(range = c(0, 6),breaks = c(0.01,0.001,0.0001,0.00001)) +
  xlab("") + ylab("") + scale_y_continuous(limits = c(1, 45), breaks=seq(1,45,1)) +
  theme_bw()

ggplot(df, aes(y = PC, x = var)) +
  geom_point(aes(colour = pval, size = -pval),na.rm = T) +
  scale_colour_gradient(high = "white", low = "blue") +
  scale_size(range = c(0, 6),breaks = c(0.01,0.001,0.0001,0.00001)) +
  xlab("") + ylab("") + scale_y_continuous(limits = c(1, 45), breaks=seq(1,45,1)) +
  theme_bw()
