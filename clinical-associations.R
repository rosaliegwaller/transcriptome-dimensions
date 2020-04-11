#!/usr/bin/env Rscript

## Libraries ------------------------------
library(dplyr)
library(data.table)
library(survival)
library(survminer)
library(jtools)
library(MASS)

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
  fom <- as.formula(paste(
    paste(Y," ~ ",sep = ""),
    paste(colnames(dplyr::select(dt,starts_with("PC"))),collapse = "+")
  ))
  mod <- glm(
    formula = fom,
    data = dt,
    family = "binomial",
    na.action = na.omit
  )
  stp <- stepAIC(mod, direction = "both", trace = FALSE)
  sum <- summ(stp,confint = TRUE, digits = 5)
  out <- list("formula"=fom,"glm"=mod,"stp"=stp,"sum"=sum)
  return(out)
}


## Data ------------------------------
setwd("/Users/rosal/OneDrive - University of Utah/2020/career/analyze/data/transcriptome-dimensions/")
load(file = "rdata/combat_pcs_20200411.rdata") #pcs and clinical data in dt

## SURVIVAL ------------------------------
# overall survival
dt$censos <- as.integer(dt$censos)
os <- analyse_multivariate(data = dt,
                           time_status = vars(ttcos,censos),
                           covariates = colnames(dplyr::select(dt,starts_with("PC")))[1:10]
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
