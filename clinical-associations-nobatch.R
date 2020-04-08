#!/usr/bin/env Rscript

# Libraries
library(dplyr)
library(data.table)
library(aplpack)
library(survival)
library(survminer)
library(MASS)

# Functions

# Data
setwd("/Users/rosal/OneDrive - University of Utah/2020/career/analyze/data/transcriptome-dimensions/")
load(file='rdata/process_20200402.RData')
load(file = "rdata/clin_20200405.rdata")

# Setup
## merge clinical data and baseline sample pcs
dt <- merge(clin,pc.scores,by="sample_id")
rm(all.pc.scores,all.qc.dt,clin,elbow.dt,pc.scores,pca)

## create per-standard deviation pcs
dt$PC1_sd <- dt$PC1/sd(dt$PC1)
dt$PC2_sd <- dt$PC2/sd(dt$PC2)
dt$PC3_sd <- dt$PC3/sd(dt$PC3)
dt$PC4_sd <- dt$PC4/sd(dt$PC4)
dt$PC5_sd <- dt$PC5/sd(dt$PC5)
dt$PC6_sd <- dt$PC6/sd(dt$PC6)
dt$PC7_sd <- dt$PC7/sd(dt$PC7)
dt$PC8_sd <- dt$PC8/sd(dt$PC8)
dt$PC9_sd <- dt$PC9/sd(dt$PC9)
dt$PC10_sd <- dt$PC10/sd(dt$PC10)
dt$PC11_sd <- dt$PC11/sd(dt$PC11)
dt$PC12_sd <- dt$PC12/sd(dt$PC12)
dt$PC13_sd <- dt$PC13/sd(dt$PC13)
dt$PC14_sd <- dt$PC14/sd(dt$PC14)
dt$PC15_sd <- dt$PC15/sd(dt$PC15)
dt$PC16_sd <- dt$PC16/sd(dt$PC16)
dt$PC17_sd <- dt$PC17/sd(dt$PC17)
dt$PC18_sd <- dt$PC18/sd(dt$PC18)
dt$PC19_sd <- dt$PC19/sd(dt$PC19)
dt$PC20_sd <- dt$PC20/sd(dt$PC20)
dt$PC21_sd <- dt$PC21/sd(dt$PC21)
dt$PC22_sd <- dt$PC22/sd(dt$PC22)
dt$PC23_sd <- dt$PC23/sd(dt$PC23)
dt$PC24_sd <- dt$PC24/sd(dt$PC24)
dt$PC25_sd <- dt$PC25/sd(dt$PC25)
dt$PC26_sd <- dt$PC26/sd(dt$PC26)
dt$PC27_sd <- dt$PC27/sd(dt$PC27)
dt$PC28_sd <- dt$PC28/sd(dt$PC28)
dt$PC29_sd <- dt$PC29/sd(dt$PC29)
dt$PC30_sd <- dt$PC30/sd(dt$PC30)
dt$PC31_sd <- dt$PC31/sd(dt$PC31)

# Progression free survival
pfs.surv<-Surv(dt$ttcpfs, dt$censpfs) # object

# step-wise cox
stepvars<-c(colnames(dplyr::select(dt,ends_with('_sd'))),'ttcpfs','censpfs')
dt_reduced<-na.omit(dt[,stepvars])
pfs.step.cox<-stepAIC(coxph(pfs.surv~.,data=dt_reduced,id=dt$sample_id),
  direction='both',trace=0,k=log(sum(dt_reduced$censpfs)))
summary(pfs.step.cox)


# Translocations
Y <- colnames(dplyr::select(dt,ends_with('_sd')))
X <- colnames(dplyr::select(dt,contains('ABNORMALITY')))
DATA <- dt

mod <- expand.grid(y=Y, x=X, stringsAsFactors = T) %>%
  mutate(formula = paste(y,"~",x)) %>%
  group_by(formula) %>%
  mutate(F = summary(lm(formula, data=DATA))$fstatistic[1]) %>%
  mutate(DF1 = summary(lm(formula, data=DATA))$fstatistic[2]) %>%
  mutate(DF2 = summary(lm(formula, data=DATA))$fstatistic[3]) %>%
  mutate(R2 = summary(lm(formula, data=DATA))$r.squared) %>%
  mutate(BETA_SD = summary(lm(formula,data=DATA))$coefficients[2,1]) %>%
  ungroup()
mod$P <- pf(mod$F,mod$DF1,mod$DF2,lower.tail=F)
setDT(mod)

tsl.mod <- mod[,bon_sig:=as.factor('No')]
p <- 0.05/length(Y)
tsl.mod[P<p]$bon_sig <- as.factor('Yes')
setorder(tsl.mod,P)
rm(DATA,X,Y,mod,p)

# Patient demographics
Y <- colnames(dplyr::select(dt,ends_with('_sd')))
X <- colnames(dplyr::select(dt,starts_with('D_PT_')))
DATA <- dt

mod <- expand.grid(y=Y, x=X, stringsAsFactors = T) %>%
  mutate(formula = paste(y,"~",x)) %>%
  group_by(formula) %>%
  mutate(F = summary(lm(formula, data=DATA))$fstatistic[1]) %>%
  mutate(DF1 = summary(lm(formula, data=DATA))$fstatistic[2]) %>%
  mutate(DF2 = summary(lm(formula, data=DATA))$fstatistic[3]) %>%
  mutate(R2 = summary(lm(formula, data=DATA))$r.squared) %>%
  mutate(BETA_SD = summary(lm(formula,data=DATA))$coefficients[2,1]) %>%
  ungroup()
mod$P <- pf(mod$F,mod$DF1,mod$DF2,lower.tail=F)
setDT(mod)

dpt.mod <- mod[,bon_sig:=as.factor('No')]
p <- 0.05/length(Y)
dpt.mod[P<p]$bon_sig <- as.factor('Yes')
setorder(dpt.mod,P)
rm(DATA,X,Y,mod,p)

# Lab values
Y <- colnames(dplyr::select(dt,ends_with('_sd')))
X <- colnames(dplyr::select(dt,starts_with('D_LAB_')))
DATA <- dt

mod <- expand.grid(y=Y, x=X, stringsAsFactors = T) %>%
  mutate(formula = paste(y,"~",x)) %>%
  group_by(formula) %>%
  mutate(F = summary(lm(formula, data=DATA))$fstatistic[1]) %>%
  mutate(DF1 = summary(lm(formula, data=DATA))$fstatistic[2]) %>%
  mutate(DF2 = summary(lm(formula, data=DATA))$fstatistic[3]) %>%
  mutate(R2 = summary(lm(formula, data=DATA))$r.squared) %>%
  mutate(BETA_SD = summary(lm(formula,data=DATA))$coefficients[2,1]) %>%
  ungroup()
mod$P <- pf(mod$F,mod$DF1,mod$DF2,lower.tail=F)
setDT(mod)

lab.mod <- mod[,bon_sig:=as.factor('No')]
p <- 0.05/length(Y)
lab.mod[P<p]$bon_sig <- as.factor('Yes')
setorder(lab.mod,P)
rm(DATA,X,Y,mod,p)
