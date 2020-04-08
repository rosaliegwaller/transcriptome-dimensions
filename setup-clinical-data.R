#!/usr/bin/env Rscript

# Description
# Choose clinical data to investigate correlations

# Notes
# Download data from https://research.themmrf.org/
## CoMMpass_IA14_FlatFiles.tar.gz
## CoMMpass_IA14_FlatFile_Dictionaries.tar.gz

library(data.table)
library(naniar)
library(ggpubr)
library(dplyr)
library(tidyr)

# Functions
seq_samples <- function(counts.file){
  samples<-read.csv(file=counts.file,header=F,nrows=1,sep="\t")
  ids<-data.table(t(samples))[-1]
  #seq.ids<-tmp[V1%like%"_1_BM"] #select baseline samples only
  ids$public_id<-gsub("_([1-9])_(BM|PB)","",ids$V1)
  colnames(ids)[1]<-"seq_id"
  return(ids)
}

clinical_data <- function(clin.file,seq.ids){
  clin.dt<-data.table(read.csv(paste("CoMMpass_IA14_FlatFiles/",clin.file,sep="")))
  colnames(clin.dt)[1]<-"public_id"
  dt.seq<-merge(seq.ids,clin.dt,by='public_id')
  missing<-rbind(dt.seq[,lapply(.SD,function(x)sum(is.na(x)))],
    dt.seq[,lapply(.SD,function(x)sum(x==""))])
  m.dt<-data.table(t(setDT(missing)[,lapply(.SD,sum,na.rm=T)]))
  m.dt$name<-colnames(missing)
  colnames(m.dt)[1]<-"n_na_missing"
  key<-data.table(read.csv(paste("CoMMpass_IA14_FlatFile_Dictionaries/",clin.file,sep=""))[2:3])
  key.na<-merge.data.table(key,m.dt,by='name')
  key.na$n_data<-nrow(dt.seq)-key.na$n_na_missing
  key.na$prop_data<-key.na$n_data/nrow(dt.seq)
  setorder(key.na,n_data)
  out<-list("key"=key.na,"dt"=dt.seq)
  return(out)
}

count_missing <- function(key.dt,clin.dt){
  missing<-rbind(clin.dt[,lapply(.SD,function(x)sum(is.na(x)))],
    clin.dt[,lapply(.SD,function(x)sum(x==""))])
  m.dt<-data.table(t(setDT(missing)[,lapply(.SD,sum,na.rm=T)]))
  m.dt$name<-colnames(missing)
  colnames(m.dt)[1]<-"n_na_missing"
  key.na<-merge.data.table(key.dt,m.dt,by='name')
  key.na$n_data<-nrow(clin.dt)-key.na$n_na_missing
  key.na$prop_data<-key.na$n_data/nrow(clin.dt)
  setorder(key.na,-n_data)
  return(key.na[,-c(4:5)])
}

# Process

## sequencing quality control data
seq_qc <- data.table(read.csv("MMRF_CoMMpass_IA14_Seq_QC_Summary.csv"))[MMRF_Release_Status=="RNA-Yes",c(1:5)]
colnames(seq_qc) <- c("public_id","visit_id","collection_reason","sample_id","batch")
seq_qc$sample_id <- gsub("_CD138pos_T([1-9]{1,2})_TSMRU_([L|K0-9]{6})","",seq_qc$sample_id)

## merge with per-patient data
pat<-data.table(read.csv("CoMMpass_IA14_FlatFiles/MMRF_CoMMpass_IA14_PER_PATIENT.csv"))
colnames(pat)[1]<-"public_id"
seq.pat<-inner_join(seq_qc,pat,by='public_id')

## merge with per-visit data
vis<-data.table(read.csv("CoMMpass_IA14_FlatFiles/MMRF_CoMMpass_IA14_PER_PATIENT_VISIT.csv"))
colnames(vis)[1:2]<-c("public_id","visit_id")
seq.pat.vis<-inner_join(seq.pat,vis,by=c('public_id','visit_id'))

## merge with survival data
sur<-data.table(read.csv("CoMMpass_IA14_FlatFiles/MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL.csv"))
seq.pat.vis.sur<-inner_join(seq.pat.vis,sur,by='public_id')

clin.dt <- setDT(seq.pat.vis.sur)

## data descriptions
key.pat <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_PER_PATIENT.csv")
key.vis <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_PER_PATIENT_VISIT.csv")
key.sur <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL.csv")

key.dt <- unique(data.table(rbind(key.pat,key.vis,key.sur))[,-1])

## count missing data and annotate to keys
key<-count_missing(key.dt,clin.dt)

## and save
#save(clin.dt,key,file = "rdata/setup-clinical-data-20200402.rdata")

## subset clinical data
vars <- c(
  "sample_id",
  "collection_reason",
  "batch",
  "D_PT_age",
  "D_PT_gender",
  "D_PT_race",
  "D_PT_ethnic",
  "ttcos", #Time to OS event (censored)
  "censos", #Censor flag: overall survival
  "ttcpfs", #Time to PFS event (censored)
  "censpfs", #Censor flag: progression-free survival
  "ttctf1", #Line 1 time to treatment failure (censored)
  "censtf1", #Line 1 censor flag: time to treatment failure
  #"lstalive", #Last known alive
  "D_PT_iss", #ISS Disease Stage
  "D_LAB_serum_m_protein", #Serum M Protein (g/dL)
  "D_LAB_serum_kappa", #Serum Kappa (mg/dL)
  "D_LAB_serum_lambda", #Serum Lambda (mg/dL)
  "D_LAB_serum_beta2_microglobulin", #Beta 2 Microglobulin (mcg/mL)
  "D_TRI_CF_ABNORMALITYPR6",	#t(11;14) abnormality present standard risk
  "D_TRI_CF_ABNORMALITYPR4",	#t(6;14) abnormality present standard risk
  "D_TRI_CF_ABNORMALITYPR3",	#t(4;14) abnormality present intermediate risk
  "D_TRI_CF_ABNORMALITYPR13",	#1q amplification abnormality present intermediate risk
  "D_TRI_CF_ABNORMALITYPR8",	#t(14;16) abnormality present high risk
  "D_TRI_CF_ABNORMALITYPR9",	#t(14;20) abnormality present high risk
  "D_TRI_CF_ABNORMALITYPR11"	#del 17p abnormality present high risk
)

clin <- select(clin.dt,all_of(vars))

## clean up selected data
clin$collection_reason <- as.factor(clin$collection_reason)
clin$batch <- as.factor(clin$batch)
clin$D_PT_gender <- as.factor(clin$D_PT_gender)
clin$D_PT_race <- as.factor(clin$D_PT_race)
clin$D_PT_ethnic <- as.factor(clin$D_PT_ethnic)
clin$censos <- as.factor(clin$censos)
clin$censpfs <- as.factor(clin$censpfs)
clin$censtf1 <- as.factor(clin$censtf1)
clin$D_PT_iss <- as.factor(clin$D_PT_iss)
clin[clin==""] <- NA
clin[clin=="Not Done"] <- NA
clin$D_TRI_CF_ABNORMALITYPR6 <- as.factor(clin$D_TRI_CF_ABNORMALITYPR6)
clin$D_TRI_CF_ABNORMALITYPR4 <- as.factor(clin$D_TRI_CF_ABNORMALITYPR4)
clin$D_TRI_CF_ABNORMALITYPR3 <- as.factor(clin$D_TRI_CF_ABNORMALITYPR3)
clin$D_TRI_CF_ABNORMALITYPR13 <- as.factor(clin$D_TRI_CF_ABNORMALITYPR13)
clin$D_TRI_CF_ABNORMALITYPR8 <- as.factor(clin$D_TRI_CF_ABNORMALITYPR8)
clin$D_TRI_CF_ABNORMALITYPR9 <- as.factor(clin$D_TRI_CF_ABNORMALITYPR9)
clin$D_TRI_CF_ABNORMALITYPR11 <- as.factor(clin$D_TRI_CF_ABNORMALITYPR11)

save(clin,file = "rdata/clin_20200405.rdata")

## seperate factors to individual columns
load(file = "rdata/clin_20200405.rdata")
clin[,baseline:=0]
clin[collection_reason=="Baseline"]$baseline <- 1
clin[,confirm_progression:=0]
clin[collection_reason=="Confirm Progression"]$confirm_progression <- 1
clin[,post_transplant:=0]
clin[collection_reason=="Post Transplant"]$post_transplant <- 1
clin[,confirm_response:=0]
clin[collection_reason=="Confirm Response"]$confirm_response <- 1

clin[,female:=0]
clin[D_PT_gender==2]$female <- 1

clin[,race.1:=0]
clin[D_PT_race==1]$race.1 <- 1
clin[is.na(D_PT_race)]$race.1 <- -9
clin[,race.2:=0]
clin[D_PT_race==2]$race.2 <- 1
clin[is.na(D_PT_race)]$race.2 <- -9
clin[,race.4:=0]
clin[D_PT_race==4]$race.4 <- 1
clin[is.na(D_PT_race)]$race.4 <- -9
clin[,race.6:=0]
clin[D_PT_race==6]$race.6 <- 1
clin[is.na(D_PT_race)]$race.6 <- -9

clin[,ethnic.1:=0]
clin[D_PT_ethnic==1]$ethnic.1 <- 1
clin[is.na(D_PT_ethnic)]$ethnic.1 <- -9
clin[,ethnic.2:=0]
clin[D_PT_ethnic==2]$ethnic.2 <- 1
clin[is.na(D_PT_ethnic)]$ethnic.2 <- -9
clin[,ethnic.3:=0]
clin[D_PT_ethnic==3]$ethnic.3 <- 1
clin[is.na(D_PT_ethnic)]$ethnic.3 <- -9

clin[,iss.1:=0]
clin[D_PT_iss==1]$iss.1 <- 1
clin[is.na(D_PT_iss)]$iss.1 <- -9
clin[,iss.2:=0]
clin[D_PT_iss==2]$iss.2 <- 1
clin[is.na(D_PT_iss)]$iss.2 <- -9
clin[,iss.3:=0]
clin[D_PT_iss==3]$iss.3 <- 1
clin[is.na(D_PT_iss)]$iss.3 <- -9

clin[,t_11_14:=0]
clin[D_TRI_CF_ABNORMALITYPR6=='Yes']$t_11_14 <- 1
clin[is.na(D_TRI_CF_ABNORMALITYPR6)]$t_11_14 <- -9

clin[,t_6_14:=0]
clin[D_TRI_CF_ABNORMALITYPR4=='Yes']$t_6_14 <- 1
clin[is.na(D_TRI_CF_ABNORMALITYPR4)]$t_6_14 <- -9

clin[,t_4_14:=0]
clin[D_TRI_CF_ABNORMALITYPR3=='Yes']$t_4_14 <- 1
clin[is.na(D_TRI_CF_ABNORMALITYPR3)]$t_4_14 <- -9

clin[,amp_1q:=0]
clin[D_TRI_CF_ABNORMALITYPR13=='Yes']$amp_1q <- 1
clin[is.na(D_TRI_CF_ABNORMALITYPR13)]$amp_1q <- -9

clin[,t_14_16:=0]
clin[D_TRI_CF_ABNORMALITYPR8=='Yes']$t_14_16 <- 1
clin[is.na(D_TRI_CF_ABNORMALITYPR8)]$t_14_16 <- -9

clin[,t_14_20:=0]
clin[D_TRI_CF_ABNORMALITYPR9=='Yes']$t_14_20 <- 1
clin[is.na(D_TRI_CF_ABNORMALITYPR9)]$t_14_20 <- -9

clin[,del_17p:=0]
clin[D_TRI_CF_ABNORMALITYPR11=='Yes']$del_17p <- 1
clin[is.na(D_TRI_CF_ABNORMALITYPR11)]$del_17p <- -9

clin$censos <- as.integer(clin$censos)
clin$censpfs <- as.integer(clin$censpfs)
clin$censtf1 <- as.integer(clin$censtf1)

clin_cols <- clin[,c(1,3:4,8:13,15:18,26:47)]
clin_cols[is.na(clin_cols)] <- -9

save(clin_cols,file="rdata/clinical_data_factors_as_columns.rdata")
