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
save(clin.dt,key,file = "rdata/setup-clinical-data-20200402.rdata")


## old
## get list of samples with seq data avaliable
#seq.ids<-seq_samples(
#  "MMRF_CoMMpass_IA14a_E74GTF_Salmon_V7.2_Filtered_Transcript_Counts.txt.gz")

## read in clinical data and count number missing in seq data
#pat<-clinical_data("MMRF_CoMMpass_IA14_PER_PATIENT.csv",seq.ids)
#sur<-clinical_data("MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL.csv",seq.ids)

## save counts of missing data in seq samples
#write.csv(pat$key,file="rdata/MMRF_CoMMpass_IA14_PER_PATIENT_KEY.csv")
#write.csv(pat$dt,file="rdata/MMRF_CoMMpass_IA14_PER_PATIENT_DT.csv")
#write.csv(sur$key,file="rdata/MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL_KEY.csv")
#write.csv(sur$dt,file="rdata/MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL_DT.csv")

## read in all data keys and merge
#key.pat <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_PER_PATIENT.csv")
#key.vis <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_PER_PATIENT_VISIT.csv")
#key.adm <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_ADMISSIONS.csv")
#key.ae <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_AE.csv")
#key.eme <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_EMERGENCY_DEPT.csv")
#key.fam <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_FAMHX.csv")
#key.med <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_MEDHX.csv")
#key.sur <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL.csv")
#key.reg <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_TREATMENT_REGIMEN.csv")
#key.res <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_TRTRESP.csv")

#keys <- unique(data.table(rbind(key.pat,key.vis,key.adm,key.ae,key.eme,key.fam,key.med,key.sur,key.reg,key.res))[,-1])
