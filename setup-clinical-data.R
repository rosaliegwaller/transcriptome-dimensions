#!/usr/bin/env Rscript

# Description
# Choose clinical data to investigate correlations

# Notes
# Download data from https://research.themmrf.org/
## CoMMpass_IA14_FlatFiles.tar.gz
## CoMMpass_IA14_FlatFile_Dictionaries.tar.gz

library(data.table)
library(naniar)

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

# Process
## get list of samples with seq data avaliable
seq.ids<-seq_samples(
  "MMRF_CoMMpass_IA14a_E74GTF_Salmon_V7.2_Filtered_Transcript_Counts.txt.gz")

## read in clinical data and count number missing in seq data
pat<-clinical_data("MMRF_CoMMpass_IA14_PER_PATIENT.csv",seq.ids)
sur<-clinical_data("MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL.csv",seq.ids)

## save counts of missing data in seq samples
#write.csv(pat$key,file="rdata/MMRF_CoMMpass_IA14_PER_PATIENT_KEY.csv")
#write.csv(pat$dt,file="rdata/MMRF_CoMMpass_IA14_PER_PATIENT_DT.csv")
#write.csv(sur$key,file="rdata/MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL_KEY.csv")
#write.csv(sur$dt,file="rdata/MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL_DT.csv")

## read in all data keys and merge
key.pat <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_PER_PATIENT.csv")
key.vis <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_PER_PATIENT_VISIT.csv")
key.adm <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_ADMISSIONS.csv")
key.ae <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_AE.csv")
key.eme <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_EMERGENCY_DEPT.csv")
key.fam <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_FAMHX.csv")
key.med <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_MEDHX.csv")
key.sur <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL.csv")
key.reg <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_TREATMENT_REGIMEN.csv")
key.res <- read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_STAND_ALONE_TRTRESP.csv")

keys <- unique(data.table(rbind(key.pat,key.vis,key.adm,key.ae,key.eme,key.fam,key.med,key.sur,key.reg,key.res))[,-1])

# SEQ QC Data
seq_qc <- data.table(read.csv("MMRF_CoMMpass_IA14_Seq_QC_Summary.csv"))[MMRF_Release_Status=="RNA-Yes"]

## TO DO: subset clinical data to seq samples from QC_Summary and select relevant data, save as rdata to load in pc analysis
seq_qc[,c("ï..Patients..KBase_Patient_ID","Visits..Study.Visit.ID",
  "Visits..Reason_For_Collection","QC.Link.SampleName","Batch","Creation.Date",
  "Process.Date","MMRF_Release_Status",
  "QC_Percent_IgH","QC_Percent_IgK","QC_Percent_IgL",
  "QC_Percent_Immunoglobulin","QC_Percent_Mapped","QC_Percent_Mitochondrial",
  "QC_Percent_Q20_Bases","QC_Percent_Uniq_Mapping","QC_Read_Counts")]
