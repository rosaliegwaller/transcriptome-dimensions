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
  tmp<-data.table(t(samples))[-1]
  seq.ids<-tmp[V1%like%"_1_BM"] #select baseline samples only
  seq.ids$public_id<-gsub("_1_BM","",samples.dt$V1)
  colnames(seq.ids)[1]<-"seq_id"
  return(seq.ids)
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
write.csv(pat$key,file="rdata/MMRF_CoMMpass_IA14_PER_PATIENT_KEY.csv")
write.csv(pat$dt,file="rdata/MMRF_CoMMpass_IA14_PER_PATIENT_DT.csv")
write.csv(sur$key,file="rdata/MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL_KEY.csv")
write.csv(sur$dt,file="rdata/MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL_DT.csv")

## TO DO: additional clinical files (not per-patient)
