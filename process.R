#!/usr/bin/env Rscript

# Libraries
library(rtracklayer)
library(dplyr)
library(data.table)

# Functions
get_data <- function(gtf.path,counts.path){
  gtf <- readGFF(filepath=gtf.path) #read in GTF
  gtf$length <- gtf$end - gtf$start + 1 #compute length of each feature
  tmp <- subset(gtf,type=='exon',)
  sum_exons<-aggregate(length ~ seqid + gene_name + gene_id + gene_biotype + transcript_id,
    data=tmp, FUN=sum)
  n_name <- sum_exons %>% count(gene_name) #count number transcripts/gene_name
  colnames(n_name)[2] <- "name_n_transcripts"
  gtf.dt <- data.table(inner_join(n_name,sum_exons,by='gene_name'))

  # read in rnaseq transcript counts
  counts <- read.csv(file = counts.path,header = T,sep = "\t") #read in counts
  setDT(counts)
  colnames(counts)[1]<-"transcript_id"

  # merge with gtf data
  counts.gtf <- inner_join(gtf.dt,counts,by="transcript_id")
  setDT(counts.gtf)
  return(counts.gtf)
}

quality_control <- function(counts.gtf,baseline.samples){
  # select protein coding exons on autosomal chromosomes and baseline samples
  tmp <- counts.gtf[seqid %in% c(1:22) & gene_biotype=='protein_coding'] %>% select("gene_name",contains("MMRF"))
  # aggregate trasncript counts to gene_name counts
  gene <- data.table(aggregate(. ~ gene_name,data=tmp,FUN=sum))
  # find proportion of genes with <100 counts in baseline samples only
  base <- gene %>% select("gene_name",all_of(baseline.samples))
  base$prop_base_lt100 <- rowSums(base[,-1] < 100) / length(baseline.samples)

  keep.genes <- base[prop_base_lt100 < 0.05]$gene_name

  # count all samples total reads and proportion of qc genes with > 100 reads
  total_reads <- gene[,-1][,lapply(.SD,sum)]
  total_reads_keep_genes <- gene[gene_name %in% keep.genes,-1][,lapply(.SD,sum)]
  prop_gene_lt100 <- gene[gene_name %in% keep.genes,-1][,lapply(.SD,function(x)sum(x<100)/length(keep.genes))]
  samples<-data.table(t((rbind(total_reads,total_reads_keep_genes,prop_gene_lt100))))
  colnames(samples)<-c("total_reads","total_reads_keep_genes","prop_gene_lt100")
  samples$sample_id <- colnames(total_reads)

  out <- list("keep.genes"=keep.genes,"samples"=samples)
  return(out)
}

process_transcripts <- function(qc.melt.dt){
  temp.dt <- qc.melt.dt[,total_raw_counts:=sum(count),by=c('sample_id','gene_name')]
  temp.dt[,kb_length:=length/1000]
  final.dt<-temp.dt[,list(cpk=sum((count+1/name_n_transcripts)/kb_length)),by=c('sample_id','gene_name','total_raw_counts')]

  final.dt[,size_factor:=median(cpk),by='sample_id']
  final.dt[,cpkmed:=cpk/size_factor]
  final.dt[,logcpkmed:=log2(cpkmed)]

  med.molten <- final.dt[,mean:=mean(logcpkmed),by='gene_name']
  med.molten[,sd:=sd(logcpkmed),by='gene_name']
  med.molten[,adjlogcpkmed:=logcpkmed]
  med.molten[(logcpkmed-mean)/sd>=5,adjlogcpkmed:=mean+5*sd]
  med.molten[(logcpkmed-mean)/sd<= -5,adjlogcpkmed:=mean-5*sd]

  out<-list("melt"=med.molten)
  return(out)
}

elbow_finder<-function(data) {
  elbow.dt<-data[order(-value)][,idx:=.I-1]
  elbow.dt[,selected:=as.factor('Not used')]

  slope<-(min(elbow.dt$value)-max(elbow.dt$value))/(nrow(elbow.dt)-1)
  perpslope<-(-1/slope)
  intercept<-max(elbow.dt$value) # Is this accurate?
  elbow.dt[,perpcept:=value - perpslope*idx]

  elbow.dt[,y:=(perpcept*slope - intercept*perpslope)/(slope-perpslope)]
  elbow.dt[,x:=(intercept-perpcept)/(perpslope - slope)]
  elbow.dt[,dist:=sqrt((value-y)^2 + (idx-x)^2)]

  maxidx<-which.max(elbow.dt$dist)
  elbow.dt[idx<maxidx]$selected<-as.factor('Selected')
  elbow.dt[,propvar:=cumsum(value)/sum(value)]
  elbow.dt
}

# Step 0: Get transcript-based counts and GTF data
setwd("/Users/rosal/OneDrive - University of Utah/2020/career/analyze/data/transcriptome-dimensions/")
#counts.gtf <- get_data("Homo_sapiens.GRCh37.74.gtf.gz","MMRF_CoMMpass_IA14a_E74GTF_Salmon_V7.2_Filtered_Transcript_Counts.txt.gz")
#save(counts.gtf,file="rdata/counts_gtf_20200402.RData")
load("rdata/counts_gtf_20200402.RData")

# Step 1: Quality control
## Get list of baseline samples
load("rdata/setup-clinical-data-20200402.rdata")
baseline <- clin.dt[collection_reason=="Baseline"]$sample_id
rm(clin.dt,key)
## Select baseline samples and remove low count genes and samples
qc <- quality_control(counts.gtf,baseline)
qc.counts <- counts.gtf[gene_name %in% qc$keep.genes] %>% select(-seqid,-gene_id,-gene_biotype)
qc.melt <- data.table::melt(qc.counts,
  id.vars=c("gene_name","name_n_transcripts","transcript_id","length"),
  variable.name="sample_id",
  value.name="count")

# Step 2: Normalize and truncate
norm <- process_transcripts(qc.melt)
## Convert to wide format of normalized values
all.dt<-dcast(norm$melt,sample_id+size_factor~gene_name,value.var='adjlogcpkmed')
## Annotate sample qc data
all.qc.dt<-data.table(inner_join(qc$samples,all.dt,by='sample_id'))

# Step 3: PCA on selected genes and samples
## select baseline samples with >100 reads in 90% of selected genes
base.dt <- all.qc.dt[sample_id %in% baseline & prop_gene_lt100 < 0.1][,-c(1:3,5)]
pca <- prcomp(base.dt[,-1,with=FALSE],center=TRUE,scale=FALSE,retx=TRUE)
pcvar <- data.table(pc=colnames(pca$x),value=pca$sdev^2)
elbow.dt <- elbow_finder(pcvar)
pc.scores <- cbind(base.dt[,1],pca$x[,elbow.dt[selected=='Selected']$pc])

## apply to all samples
all.pc.scores<-cbind(all.qc.dt[,c(1:5)],scale(data.matrix(all.qc.dt[,-c(1:5)]),
  center=colMeans(base.dt[,-1]),scale=FALSE)%*%pca$rotation[,1:31])

# now save these
save(pca,elbow.dt,pc.scores,all.qc.dt,all.pc.scores,file='rdata/process_20200402.RData')
