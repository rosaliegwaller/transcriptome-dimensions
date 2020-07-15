#!/usr/bin/env Rscript

# Libraries
library(rtracklayer)
library(dplyr)
library(data.table)
library(sva)

# Functions ----
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
  tmp <- counts.gtf[seqid %in% c(1:22) & gene_biotype=='protein_coding'] %>% dplyr::select("gene_name",contains("MMRF"))
  # aggregate trasncript counts to gene_name counts
  gene <- data.table(aggregate(. ~ gene_name,data=tmp,FUN=sum))
  # find proportion of genes with <100 counts in baseline samples only
  base <- gene %>% dplyr::select("gene_name",all_of(baseline.samples))
  base$prop_base_lt100 <- rowSums(base[,-1] < 100) / length(baseline.samples)

  keep.genes <- base[prop_base_lt100 < 0.05]$gene_name

  # count all samples total reads and proportion of qc genes with > 100 reads
  total_reads <- gene[,-1][,lapply(.SD,sum)]
  total_reads_keep_genes <- gene[gene_name %in% keep.genes,-1][,lapply(.SD,sum)]
  prop_gene_lt100 <- gene[gene_name %in% keep.genes,-1][,lapply(.SD,function(x)sum(x<100)/length(keep.genes))]
  samples<-data.table(t((rbind(total_reads,total_reads_keep_genes,prop_gene_lt100))))
  colnames(samples)<-c("total_reads","total_reads_keep_genes","prop_gene_lt100")
  samples$sample_id <- colnames(total_reads)
  
  remove.samples <- samples[prop_gene_lt100>0.1]$sample_id

  out <- list("keep.genes"=keep.genes,"samples"=samples,"remove.samples"=remove.samples)
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

find_pcs <- function(cbat.dt) {
  # run pca on combat data
  pca <- prcomp(cbat.dt[,-"sample_id",with=FALSE],center=TRUE,scale=FALSE,retx=TRUE)
  pcvar <- data.table(pc=colnames(pca$x),value=pca$sdev^2)
  elbow <- elbow_finder(pcvar)
  score <- cbind(cbat.dt[,"sample_id"],pca$x[,elbow[selected=='Selected']$pc])
  
  out <- list("score"=score,"elbow"=elbow,"pca"=pca)
  return(out)
}

# Step 0: Get transcript-based counts and GTF data
setwd("/Users/rosal/OneDrive - University of Utah/2020/analyze/data/transcriptome-dimensions/")
#counts.gtf <- get_data("Homo_sapiens.GRCh37.74.gtf.gz","MMRF_CoMMpass_IA14a_E74GTF_Salmon_V7.2_Filtered_Transcript_Counts.txt.gz")
#save(counts.gtf,file="rdata/counts_gtf_20200402.RData")
load("rdata/counts_gtf_20200402.RData")

## Select samples
load(file = "rdata/clin_20200405.rdata") #load clinical data
bm.clin <- clin[grep("_BM",clin$sample_id),]
bm.clin$patient_id <- gsub("_5_BM","",gsub("_4_BM","",gsub("_3_BM","",gsub("_2_BM","",gsub("_1_BM","",gsub("MMRF_","",bm.clin$sample_id))))))
base.samples <- bm.clin[collection_reason=="Baseline"][,c('sample_id','patient_id','collection_reason')]
prog.samples <- bm.clin[patient_id%in%base.samples$patient_id & collection_reason=="Confirm Progression"][,c('sample_id','patient_id','collection_reason')]
tmp <- rbind(base.samples,prog.samples)
samples <- droplevels(merge(tmp,clin)) #baseline and progression samples with clinical data
rm(tmp,base.samples,bm.clin,clin,prog.samples)
baseline <- samples[collection_reason=="Baseline"]$sample_id #list of baseline samples

express <- counts.gtf %>% select(1:7,all_of(samples$sample_id)) #expression estimates for samples and transcript info
rm(counts.gtf)

# Step 1: Quality control ---------------
## Select baseline samples and remove low count genes and samples
qc <- quality_control(express,baseline)
qc.counts <- express[gene_name %in% qc$keep.genes] %>% dplyr::select(-seqid,-gene_id,-gene_biotype,-qc$remove.samples)
qc.melt <- data.table::melt(qc.counts,
  id.vars=c("gene_name","name_n_transcripts","transcript_id","length"),
  variable.name="sample_id",
  value.name="count")
rm(qc.counts)

# Step 2: Normalize and truncate ---------------
norm <- process_transcripts(qc.melt)
rm(qc.melt)

## Convert to wide format of normalized values
all.dt <- dcast(norm$melt,sample_id+size_factor~gene_name,value.var='adjlogcpkmed')
norm.dt <- merge(samples,all.dt)
rm(all.dt)

# Step 3: Batch correction ---------------
DAT = t(norm.dt[,-c(1:27)])
colnames(DAT) <- norm.dt$sample_id
BATCH = as.numeric(norm.dt$batch)

CLIN = norm.dt[,c(1:27)]
sapply(CLIN, function(x) sum(is.na(x))) #count missing
MOD <- data.matrix(CLIN[,c("D_PT_age","D_PT_gender","ttcos","censos","ttcpfs","censpfs","ttctf1","censtf1")])
cbat <- ComBat(dat = DAT, batch = BATCH, mod = MOD) #run combat

## convert combat output to data table
cbat.dt <- data.table(t(cbat)) #sample x gene data table
cbat.dt$sample_id <- colnames(cbat) #annotate sample ids
rm(DAT,BATCH,MOD,cbat,CLIN)

# Step 4: PCA on selected genes and samples - batch corrected ---------------
## select baseline samples with >100 reads in 90% of selected genes
pc <- find_pcs(cbat.dt[sample_id %in% baseline])
#View(pc$elbow)

# standardize pc scores ----
dt <- pc$score
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
dt$PC32_sd <- dt$PC32/sd(dt$PC32)
dt$PC33_sd <- dt$PC33/sd(dt$PC33)
dt$PC34_sd <- dt$PC34/sd(dt$PC34)
dt$PC35_sd <- dt$PC35/sd(dt$PC35)
dt$PC36_sd <- dt$PC36/sd(dt$PC36)
dt$PC37_sd <- dt$PC37/sd(dt$PC37)
dt$PC38_sd <- dt$PC38/sd(dt$PC38)
dt$PC39_sd <- dt$PC39/sd(dt$PC39)
dt$PC40_sd <- dt$PC40/sd(dt$PC40)
#dt$PC41_sd <- dt$PC41/sd(dt$PC41)
#dt$PC42_sd <- dt$PC42/sd(dt$PC42)
#dt$PC43_sd <- dt$PC43/sd(dt$PC43)
#dt$PC44_sd <- dt$PC44/sd(dt$PC44)
#dt$PC45_sd <- dt$PC45/sd(dt$PC45)

score_sd <- dt %>% dplyr::select("sample_id",ends_with("_sd"))

# now save these
save(samples,score_sd,pc,cbat.dt,file='rdata/process_20200428.RData')
