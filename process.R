#!/usr/bin/env Rscript

# Libraries
library(rtracklayer)
library(dplyr)
library(data.table)

# Functions
get_data <- function(gtf.path,counts.path){
  gtf <- readGFF(filepath=gtf.path) #read in GTF
  gtf$TRANSCRIPT_LENGTH <- gtf$end-gtf$start+1 #compute length of each feature
  #select protein coding exons
  tmp <- subset(gtf,gene_biotype=='protein_coding' & type=='exon',)
  #sum exon lengths for each transcript
  sum_exons<-aggregate(TRANSCRIPT_LENGTH~seqid+gene_name+gene_id+transcript_id,
    data = tmp, FUN=sum)
  colnames(sum_exons)<-c("CHR","GENE_NAME","GENE_ID",
                          "TRANSCRIPT_ID","TRANSCRIPT_LENGTH")
  #count number of transcripts for each gene name
  n_name <- sum_exons %>% count(GENE_NAME)
  colnames(n_name)[2] <- "GENE_NAME_N_TRANSCRIPTS"
  #count number of transcripts for each ensemble gene id
  n_id <- sum_exons %>% count(GENE_ID)
  colnames(n_id)[2] <- "GENE_ID_N_TRANSCRIPTS"
  #join
  gtf.dt <- data.table(
    inner_join(
      inner_join(n_name,sum_exons,by='GENE_NAME'),
      n_id,by='GENE_ID')
    )
  gtf.dt <- gtf.dt[,c("GENE_NAME","GENE_NAME_N_TRANSCRIPTS",
                  "GENE_ID","GENE_ID_N_TRANSCRIPTS",
                  "TRANSCRIPT_ID","TRANSCRIPT_LENGTH")]
  # read in rnaseq transcript counts
  counts <- read.csv(file = counts.path,header = T,sep = "\t") #read in counts
  counts.dt <- data.table(counts)
  # merge with gtf data
  counts.gtf <- inner_join(gtf.dt,counts.dt,by="TRANSCRIPT_ID")
  setDT(counts.gtf)
  samples <- colnames(counts.gtf)[-1:-6] #get list of sample ids
  molten.counts<-data.table::melt(counts.gtf,
    id.vars=c("GENE_NAME","GENE_ID","TRANSCRIPT_ID"),
    measure.vars=list("GENE_NAME_N_TRANSCRIPTS","GENE_ID_N_TRANSCRIPTS",
      "TRANSCRIPT_LENGTH",samples),
    variable.name="SAMPLE_ID",
    value.name=c("GENE_NAME_N_TRANSCRIPTS","GENE_ID_N_TRANSCRIPTS",
      "TRANSCRIPT_LENGTH","SAMPLE_COUNT"))
  setattr(molten.counts[["SAMPLE_ID"]],"levels",samples)
  return(molten.counts)
}

process_transcripts<-function(qc.counts) {
  setDT(qc.counts)
  temp.dt<-[,total_raw_counts:=sum(SAMPLE_COUNT),by=c('SAMPLE_ID','GENE_NAME')]
  temp.dt[,kb_length:=TRANSCRIPT_LENGTH/1000]
  final.dt<-temp.dt[,list(cpk=sum((SAMPLE_COUNT+1/GENE_NAME_N_TRANSCRIPTS)/kb_length)),by=c('SAMPLE_ID','GENE_NAME','total_raw_counts')]

  final.dt[,size_factor:=median(cpk),by='SAMPLE_ID']
  final.dt[,cpkmed:=cpk/size_factor]
  final.dt[,logcpkmed:=log2(cpkmed)]
  return(final.dt)
}

# Step 0: Get transcript-based counts and GTF data
setwd("/Users/rosal/OneDrive - University of Utah/2020/career/analyze/data/transcriptome-dimensions/")
#counts.melt<-get_data("Homo_sapiens.GRCh37.74.gtf.gz","MMRF_CoMMpass_IA14a_E74GTF_Salmon_V7.2_Filtered_Transcript_Counts.txt.gz")
load("rdata/counts.melt")

# Step 1: Quality control
## Select baseline samples and remove low count genes in baseline samples
low.genes<-counts.melt[SAMPLE_ID %like% "_1_BM",list(GENE_COUNT=sum(SAMPLE_COUNT)),by=c('GENE_NAME','SAMPLE_ID')][,quantile(GENE_COUNT,0.05,type=3),by='GENE_NAME'][V1<100]$GENE_NAME
qc.melt<-counts.melt[SAMPLE_ID%like%"_1_BM"&GENE_NAME%in%low.genes]

# Step 2: Normalize and truncate
med.molten<-process_transcripts(qc.melt)
med.molten[,mean:=mean(logcpkmed),by='GENE_NAME']
med.molten[,sd:=sd(logcpkmed),by='GENE_NAME']
med.molten[,adjlogcpkmed:=logcpkmed]
med.molten[(logcpkmed-mean)/sd>=5,adjlogcpkmed:=mean+5*sd]
med.molten[(logcpkmed-mean)/sd<= -5,adjlogcpkmed:=mean-5*sd]

# Sample QC?? -------------
sample.qc.dt<-merge(med.molten[total_raw_counts>=100,list(prop100plus=.N/10188),by='sample_id'],
      med.molten[,list(total_reads=sum(total_raw_counts)),by='sample_id'],by='sample_id')
sample.qc.dt$replicate<-as.factor('no')
sample.qc.dt[sample_id%in%a5581.rsem.iso.rep.molten$sample_id,replicate:=as.factor('yes')]

# Step 3: PCA on all features -----
temp1.dt<-dcast(med.molten,sample_id+size_factor~geneid,value.var='adjlogcpkmed')
pca<-prcomp(temp1.dt[,-c(1:2),with=FALSE],center=TRUE,scale=FALSE,retx=TRUE)
pcvar<-data.table(pc=colnames(pca$x),value=pca$sdev^2)
elbow.dt<-elbow_finder(pcvar)
pc.scores<-cbind(temp1.dt[,c(1:2)],pca$x[,elbow.dt[selected=='Selected']$pc])
final.dt<-merge(tissue.map,pc.scores,by='sample_id')

# Apply dimensions to tumor samples
all.molten<-process_rsem_isoforms(rsem.molten,good.genes)
temp1.dt<-dcast(all.molten,sample_id~geneid,value.var='logcpkmed')
all.pc.scores<-cbind(temp1.dt[,1],scale(data.matrix(temp1.dt[,-1]),center=colMeans(temp1.dt[,-1]),scale=FALSE)%*%pca$rotation[,1:12])
all.final.dt<-merge(tissue.map,all.pc.scores,by='sample_id')

# now save these
save(pca,elbow.dt,final.dt,all.final.dt,good.genes,file='rsem_iso_median_normalized_rep_summed_20200304.RData')





# working

qc.counts <- counts.gtf[!GENE_NAME%in%low.genes] %>% dplyr::select(contains("GENE_NAME") | contains("TRANSCRIPT_") | contains("_1_BM"))

norm.dt<-process_transcripts(qc.counts)

samples <- colnames(qc.counts )[-1:-4] #get list of sample ids
qc.melt<-data.table::melt(qc.counts,
  id.vars=c("GENE_NAME","TRANSCRIPT_ID"),
  measure.vars=list("GENE_NAME_N_TRANSCRIPTS","TRANSCRIPT_LENGTH",samples),
  variable.name="SAMPLE_ID",
  value.name=c("GENE_NAME_N_TRANSCRIPTS","TRANSCRIPT_LENGTH","SAMPLE_COUNT"))
setattr(qc.melt[["SAMPLE_ID"]],"levels",samples)

dt.seq[,lapply(.SD,function(x)sum(is.na(x)))],

keep.genes<-
keep.samples<-

## method to aggregate in long format
#dt2 <- counts.gtf[!GENE_NAME%in%low.genes] %>% dplyr::select(contains("GENE_NAME") | contains("_1_BM"))
#dt.gene <- aggregate(.~GENE_NAME+GENE_NAME_N_TRANSCRIPTS,data=dt2,FUN=sum)

#qc.counts<-molten.counts[SAMPLE_ID%like%"_1_BM" & !GENE_NAME%in%low.genes] #na bug
