#!/usr/bin/env Rscript

# Libraries
library(rtracklayer)
library(dplyr)
library(data.table)

# Functions
process_gtf<-function(gtf.path){
  gtf <- readGFF(filepath=gtf.path) #read in GTF
  gtf$length <- gtf$end - gtf$start + 1 #compute length of each feature
  tmp <- subset(gtf,gene_biotype=='protein_coding' & type=='exon',) #select protein coding exons
  sum_exons <- aggregate(length ~ gene_id + gene_name + transcript_id, data = tmp, FUN=sum) #sum exon lengths for each transcript
  n_transcripts <- sum_exons %>% count(gene_id) #count number of transcripts for each gene
  gtf.dt <- data.table(merge(n_transcripts,sum_exons,by='gene_id')) #merge tables
  return(gtf.dt)
}

process_counts<-function(counts.path){
  counts <- read.csv(file = counts.path,header = T,sep = "\t") #read in feature counts
  counts.dt <- data.table(counts)
  return(counts.dt)
}

select_baseline<-function(counts.dt){
  baseline.dt <- select(counts.dt,,c(TRANSCRIPT_ID,contains("_1_BM"))) #select baseline samples
  return(baseline.dt)
}

## TO DO - EDIT THIS FUNCTION
select_features<-function(gtf.dt,baseline.dt){
  genes <- select(gtf.dt,,c(transcript_id,gene_id))
  tmp <- merge(genes,baseline.dt,by.x='transcript_id',by.y='TRANSCRIPT_ID')
  samples <- colnames(tmp)[-1:-2]
  gene.count <- aggregate(samples~gene_id,data=tmp,FUN=sum)

  aggregate(data)

  tmp %>%

  aggregate(patterns("^MMRF")

  tmp.melt <- melt(tmp)


  tmp.dt$count <- tmp.dt %>% dplyr::select(contains("MMRF")) %>% rowSums()
  gene.count <- aggregate(count ~ gene_id, data = tmp.dt, FUN=sum)
  good.genes <- gene.count %>% dplyr::filter(count >= 100)

  features$N_LT_100 <- rowSums(tmp < 100)
  features$PROP_LT_100 <- features$N_LT_100 / ncol(tmp)

  counts_LT5_LT100 <- counts[features$PROP_LT_100 < 0.05,]
  return()
}

# Step 0: get data - takes some time
## feature annotation
gtf.dt <- process_gtf("/Users/rosal/OneDrive - University of Utah/2020/career/analyze/data/transcriptome-dimensions/Homo_sapiens.GRCh37.74.gtf.gz")

## feature counts from Salmon
counts.dt <- process_counts("/Users/rosal/OneDrive - University of Utah/2020/career/analyze/data/transcriptome-dimensions/MMRF_CoMMpass_IA14a_E74GTF_Salmon_V7.2_Filtered_Transcript_Counts.txt.gz")

molten.counts <- melt(counts.dt,id.vars="TRANSCRIPT_ID",variable.name="SAMPLE_ID",value.name="COUNT")
molten.data <- dplyr::inner_join(gtf.df,molten.counts,by.x="transcript_id",by.y="TRANSCRIPT_ID")

# Step 1: quality control
baseline.dt <- select_baseline(counts) #select baseline samples
coding <- select_features() #select coding genes with >5% of samples >100 reads



med.molten[,sd:=sd(logcpkmed),by='geneid']

tmp$N_Counts <- colSums(tmp) #count total reads
dt <- data.table(tmp[tmp$N_Counts < 10000000,] #remove sample if reads < 10M
