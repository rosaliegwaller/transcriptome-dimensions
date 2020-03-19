#!/usr/bin/env Rscript

## collect arguments
#args <- commandArgs(trailingOnly = TRUE)

## parse arguments (we expect the form --arg=value)
#parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
#argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
#names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
#args <- argsL
#rm(argsL)

## give some value to options if not provided 
#if(is.null(args$opt_arg1)) {args$opt_arg1="default_option1"}
#if(is.null(args$opt_arg2)) {args$opt_arg2="default_option1"} else {args$opt_arg2=as.numeric(args$opt_arg2)}

## default setting when no all arguments passed or help needed
#if("--help" %in% args | is.null(args$arg1) | is.null(args$arg2)) {
#	cat("The R Script arguments_section.R
#	    
#	    Mandatory arguments:
#	    --arg1=type           - description
#	    --arg2=type           - description
#	    --help                - print this text
#	   
#	    Optionnal arguments:
#	    --opt_arg1=String          - example:an absolute path, default:default_option1
#	    --opt_arg2=Value           - example:a threshold, default:10
#	    WARNING : here put all the things the user has to know
#	   
#	    Example:
#	    ./arguments_section.R --arg1=~/Documents/ --arg2=10 --opt_arg2=8 \n\n")
#	q(save="no")
#}
#
#cat("first mandatory argument : ", args$arg1,"\n",sep="")
#cat("second mandatory argument : ", args$arg2,"\n",sep="")
#cat("first optional argument : ", args$opt_arg1,"\n",sep="")
#cat("second optional argument : ", args$opt_arg2,"\n",sep="")

setwd("/home/rose/data/transcriptome-dimensions")

# setup
library(dplyr)

# load data
load("MMRF_CoMMpass_IA14a_E74GTF_Salmon_V7.2_Filtered_Transcript_Counts_Coding.rdata")

# remove feature if >5% samples <100 reads
features <- as.data.frame(counts[,1])
colnames(features) <- "TRANSCRIPT_ID"

tmp <- counts %>% dplyr::select(contains("_1_BM"))
rownames(tmp) <- features$TRANSCRIPT_ID

features$N_LT_100 <- rowSums(tmp < 100)
features$PROP_LT_100 <- features$N_LT_100 / ncol(tmp)

counts_LT5_LT100 <- counts[features$PROP_LT_100 < 0.05,]
save(counts_LT5_LT100,file="MMRF_CoMMpass_IA14a_E74GTF_Salmon_V7.2_Filtered_Transcript_Counts_Coding_100.rdata")
rm(tmp)

# remove sample if >10% remaining transcripts have <100 counts
tmp <- counts_LT5_LT100 %>% dplyr::select(contains("MMRF"))
samples <- as.data.frame(colnames(tmp))
colnames(samples) <- "MMRF_ID"

tmp <- as.data.frame(t(tmp))
colnames(tmp) <- counts_LT5_LT100$TRANSCRIPT_ID

samples$N_LT_100 <- rowSums(tmp < 100)
samples$P_LT_100 <- samples$N_LT_100 / ncol(tmp)

rm(tmp)

keep <- subset(samples,P_LT_100<=.1,"MMRF_ID")
keep <- as.character(keep$MMRF_ID)

# select samples and features and save data
df <- counts_LT5_LT100[,c("TRANSCRIPT_ID",keep)]
df <- df %>% dplyr::select(contains("_1_BM"))

save(df,file="MMRF_CoMMpass_IA14a_E74GTF_Salmon_V7.2_Filtered_Transcript_Counts_Coding_QC_Baseline.rdata")

# TO DO: repeat remove features with samples removed
# TO DO: gene-based QC for samples...
