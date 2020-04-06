#!/usr/bin/env Rscript

# Libraries
library(dplyr)
library(data.table)
library(sva)

# Functions
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

# Process
setwd("/Users/rosal/OneDrive - University of Utah/2020/career/analyze/data/transcriptome-dimensions/")
## Load PC data
load("rdata/process_20200402.RData")

## Load batch data
#load("rdata/setup-clinical-data-20200402.rdata")
#batch <- clin.dt %>% select("public_id","sample_id","collection_reason","batch")
#rm(clin.dt,key)
load(file = "rdata/clin_20200405.rdata")

## Merge and remove extra objects
dt <- merge(clin,pc.scores)
rm(all.pc.scores,elbow.dt,pc.scores,pca)

# Logistic regression
Y <- colnames(dt[,c(4,8,10,12,15:18,26:56)])
X <- "batch"
DATA <- dt

lm_batch <- expand.grid(y=Y, x=X, stringsAsFactors = T) %>%
  mutate(formula = paste(y,"~",x)) %>%
  group_by(formula) %>%
  mutate(F = summary(lm(formula, data=DATA))$fstatistic[1]) %>%
  mutate(DF1 = summary(lm(formula, data=DATA))$fstatistic[2]) %>%
  mutate(DF2 = summary(lm(formula, data=DATA))$fstatistic[3]) %>%
  mutate(R2 = summary(lm(formula, data=DATA))$r.squared) %>%
  ungroup()
lm_batch$P <- pf(lm_batch$F,lm_batch$DF1,lm_batch$DF2,lower.tail=F)
rm(DATA,X,Y)

# Categorical variables
Y <- colnames(dt[,c(5:7,9,11,13:14,19:25)])
DATA <- select(dt,Y)
p <- 1
for(col in DATA){
  p <- rbind(p,fisher.test(x=dt$batch,y=col,simulate.p.value=TRUE,B=5000)$p.value)
}
fisher.p <- data.table(cbind(Y,p[-1,]))
colnames(fisher.p) <- c("variable","p-value")
rm(col,p,DATA,Y)

# Remove batch effects
## transpose expression data to genes x samples
## select samples with clinical data
## remove non-measure values
## adjusted and normalized values
#DAT = data.table(t(all.qc.dt[sample_id %in% clin$sample_id][,-c(1:5)]))
#rownames(DAT) <- colnames(all.qc.dt)[-c(1:5)]
DAT = t(all.qc.dt[sample_id %in% clin$sample_id][,-c(1:5)])
colnames(DAT) <- clin$sample_id

## select clinical variables that are significantly associated with batch
sb <- c("ttcos","censos","ttcpfs","censpfs","ttctf1","censtf1","D_PT_age")
MOD <- as.numeric(as.matrix(select(clin,all_of(sb))))

BATCH = as.numeric(clin$batch)

## run combat
dat_com <- ComBat(dat = DAT, batch = clin$batch, mod = MOD)

## transpose data
all.combat <- data.table(t(dat_com))
all.combat$sample_id <- clin$sample_id
#all.combat$prop_gene_lt100 <- clin$prop_gene_lt100

# run pca on combat corrected data
baseline <- clin[collection_reason=="Baseline"]$sample_id
base.combat <- all.combat[sample_id %in% baseline]
pca.combat <- prcomp(base.combat[,-"sample_id",with=FALSE],center=TRUE,scale=FALSE,retx=TRUE)
pcvar.combat <- data.table(pc=colnames(pca.combat$x),value=pca.combat$sdev^2)
elbow.combat <- elbow_finder(pcvar.combat)
pc.scores.combat <- cbind(base.combat[,"sample_id"],
  pca.combat$x[,elbow.combat[selected=='Selected']$pc])

## apply to all samples
all.pc.scores.combat <- cbind(all.combat[,"sample_id"],
  scale(data.matrix(all.combat[,-"sample_id"]),
  center=colMeans(base.combat[,-"sample_id"]),
  scale=FALSE)%*%pca.combat$rotation[,1:37])

## and save
save(all.combat,base.combat,
  pca.combat,pcvar.combat,elbow.combat,
  pc.scores.combat,all.pc.scores.combat,
  file='rdata/batch_pca_20200405.RData')


# test for batch sig associations with new pcs
## Merge and remove extra objects
dt2 <- merge(clin,pc.scores.combat)

# Logistic regression
Y <- colnames(dt2[,c(26:62)])
X <- "batch"
DATA <- dt2

pc2.batch <- expand.grid(y=Y, x=X, stringsAsFactors = T) %>%
  mutate(formula = paste(y,"~",x)) %>%
  group_by(formula) %>%
  mutate(F = summary(lm(formula, data=DATA))$fstatistic[1]) %>%
  mutate(DF1 = summary(lm(formula, data=DATA))$fstatistic[2]) %>%
  mutate(DF2 = summary(lm(formula, data=DATA))$fstatistic[3]) %>%
  mutate(R2 = summary(lm(formula, data=DATA))$r.squared) %>%
  ungroup()
pc2.batch$P <- pf(pc2.batch$F,pc2.batch$DF1,pc2.batch$DF2,lower.tail=F)
rm(DATA,X,Y)


### OLD ###

# D_PT_gender
chisq.test(x = dt$batch, y = dt$D_PT_gender)
fisher.test(x = dt$batch, y = dt$D_PT_gender, simulate.p.value=TRUE)

# D_PT_race
chisq.test(x = dt$batch, y = dt$D_PT_race)
fisher.test(x = dt$batch, y = dt$D_PT_race, simulate.p.value=TRUE)

# D_PT_ethnic
chisq.test(x = dt$batch, y = dt$D_PT_ethnic)
fisher.test(x = dt$batch, y = dt$D_PT_ethnic, simulate.p.value=TRUE)

# censos
chisq.test(x = dt$batch, y = dt$censos)
fisher.test(x = dt$batch, y = dt$censos, simulate.p.value=TRUE)

# censpfs
chisq.test(x = dt$batch, y = dt$censpfs)
fisher.test(x = dt$batch, y = dt$censpfs, simulate.p.value=TRUE)

# censtf1
chisq.test(x = dt$batch, y = dt$censtf1)
fisher.test(x = dt$batch, y = dt$censtf1, simulate.p.value=TRUE)

# D_PT_iss
chisq.test(x = dt$batch, y = dt$D_PT_iss)
fisher.test(x = dt$batch, y = dt$D_PT_iss, simulate.p.value=TRUE)

# D_TRI_CF_ABNORMALITYPR6
chisq.test(x = dt$batch, y = dt$D_TRI_CF_ABNORMALITYPR6)
fisher.test(x = dt$batch, y = dt$D_TRI_CF_ABNORMALITYPR6, simulate.p.value=TRUE)

# D_TRI_CF_ABNORMALITYPR4
chisq.test(x = dt$batch, y = dt$D_TRI_CF_ABNORMALITYPR4)
fisher.test(x = dt$batch, y = dt$D_TRI_CF_ABNORMALITYPR4, simulate.p.value=TRUE)

# D_TRI_CF_ABNORMALITYPR3
chisq.test(x = dt$batch, y = dt$D_TRI_CF_ABNORMALITYPR3)
fisher.test(x = dt$batch, y = dt$D_TRI_CF_ABNORMALITYPR3, simulate.p.value=TRUE)

# D_TRI_CF_ABNORMALITYPR13
chisq.test(x = dt$batch, y = dt$D_TRI_CF_ABNORMALITYPR13)
fisher.test(x = dt$batch, y = dt$D_TRI_CF_ABNORMALITYPR13, simulate.p.value=TRUE)

# D_TRI_CF_ABNORMALITYPR8
chisq.test(x = dt$batch, y = dt$D_TRI_CF_ABNORMALITYPR8)
fisher.test(x = dt$batch, y = dt$D_TRI_CF_ABNORMALITYPR8, simulate.p.value=TRUE)

# D_TRI_CF_ABNORMALITYPR9
chisq.test(x = dt$batch, y = dt$D_TRI_CF_ABNORMALITYPR9)
fisher.test(x = dt$batch, y = dt$D_TRI_CF_ABNORMALITYPR9, simulate.p.value=TRUE)

# D_TRI_CF_ABNORMALITYPR11
chisq.test(x = dt$batch, y = dt$D_TRI_CF_ABNORMALITYPR11)
fisher.test(x = dt$batch, y = dt$D_TRI_CF_ABNORMALITYPR11, simulate.p.value=TRUE)




ct <- as.data.frame.matrix(table(dt$batch,dt$D_PT_race))
ctz <- setDT(ct)[rowSums(ct) > 0]
chisq.test(ctz)


Y <- colnames(dt[,c(5:7,9,11,13:14,19:25)])
X <- "batch"
DATA <- select(dt,Y)
output <- matrix(ncol=4, nrow=length(Y))
row <- 1

for(col in DATA){
  ct <- as.data.frame.matrix(table(dt$batch,col))
  ctz <- setDT(ct)[rowSums(ct) > 0]
  cs <- chisq.test(ctz)
  print(c(col,cs$statistic,cs$parameter,cs$p.value))
  #row = row + 1
}

x_batch <- expand.grid(y=Y, x=X) %>%
  ct <- as.data.frame.matrix(table(dt$Y, dt$X)) %>%
  ctz <- setDT(ct)[rowSums(ct) > 0]
  mutate(x <- chisq.test(ctz)$statistic) %>%
  mutate(df <- chisq.test(ctz)$parameter) %>%
  mutate(p <- chisq.test(ctz)$p.value)
