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

## Load normalized expression data
load(file = "rdata/norm_data_20200406.rdata")

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
DATA <- dt %>% dplyr::select(all_of(Y))
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

# COMBAT with all clinical variables of interest
load(file = "rdata/norm_data_20200406.rdata")
load(file = "rdata/clin_20200405.rdata")
dt <- merge(clin,base.dt,by="sample_id")
DAT = t(dt[,-c(1:25)])
colnames(DAT) <- dt$sample_id
MOD <- data.matrix(dt[,c(4:25)])
MOD[is.na(MOD)] <- -9
BATCH = as.numeric(dt$batch)

## run combat
cbat <- ComBat(dat = DAT, batch = BATCH, mod = MOD)

## transpose data
cbat.dt <- data.table(t(cbat))
cbat.dt$sample_id <- dt$sample_id
rm(DAT,MOD,BATCH,cbat,clin,base.dt,all.dt)

## run pca on combat corrected data
pca <- prcomp(cbat.dt[,-"sample_id",with=FALSE],center=TRUE,scale=FALSE,retx=TRUE)
pcvar <- data.table(pc=colnames(pca$x),value=pca$sdev^2)
elbow <- elbow_finder(pcvar)
score <- cbind(cbat.dt[,"sample_id"],
  pca$x[,elbow[selected=='Selected']$pc])

## test for batch sig associations with new pcs
DATA <- merge(dt[,c(1:25)],score)
Y <- colnames(dplyr::select(DATA,starts_with('PC')))
X <- "batch"

batch <- expand.grid(y=Y, x=X, stringsAsFactors = T) %>%
  mutate(formula = paste(y,"~",x)) %>%
  group_by(formula) %>%
  mutate(F = summary(lm(formula, data=DATA))$fstatistic[1]) %>%
  mutate(DF1 = summary(lm(formula, data=DATA))$fstatistic[2]) %>%
  mutate(DF2 = summary(lm(formula, data=DATA))$fstatistic[3]) %>%
  mutate(R2 = summary(lm(formula, data=DATA))$r.squared) %>%
  ungroup()
batch$P <- pf(batch$F,batch$DF1,batch$DF2,lower.tail=F)

# -----------------------------------------------------------------------------
# COMBAT with subset of clinical variables of interest
load(file = "rdata/norm_data_20200406.rdata")
load(file = "rdata/clin_20200405.rdata")
dt <- merge(clin,base.dt,by="sample_id")
DAT = t(dt[,-c(1:25)])
colnames(DAT) <- dt$sample_id
MOD <- data.matrix(dt[,c(4)])
MOD[is.na(MOD)] <- -9
BATCH = as.numeric(dt$batch)

## run combat
cbat <- ComBat(dat = DAT, batch = BATCH, mod = MOD)

## transpose data
cbat.dt <- data.table(t(cbat))
cbat.dt$sample_id <- dt$sample_id
rm(DAT,MOD,BATCH,cbat,clin,base.dt,all.dt)

## run pca on combat corrected data
pca <- prcomp(cbat.dt[,-"sample_id",with=FALSE],center=TRUE,scale=FALSE,retx=TRUE)
pcvar <- data.table(pc=colnames(pca$x),value=pca$sdev^2)
elbow <- elbow_finder(pcvar)
score <- cbind(cbat.dt[,"sample_id"],
  pca$x[,elbow[selected=='Selected']$pc])

## test for batch sig associations with new pcs
DATA <- merge(dt[,c(1:25)],score)
Y <- colnames(dplyr::select(DATA,starts_with('PC')))
X <- "batch"

batch <- expand.grid(y=Y, x=X, stringsAsFactors = T) %>%
  mutate(formula = paste(y,"~",x)) %>%
  group_by(formula) %>%
  mutate(F = summary(lm(formula, data=DATA))$fstatistic[1]) %>%
  mutate(DF1 = summary(lm(formula, data=DATA))$fstatistic[2]) %>%
  mutate(DF2 = summary(lm(formula, data=DATA))$fstatistic[3]) %>%
  mutate(R2 = summary(lm(formula, data=DATA))$r.squared) %>%
  ungroup()
batch$P <- pf(batch$F,batch$DF1,batch$DF2,lower.tail=F)
