# ComBat on subset of clinical variables of interest

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

find_pcs <- function(DAT) {
  dt <- data.table(t(DAT)) #sample x gene data table
  dt$sample_id <- colnames(DAT) #annotate sample ids
  
  # run pca on combat data
  pca <- prcomp(dt[,-"sample_id",with=FALSE],center=TRUE,scale=FALSE,retx=TRUE)
  pcvar <- data.table(pc=colnames(pca$x),value=pca$sdev^2)
  elbow <- elbow_finder(pcvar)
  score <- cbind(dt[,"sample_id"],pca$x[,elbow[selected=='Selected']$pc])
  
  out <- list("score"=score,"elbow"=elbow)
  return(out)
}

find_covariates <- function(CLIN, Y1, Y2){
  DATA <- CLIN
  X <- "batch"
  
  # Continuous variables
  lmod <- expand.grid(y=Y1, x=X, stringsAsFactors = T) %>%
    mutate(formula = paste(y,"~",x)) %>%
    group_by(formula) %>%
    mutate(F = summary(lm(formula, data=DATA))$fstatistic[1]) %>%
    mutate(DF1 = summary(lm(formula, data=DATA))$fstatistic[2]) %>%
    mutate(DF2 = summary(lm(formula, data=DATA))$fstatistic[3]) %>%
    mutate(R2 = summary(lm(formula, data=DATA))$r.squared) %>%
    ungroup()
  lmod$P <- pf(lmod$F,lmod$DF1,lmod$DF2,lower.tail=F)
  setDT(lmod)
  lmod[,bon_sig:=as.factor('No')]
  lmod[P < 0.05/length(Y1)]$bon_sig <- as.factor('Yes')
  setorder(lmod,P)
  
  # Categorical variables
  DATA <- CLIN %>% dplyr::select(all_of(Y2))
  p <- 1
  for(col in DATA){
    p <- rbind(p,fisher.test(x=CLIN$batch,y=col,simulate.p.value=TRUE,B=5000)$p.value)
  }
  fisher.p <- data.table(cbind(Y2,p[-1,]))
  colnames(fisher.p) <- c("variable","p-value")

  out <- list("con"=lmod,"cat"=fisher.p)
  return(out)
}

batch_correction <- function(DAT, BATCH, MOD) {
  cbat <- ComBat(dat = DAT, batch = BATCH, mod = MOD) #run combat
  cbat.dt <- data.table(t(cbat)) #sample x gene data table
  cbat.dt$sample_id <- colnames(DAT) #annotate sample ids
  
  # run pca on combat data
  pca <- prcomp(cbat.dt[,-"sample_id",with=FALSE],center=TRUE,scale=FALSE,retx=TRUE)
  pcvar <- data.table(pc=colnames(pca$x),value=pca$sdev^2)
  elbow <- elbow_finder(pcvar)
  score <- cbind(cbat.dt[,"sample_id"],pca$x[,elbow[selected=='Selected']$pc])
  
  out <- list("score"=score,"elbow"=elbow)
  return(out)
}

find_association <- function(CLIN, SCORE){
  DATA <- merge(CLIN,SCORE)
  Y <- colnames(dplyr::select(SCORE,starts_with('PC')))
  X <- colnames(CLIN[,-"sample_id"])
  
  lmod <- expand.grid(y=Y, x=X, stringsAsFactors = T) %>%
    mutate(formula = paste(y,"~",x)) %>%
    group_by(formula) %>%
    mutate(F = summary(lm(formula, data=DATA))$fstatistic[1]) %>%
    mutate(DF1 = summary(lm(formula, data=DATA))$fstatistic[2]) %>%
    mutate(DF2 = summary(lm(formula, data=DATA))$fstatistic[3]) %>%
    mutate(R2 = summary(lm(formula, data=DATA))$r.squared) %>%
    ungroup()
  lmod$P <- pf(lmod$F,lmod$DF1,lmod$DF2,lower.tail=F)
  setDT(lmod)
  lmod[,bon_sig:=as.factor('No')]
  lmod[P < 0.05/length(Y)]$bon_sig <- as.factor('Yes')
  setorder(lmod,P)
  
  return(lmod)
}

# Data
setwd("/Users/rosal/OneDrive - University of Utah/2020/career/analyze/data/transcriptome-dimensions/")
load(file = "rdata/norm_data_20200406.rdata")
load(file = "rdata/clin_20200405.rdata")

dt <- merge(clin,base.dt,by="sample_id")
rm(clin,base.dt,all.dt)

DAT = t(dt[,-c(1:25)])
colnames(DAT) <- dt$sample_id
BATCH = as.numeric(dt$batch)

CLIN = dt[,c(1,3:25)]

# no batch correction
nbat <- find_pcs(DAT)
lmod <- find_association(CLIN,nbat$score)

# find sig associations of clin data and batch
batch_clin <- find_covariates(CLIN, 
                Y1 = colnames(dplyr::select(CLIN,"D_PT_age",starts_with("ttc"),starts_with("D_LAB"))),
                Y2 = colnames(dplyr::select(CLIN,starts_with("cens"),"D_PT_gender","D_PT_race","D_PT_ethnic","D_PT_iss",contains("ABNORMALITY")))
                )

# 0: no clinical data in combat covariate model
MOD = NULL
cbat0 <- batch_correction(DAT,BATCH,MOD)
lmod0 <- find_association(CLIN,cbat0$score)

# 1: some clinical data
MOD <- data.matrix(CLIN[,c("D_PT_age","D_PT_gender","ttcos","ttcpfs","ttctf1")])
cbat1 <- batch_correction(DAT,BATCH,MOD)
lmod1 <- find_association(CLIN,cbat1$score)

# 2: all clinical data
MOD <- data.matrix(CLIN[,-c("sample_id","batch")])
MOD[is.na(MOD)] <- -9
cbat2 <- batch_correction(DAT,BATCH,MOD)
lmod2 <- find_association(CLIN,cbat2$score)

# merge PC association models together
pc_lm <- full_join(
  full_join(lmod[,c("formula","F","P")],lmod0[,c("formula","F","P")],by="formula"),
  full_join(lmod1[,c("formula","F","P")],lmod2[,c("formula","F","P")],by="formula"),
  by="formula")
colnames(pc_lm) <- c("formula","nbat.F","nbat.P","nmod.F","nmod.P","smod.F","smod.P","amod.F","amod.P")
pc_lm <- dplyr::select(pc_lm,"formula",ends_with("F"),ends_with("P"))
View(pc_lm %>% select("formula",ends_with("P")))
