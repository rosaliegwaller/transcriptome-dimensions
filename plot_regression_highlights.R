#Plot highlights from dimension associations

#Packages
library("ggplot2")
library("survminer")
library("survival")

#Load data
setwd("/Users/rosal/OneDrive - University of Utah/2020/analyze/data/transcriptome-dimensions/")
load(file = "rdata/regression_models-20200429.RData") #regression results from "regression_overview_plot.R"
load(file='rdata/process_20200428.RData') #pcs and clinical data in dt
data <- merge(samples,score_sd,by='sample_id')
rm(samples,score_sd,pc,cbat.dt)

#Ploting setup
theme_set(theme_minimal() + theme(legend.position = "right"))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
PC <- colnames(dplyr::select(data,starts_with("PC")))

#ISS score grouped by dimenions ----
data$D_PT_iss = factor(data$D_PT_iss, levels = c(1,2,3), ordered = TRUE)
mdt <- data %>% dplyr::select("patient_id","D_PT_iss",ends_with("_sd")) %>% na.omit()
fit <- models$ISS

beta <- coef(fit)
RR <- exp(coef(fit))
CI <- exp(confint(fit))
colnames(CI) <- c("Lower", "Upper")
P <- pnorm(abs(coef(summary(fit))[-c(41,42),"t value"]),lower.tail = FALSE)*2
tb <- as.data.table(cbind(beta,RR, CI, P))
tb$PC <- PC
rm(beta,RR,CI,P)
setorder(tb,P)
tb[P<0.05]

ggplot(data = mdt,aes(x=PC3_sd,y=PC1_sd,color=D_PT_iss))+
  geom_point()

#combination score
pcs <- data.matrix(dplyr::select(data,starts_with("PC")))
bta <- coef(fit)
tmp = as.data.table(pcs %*% diag(bta)) #multiply all pcs by corresponding beta value
unique((bta[1]*pcs[,1])==tmp[,1])#check algorithm 
unique((bta[37]*pcs[,37])==tmp[,37])#check algorithm 
colnames(tmp) = PC
tmp$patient_id = data$patient_id
sig = dplyr::select(tmp,"patient_id",all_of(tb[tb$P<0.05,]$PC)) #select significant pcs in model
sig$comb = rowSums(sig[,-"patient_id"])
cmb = merge(sig[,c("patient_id","comb")],mdt,by="patient_id")
rm(fit,tbl,pcs,bta,tmp,sig)
cmb$PC = "combined"

ggboxplot(cmb,x="PC",y="comb",fill="D_PT_iss") +
  stat_compare_means(aes(group = D_PT_iss),method = "anova")#,label = "p.signif",hide.ns = TRUE)


#Kaplan-Meier curves of the weighted, combined PCs ----
#Overall survival
fit <- models$OS
tbl <- as.data.table(coef(fit$summary))
tbl$pc <- PC
pcs <- data.matrix(dplyr::select(data,starts_with("PC")))
bta <- tbl$coef
tmp = as.data.table(pcs %*% diag(bta)) #multiply all pcs by corresponding beta value
unique((bta[1]*pcs[,1])==tmp[,1])#check algorithm 
unique((bta[37]*pcs[,37])==tmp[,37])#check algorithm 
colnames(tmp) = PC
tmp$patient_id = data$patient_id
sig = dplyr::select(tmp,"patient_id",all_of(tbl[tbl$`Pr(>|z|)` < 0.05,]$pc)) #select significant pcs in model
sig$comb = rowSums(sig[,-"patient_id"])
cmb = merge(sig[,c("patient_id","comb")],data[,c("patient_id","ttcos","censos")],by="patient_id")
rm(fit,tbl,pcs,bta,tmp,sig)
cmb$censos = as.numeric(cmb$censos)

cmb$combt<-cut(cmb$comb,quantile(cmb$comb,c(0,1/3,2/3,1)),include.lowest = T) #Tertiles
ggsurvplot(survfit(Surv(time = cmb$ttcos, event = cmb$censos) ~ combt, data=cmb),palette=cbPalette[c(6,1,7)])

#PFS survival
fit <- models$PFS
tbl <- as.data.table(coef(fit$summary))
tbl$pc <- PC
pcs <- data.matrix(dplyr::select(data,starts_with("PC")))
bta <- tbl$coef
tmp = as.data.table(pcs %*% diag(bta)) #multiply all pcs by corresponding beta value
unique((bta[1]*pcs[,1])==tmp[,1])#check algorithm 
unique((bta[37]*pcs[,37])==tmp[,37])#check algorithm 
colnames(tmp) = PC
tmp$patient_id = data$patient_id
sig = dplyr::select(tmp,"patient_id",all_of(tbl[tbl$`Pr(>|z|)` < 0.05,]$pc)) #select significant pcs in model
sig$comb = rowSums(sig[,-"patient_id"])

cmb = merge(sig[,c("patient_id","comb")],data[,c("patient_id","ttcpfs","censpfs")],by="patient_id")
rm(fit,tbl,pcs,bta,tmp,sig)
cmb$censpfs = as.numeric(cmb$censpfs)

cmb$combt<-cut(cmb$comb,quantile(cmb$comb,c(0,1/3,2/3,1)),include.lowest = T) #Tertiles
ggsurvplot(survfit(Surv(time = cmb$ttcpfs, event = cmb$censpfs) ~ combt, data=cmb),palette=cbPalette[c(6,1,7)])

#Time to first-line treatment failure
fit <- models$TFTF
tbl <- as.data.table(coef(fit$summary))
tbl$pc <- PC
pcs <- data.matrix(dplyr::select(data,starts_with("PC")))
bta <- tbl$coef
tmp = as.data.table(pcs %*% diag(bta)) #multiply all pcs by corresponding beta value
unique((bta[1]*pcs[,1])==tmp[,1])#check algorithm 
unique((bta[37]*pcs[,37])==tmp[,37])#check algorithm 
colnames(tmp) = PC
tmp$patient_id = data$patient_id
sig = dplyr::select(tmp,"patient_id",all_of(tbl[tbl$`Pr(>|z|)` < 0.05,]$pc)) #select significant pcs in model
sig$comb = rowSums(sig[,-"patient_id"])

cmb = merge(sig[,c("patient_id","comb")],data[,c("patient_id","ttctf1","censtf1")],by="patient_id")
rm(fit,tbl,pcs,bta,tmp,sig)
cmb$censtf1 = as.numeric(cmb$censtf1)

cmb$combt<-cut(cmb$comb,quantile(cmb$comb,c(0,1/3,2/3,1)),include.lowest = T) #Tertiles
ggsurvplot(survfit(Surv(time = cmb$ttctf1, event = cmb$censtf1) ~ combt, data=cmb),palette=cbPalette[c(6,1,7)])

#Differences in PCs by race ----
#Create combination from regression model
fit <- models$race
tbl <- as.data.table(coef(summary(fit))[-1,])
tbl$pc <- colnames(dplyr::select(data[D_PT_race%in%c(1,2)],starts_with("PC")))
pcs <- data.matrix(dplyr::select(data[D_PT_race%in%c(1,2)],starts_with("PC")))
bta <- tbl$Estimate
tmp = as.data.table(pcs %*% diag(bta)) #multiply all pcs by corresponding beta value
unique((bta[1]*pcs[,1])==tmp[,1])#check algorithm 
unique((bta[37]*pcs[,37])==tmp[,37])#check algorithm 
colnames(tmp) = colnames(dplyr::select(data[D_PT_race%in%c(1,2)],starts_with("PC")))
tmp$patient_id = data[D_PT_race%in%c(1,2)]$patient_id
sig = dplyr::select(tmp,"patient_id",all_of(tbl[tbl$`Pr(>|z|)` < 0.05,]$pc)) #select significant pcs in model
sig$comb = rowSums(sig[,-"patient_id"])
rcmb = merge(sig[,c("patient_id","comb")],data[D_PT_race%in%c(1,2)][,c("patient_id","D_PT_race")],by="patient_id")
rcmb <- mutate(rcmb,race=if_else(D_PT_race==1,"white","black"))
rcmb$PC = "sig"
rm(fit,tbl,pcs,bta,tmp,sig)

#Plot race
race.melt <- melt(data = dplyr::select(data[D_PT_race%in%c(1,2)],"patient_id","D_PT_race",starts_with("PC")),id.vars = c("patient_id","D_PT_race"))
race.melt <- mutate(race.melt,race=if_else(D_PT_race==1,"white","black"))
race.melt <- mutate(race.melt,PC=gsub("PC","",gsub("_sd","",variable)))
#Box plot showing mean comparision
ggboxplot(race.melt[race.melt$PC==37,],x="PC",y="value",color="race",palette="jco") +
  stat_compare_means(aes(group = race),method = "anova")#label = "p.signif",hide.ns = TRUE,
ggboxplot(rcmb,x="PC",y="comb",fill="race",palett=cbPalette[c(3,8)]) +
  stat_compare_means(aes(group = race),method = "anova",label = "p.signif",hide.ns = TRUE)

#Density plot
ggplot(rcmb,aes(x=comb,fill=race)) +
  geom_density(alpha=0.3)

ggplot(race.melt[race.melt$PC==37,],aes(x=value,fill=race)) +
  geom_density(alpha=0.3)
