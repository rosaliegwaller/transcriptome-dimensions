#Disparities

#Load packages
library("data.table")
library("ggpubr")
library("tidyverse")
library("rstatix")
library("gridExtra")

#Functions
find_lm <- function(data,y,x){
  mdl <- expand.grid(y=Y, x=X, stringsAsFactors = T) %>%
    mutate(formula = paste(y,"~",x)) %>%
    group_by(formula) %>%
    mutate(F = summary(lm(formula, data=D))$fstatistic[1]) %>%
    mutate(DF1 = summary(lm(formula, data=D))$fstatistic[2]) %>%
    mutate(DF2 = summary(lm(formula, data=D))$fstatistic[3]) %>%
    mutate(R2 = summary(lm(formula, data=D))$r.squared) %>%
    ungroup()
  mdl$P <- pf(mdl$F,mdl$DF1,mdl$DF2,lower.tail=F)
  
  setorder(mdl,P)
  setDT(mdl)
  return(mdl)
}

#Inport data including pcs, race, ethnicity, gender, and age at diagnosis
setwd("/Users/rosal/OneDrive - University of Utah/2020/analyze/data/transcriptome-dimensions/")
load(file='rdata/process_20200428.RData') #pcs and clinical data in dt
data <- merge(samples[,c("sample_id","patient_id","batch","D_PT_age","D_PT_gender","D_PT_race","D_PT_ethnic")],score_sd)
rm(cbat.dt,pc,samples,score_sd)
load(file = "rdata/regression_models-20200429.RData") #regression results from "regression_overview_plot.R"

summary(data)

#Race
D <- droplevels(data[D_PT_race%in%c(1,2)] %>% dplyr::select("patient_id","D_PT_race",ends_with("_sd")))
Y <- colnames(dplyr::select(D,starts_with('PC')))
X <- "D_PT_race"
race <- find_lm(D,Y,X)
race[P<0.05/40]

#Plot race
race.melt <- melt(data = D,id.vars = c("patient_id","D_PT_race"))
race.melt <- mutate(race.melt,race=if_else(D_PT_race==1,"white","black"))
race.melt <- mutate(race.melt,PC=gsub("PC","",gsub("_sd","",variable)))
#Box plot showing mean comparision
race.p <- ggboxplot(race.melt[race.melt$PC==37,],x="PC",y="value",color="race",palette="jco") +
          stat_compare_means(aes(group = race),method = "anova")#label = "p.signif",hide.ns = TRUE,

#Ethnic
D <- droplevels(data[D_PT_ethnic%in%c(1,2)] %>% dplyr::select("sample_id","D_PT_ethnic",ends_with("_sd")))
Y <- colnames(dplyr::select(D,starts_with('PC')))
X <- "D_PT_ethnic"
ethc <- find_lm(D,Y,X)
ethc[P<0.05/40]

#Plot ethnic
ethc.melt <- melt(data = D,id.vars = c("sample_id","D_PT_ethnic"))
ethc.melt <- mutate(ethc.melt,PC=gsub("PC","",gsub("_sd","",variable)))
#Box plot showing mean comparision
ethc.p <- ggboxplot(ethc.melt[ethc.melt$PC==16,],x="PC",y="value",color="D_PT_ethnic",palette="jco") +
  stat_compare_means(aes(group = D_PT_ethnic),method = "anova",label = "p.format")
ethc.p

#Gender
D <- droplevels(data[D_PT_gender%in%c(1,2)] %>% dplyr::select("sample_id","D_PT_gender",ends_with("_sd")))
Y <- colnames(dplyr::select(D,starts_with('PC')))
X <- "D_PT_gender"
gndr <- find_lm(D,Y,X)
gndr[P<0.05/40]

#Plot gender
gndr.melt <- melt(data = D,id.vars = c("sample_id","D_PT_gender"))
gndr.melt <- mutate(gndr.melt,PC=gsub("PC","",gsub("_sd","",variable)))
#Box plot showing mean comparision
gndr.p <- ggboxplot(gndr.melt[gndr.melt$PC%in%c(39,3,11,9),],x="PC",y="value",color="D_PT_gender",palette="jco") +
  stat_compare_means(aes(group = D_PT_gender),method = "anova",label = "p.format")
gndr.p

#Age
D <- droplevels(data %>% dplyr::select("sample_id","D_PT_age",ends_with("_sd")))
Y <- colnames(dplyr::select(D,starts_with('PC')))
X <- "D_PT_age"
agdx <- find_lm(D,Y,X)
agdx[P<0.05/40]

D <- mutate(D,onset=if_else(D_PT_age<60,"dx<60yrs","dx>59yrs"))
Y <- colnames(dplyr::select(D,starts_with('PC')))
X <- "onset"
onst <- find_lm(D,Y,X)
onst[P<0.05/40]

#Plot agdx
agdx.melt <- melt(data = D,id.vars = c("sample_id","D_PT_age","onset"))
#agdx.melt <- mutate(agdx.melt,agdx=if_else(D_PT_age<60,"early-dx","standard"))
agdx.melt <- mutate(agdx.melt,PC=gsub("PC","",gsub("_sd","",variable)))
#Box plot showing mean comparision
agdx.p <- ggboxplot(agdx.melt[agdx.melt$PC%in%c(1,31),],x="PC",y="value",color="onset",palette="jco") +
  stat_compare_means(aes(group = onset),method = "anova",label = "p.format")
agdx.p

#Plot all together
grid.arrange(race.p,ethc.p,gndr.p,agdx.p,nrow=2)
