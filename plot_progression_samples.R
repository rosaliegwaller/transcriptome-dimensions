#Plot follow-up samples

#Load package
library("dplyr")
library("ggplot2")
library("data.table")
library("wesanderson")

#Read data file of pc scores 
setwd("/Users/rosal/OneDrive - University of Utah/2020/analyze/data/transcriptome-dimensions/")
load(file='rdata/process_20200428.RData') #baseline pc scores and ComBat corrected baseline and progression expression estimates

baseline = samples[collection_reason=="Baseline"]$sample_id
PC = colnames(score_sd)[-1]

#Extract clinical data
vis = data.table(read.csv("CoMMpass_IA14_FlatFiles/MMRF_CoMMpass_IA14_PER_PATIENT_VISIT.csv"))
key.vis = data.table(read.csv("CoMMpass_IA14_FlatFile_Dictionaries/MMRF_CoMMpass_IA14_PER_PATIENT_VISIT.csv"))
ids = c('MMRF_1790','MMRF_1433','MMRF_1157')
dat = vis[SPECTRUM_SEQ%like%'MMRF' & PUBLIC_ID%in%ids]
sapply(dat, function(x) sum(is.na(x)))[1:10] #count missing

missing<-rbind(dat[,lapply(.SD,function(x)sum(is.na(x)))],
               dat[,lapply(.SD,function(x)sum(x==""))])
m.dt<-data.table(t(setDT(missing)[,lapply(.SD,sum,na.rm=T)]))
m.dt$name<-colnames(missing)
colnames(m.dt)[1]<-"n_na_missing"
m.dt[n_na_missing>3]$name
dat = dat %>% dplyr::select(!m.dt[n_na_missing>3]$name)
tmp = data.table(t(dat))
key = key.vis[!name%in%m.dt[n_na_missing>3]$name,-c(1,4,5)]
vis.dat = cbind(key,tmp)

tme = dat[,c('PUBLIC_ID','SPECTRUM_SEQ','D_IM_COLLECTION_DAY')]
rm(ids,dat,missing,m.dt,dat,tmp,key,key.vis,vis,clin)
tme = mutate(tme,day:=if_else(PUBLIC_ID=='MMRF_1790',D_IM_COLLECTION_DAY+6,D_IM_COLLECTION_DAY+13))

#Check algorithm
unique(pc$score[sample_id=="MMRF_1157_1_BM"][,-"sample_id"] == scale(data.matrix(cbat.dt[sample_id=="MMRF_1157_1_BM"][,-"sample_id"]),center=colMeans(cbat.dt[sample_id%in%baseline,-"sample_id"]),scale=FALSE) %*% pc$pca$rotation[,1:40])

#Select samples
progression <- samples[collection_reason=="Confirm Progression"]$patient_id
ids <- samples[patient_id%in%progression]$sample_id

#Compute pc score
pcs <- cbind(cbat.dt[sample_id%in%ids][,'sample_id'],
             scale(data.matrix(cbat.dt[sample_id%in%ids][,-"sample_id"]),
                   center=colMeans(cbat.dt[sample_id%in%baseline,-"sample_id"]),
                   scale=FALSE) %*% pc$pca$rotation[,1:40])
#Check pcs match for baseline samples
unique(pc$score[sample_id%in%pcs$sample_id]==pcs[sample_id%in%baseline])

#Standardize ----
base_sd <- as.data.frame(t(as.data.frame(lapply(X = pc$score[,-"sample_id"],FUN = sd))))
pcs_sd <- as.data.frame(as.matrix(pcs[,-"sample_id"]) %*% diag(1/base_sd$V1)) #standardize by baseline PC SD
colnames(pcs_sd) = PC
unique(round(pcs_sd[,"PC1_sd"],digits = 5) == round(pcs[,"PC1"]/base_sd$V1[1],digits = 5))#check algorithm
unique(round(pcs_sd[,"PC20_sd"],digits = 5) == round(pcs[,"PC20"]/base_sd$V1[20],digits = 5))#check algorithm
pcs_sd$sample_id = pcs$sample_id
dt <- merge(samples[,c('sample_id','patient_id','collection_reason')],pcs_sd,by='sample_id')

#Select samples
ids <- c('1790','1433','1157')
dat = dt[patient_id%in%ids]
dat$SPECTRUM_SEQ = gsub("_BM","",dat$sample_id)
dt = merge(tme,dat,by="SPECTRUM_SEQ")
dt = dplyr::select(dt,"patient_id","sample_id","day",starts_with("PC"))

#Transform data matrix ----
molten <- melt(dt,id.vars = c("patient_id","sample_id","day"))
base.molten = melt(score_sd,id.vars = "sample_id")
#rm(cbat.dt,pcs,dt,pc,score_sd,ids,base_sd,pcs_sd,tmp)
molten$PC = gsub("PC","",gsub("_sd","",molten$variable))
molten$timepoint = gsub("_BM","",gsub("MMRF_\\d{4}_","",molten$sample_id))
head(molten)
setDT(molten)

#Annotate quartiels ----
molten$quartile = 0
tmp = molten[0,]
p = unique(molten$variable)
for (i in seq_along(p)){
  z = molten[variable==p[i]]
  a = molten[variable==p[i]]$value
  b = base.molten[variable==p[i]]$value
  c = quantile(b,probs=0:4/4)
  d = cut(a,c,include.lowest=TRUE,labels=FALSE)
  z$quartile = d
  tmp = rbind(z,tmp)
}
molten = tmp
rm(tmp,base_sd,base.molten,cbat.dt,dt,pc,pcs,pcs_sd,samples,score_sd,z,a,b,baseline,c,d,e,i,ids,p,PC,progression)
rm(dat,tme,vis.dat)

#Plot ----
#Ploting setup
theme_set(theme_minimal() + theme(legend.position = "top"))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#line and dot plot by collection day of PC value
ggplot(data=molten[PC%in%c(6,11,12,20,25)],
       aes(x=day,y=value,group=PC,color=PC))+
  geom_point(size=3) + geom_line(size=1,aes(linetype=PC)) + facet_grid(patient_id ~ .) +
  theme(panel.border = element_rect(color = "black", fill=NA, size=1),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 13,face = "bold"),
        strip.background = element_rect(fill="gray")) +
  scale_color_brewer(palette = "Dark2")



#line and dot plot
ggplot(data=molten[PC%in%c(12,20,25)],
       aes(x=timepoint,y=quartile,group=PC,color=PC))+
  geom_point(size=3) + geom_line(size=1,aes(linetype=PC)) + facet_grid(patient_id ~ .) +
  theme(panel.border = element_rect(color = "black", fill=NA, size=1),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 13,face = "bold"),
        strip.background = element_rect(fill="gray")) +
  scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00"))

#bar chart for each timepoint
dat = molten[patient_id%in%ids[3]]
smp = unique(dat$sample_id)
ggplot(data=dat[sample_id%in%smp],aes(x=variable,y=value))+
  geom_bar(stat="identity",aes(fill=timepoint)) + facet_grid(timepoint ~ .) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
  
  
  theme(panel.border = element_rect(color = "black", fill=NA, size=1),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 13,face = "bold"),
        strip.background = element_rect(fill="gray"))
  #scale_color_manual(values=c("#3B9AB2", "#EBCC2A", "#F21A00"))

ggplot(data=molten[patient_id%in%ids & PC%in%c(6,11,12,20,25)],aes(x=PC, y=quartile, group=timepoint)) +
  geom_bar(stat="identity",position="dodge",aes(fill=timepoint)) +
  facet_grid(patient_id ~ .)


dat = molten[patient_id%in%ids]
smp = unique(dat$sample_id)
for (i in seq_along(smp)){
  p <- ggplot(data=molten[sample_id==smp[i]],
              aes(x=variable, y=value, group=sample_id)) +
    geom_bar(stat="identity",color="black") +
    #ylim(min(molten[patient_id%in%ids]$value),max(molten[patient_id%in%ids]$value))+
    annotate("text",x=20,y=max(molten[patient_id%in%ids]$value),label=ids[i],size=10)+
    scale_fill_discrete(name = "Sample", labels = c("1","2","3","4","5"))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())
  print(p)
  #  ggsave(p,file=paste("plots/multi_",ids[i],"_scaled.pdf",sep=''),scale=2)
}



#old plots
ggplot(data=molten[patient_id%in%ids & PC%in%c(6,11,12,20,25)],
       aes(x=timepoint,y=quartile,group=PC,color=PC))+
  geom_point() + geom_line() + facet_grid(variable ~ .)


ggplot(data=molten[patient_id%in%ids & PC%in%c(12,20,25)],
       aes(x=timepoint,y=value,group=PC,color=PC))+
  geom_point() + geom_line() + facet_grid(patient_id ~ .)







ggplot(data=molten[patient_id%in%ids],aes(x=PC, y=quartile, group=timepoint)) +
  geom_point() + theme(panel.grid.minor = element_blank()) +
  facet_grid(patient_id ~ .)

ggplot(data=molten[patient_id%in%ids & PC%in%c(6,11,12,20,25)],aes(x=PC, y=quartile, group=timepoint)) +
  geom_bar(stat="identity",position="dodge",aes(fill=timepoint)) +
  facet_grid(patient_id ~ .)



ggplot(data=molten[patient_id==1157 & PC==1],aes(x=timepoint,y=value,group=patient_id))+
  geom_point() + geom_line()
ggplot(data=molten[patient_id%in%ids & PC%in%c(1:5)],aes(x=timepoint,y=value,group=patient_id,color=patient_id))+
  geom_point() + geom_line() + facet_grid(PC ~ .)
ggplot(data=molten[patient_id%in%ids & PC%in%c(1,12,18,22,9)],aes(x=timepoint,y=value,group=patient_id,color=patient_id))+
  geom_point() + geom_line() + facet_grid(variable ~ .)

ggplot(data=molten[patient_id%in%ids],aes(x=PC, y=value, group=timepoint)) +
  geom_bar(stat="identity",position="dodge",aes(fill=timepoint)) +
  facet_grid(patient_id ~ .)

ggplot(data=molten[patient_id%in%ids],aes(x=PC, y=quartile, group=timepoint)) +
  geom_bar(stat="identity",position="dodge",aes(fill=timepoint)) +
  facet_grid(patient_id ~ .)

ggplot(data=molten[patient_id%in%ids],aes(x=PC, y=quartile, group=timepoint)) +
  geom_point() + theme(panel.grid.minor = element_blank()) +
  facet_grid(patient_id ~ .)

ggplot(data=molten[patient_id%in%ids & PC%in%c(6,11,12,20,25)],aes(x=PC, y=quartile, group=timepoint)) +
  geom_bar(stat="identity",position="dodge",aes(fill=timepoint)) +
  facet_grid(patient_id ~ .)




ggplot(data=molten[patient_id%in%ids & PC%in%c(1,12,18,22,9)],aes(x=PC, y=value, group=timepoint)) +
  geom_bar(stat="identity",position="dodge",aes(fill=timepoint)) +
  facet_grid(patient_id ~ .)


for (i in seq_along(ids)){
  p <- ggplot(data=molten[patient_id==ids[i]],
              aes(x=variable, y=value, group=sample_id)) +
    geom_bar(stat="identity",position="dodge",aes(fill=sample_id)) +
    ylim(min(molten[patient_id%in%ids]$value),max(molten[patient_id%in%ids]$value))+
    annotate("text",x=20,y=max(molten[patient_id%in%ids]$value),label=ids[i],size=10)+
    scale_fill_discrete(name = "Sample", labels = c("1","2","3","4","5"))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())
  print(p)
#  ggsave(p,file=paste("plots/multi_",ids[i],"_scaled.pdf",sep=''),scale=2)
}
