#Install relevant packages
install.packages("dplyr")
install.packages("data.table")
install.packages("ggplot2")

#Load package
library("dplyr")
library("data.table")
library("ggplot2")

#Read data file of pc scores 
setwd("/Users/rosal/OneDrive - University of Utah/2020/analyze/data/transcriptome-dimensions/")
load(file='rdata/process_20200428.RData') #pcs and clinical data in dt
base.pc.clin.dt <- merge(samples,score_sd,by='sample_id')
data <- base.pc.clin.dt %>% dplyr::select("sample_id","patient_id",ends_with("_sd")) %>% na.omit()
head(data)
rm(base.pc.clin.dt,cbat.dt,pc,score_sd,samples)

#Tidy data
molten <- melt(data,id.vars = c("sample_id","patient_id"))
molten$sign <- ifelse(molten$value >= 0, "positive", "negative")

min(molten$value)
max(molten$value)

#Select extreem examples
data$sum <- data %>% dplyr::select(ends_with("_sd")) %>% rowSums() #sum standardized pc values
setorder(data,sum)
ids <- c(head(data)[]$patient_id,tail(data)[]$patient_id)
ex <- molten[patient_id%in%ids]

#Plot
for (i in seq_along(ids)){
  p <- ggplot(data=ex[variable%like%"_sd" & patient_id==ids[i]],
              aes(x=variable, y=value, fill=sign)) +
    geom_bar(stat="identity") +
    ylim(min(ex$value),max(ex$value))+
    #annotate("text",x=2,y=min(ex$value),label=ids[i],size=5)+
    theme(legend.position = "none")+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())
  print(p)
  ggsave(p,file=paste("plots/",ids[i],"_scaled.pdf",sep=''),scale=2)
}
