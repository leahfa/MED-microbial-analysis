#####################################################################################################
# Run this after running the scripts for  data preparation and making subsets 
# Make sure you have package 'lmic' installed and also sourced the file functions_for_downstream_analysis
#####################################################################################################

library(lmic)
library(ggplot2) # (optional)

#### Pre-arrange metadata ####
#subset data to  2 timepoints:
dfs<-df[which(df$Time!="Vis6"),]
dim(dfs)
#arrange data  first by Visit and tehn by Subject: 
dfs<-dfs[order(df$Time,dfs$PatID.full),]
identical(dfs$PatID.full[1:20],dfs$PatID.full[21:40])

#### make a dataframe of dieatary deltas: ####
delta.nut<-as.data.frame(matrix(NA,ncol=length(diet.cols),nrow=20))
for ( i in 1:length(diet.cols) ){
 delta.nut[ ,i]<-dfs[which(dfs$Time=="Vis1"),which(colnames(dfs)==diet.cols[i])]-
   dfs[which(dfs$Time=="Vis0"),which(colnames(dfs)==diet.cols[i])]
}
colnames(delta.nut)<-diet.cols
rownames(delta.nut)<-dfs$PatID.full[1:20]

#is also possible (but not advised in my opinion) to normalize by value at 1st timepoint:
delta.nut.norm<-as.data.frame(matrix(NA,ncol=length(diet.cols),nrow=20))
for ( i in 1:length(diet.cols) ){
  delta.nut.norm[ ,i]<-(dfs[which(dfs$Time=="Vis1"),which(colnames(dfs)==diet.cols[i])]-
    dfs[which(dfs$Time=="Vis0"),which(colnames(dfs)==diet.cols[i])])/dfs[which(dfs$Time=="Vis0"),
                                                  which(colnames(dfs)==diet.cols[i])]
}
#variables which have 0 at 1st timepoint cant be normalized (unless zero-imputed first, bu then, with 
#what value to impute?), so changing them from Inf (inifinte) to NA:
is.na(delta.nut.norm) <- sapply(delta.nut.norm, is.infinite)
colnames(delta.nut.norm)<-diet.cols
rownames(delta.nut.norm)<-dfs$PatID.full[1:20]

#### Make dataframe of clinical deltas: ####
#note: 1 Subject missing bloodtests
delta.clin<-as.data.frame(matrix(NA,ncol=length(clin.cols),nrow=20))
for ( i in 1:length(clin.cols) ){
  delta.clin[ ,i]<-dfs[which(dfs$Time=="Vis1"),which(colnames(dfs)==clin.cols[i])]-dfs[which(dfs$Time=="Vis0"),
                                                                which(colnames(dfs)==clin.cols[i])]
}
colnames(delta.clin)<-clin.cols
rownames(delta.clin)<-dfs$PatID.full[1:20]


#Normalized version:
delta.clin.norm<-as.data.frame(matrix(NA,ncol=length(clin.cols),nrow=20))
for ( i in 1:length(clin.cols) ){
  delta.clin.norm[ ,i]<-dfs[which(dfs$Time=="Vis1"),which(colnames(dfs)==clin.cols[i])]-
    dfs[which(dfs$Time=="Vis0"),which(colnames(dfs)==clin.cols[i])]/dfs[which(dfs$Time=="Vis0"),which(colnames(dfs)==clin.cols[i])]
    
}
colnames(delta.clin.norm)<-clin.cols
rownames(delta.clin.norm)<-dfs$PatID.full[1:20]


#### Make dataframe of bacterial deltas ####
# object dat1.short and dat0.short mafde in script make_subscripts_by_timepoint
delta.bact<-make.delta.bact(dat1.short,dat0.short,dfs)
delta.bact.norm<-make.delta.bact.norm(dat1.short,dat0.short,dfs)
# and for  absolute deltas (deltas of absolute rather tahn relative abundances):
delta.abact<-make.delta.bact(abdat1,abdat0,dfs)
delta.abact.norm<-make.delta.bact.norm(abdat1,abdat0,dfs)
#verify all is well:
identical(rownames(delta.nut),rownames(delta.bact))
identical(rownames(delta.clin),rownames(delta.bact))
hist(delta.bact$Bacteroides,breaks=20)
hist(delta.abact$Bacteroides,breaks=20)
#### Some simple exampels of correlating  delta objects ####
cor.test(delta.bact.norm$Faecalibacterium,delta.nut.norm$Daily_score,method="spearman")
cor.test(delta.bact$Faecalibacterium,delta.nut$Daily_score,method="spearman")


#### Some examples of correlation tables for delta variables ####
rdcal.bact<-cor.1dim(clean.deltas(delta.bact,6,0.005),delta.clin$Calprotectin,4)
rdcal.abact<-cor.1dim(clean.deltas(delta.abact,6,0.0005),delta.clin$Calprotectin,4)
rdcal.nut<-cor.1dim(delta.nut,delta.clin$Calprotectin,4)
res4<-cor.2dim(delta.nut,clean.deltas(delta.bact,8,0.03),4,"no")
res4<-cor.2dim(clean.deltas(delta.abact,12,0.01),delta.nut[ ,diet.cols],4,"no")



####for plotting ####
#easiest to plot whne all variables are in same object!
delta.nut.bact<-merge(delta.nut,delta.bact,by="row.names")
delta.nut.bact.norm<-merge(delta.nut.norm,delta.bact.norm,by="row.names")
delta.nut.abact<-merge(delta.nut,delta.abact,by="row.names")
delta.nut.abact.norm<-merge(delta.nut.norm,delta.abact.norm,by="row.names")
delta.clin.bact.norm<-merge(delta.clin.norm,delta.bact.norm,by="row.names")
delta.clin.bact<-merge(delta.clin,delta.bact,by="row.names")
delta.clin.abact<-merge(delta.clin,delta.abact,by="row.names")
delta.nut.clin<-merge(delta.nut,delta.clin,by="row.names")
#merge all 3 deltas to one dataframe:
rownames(delta.nut.clin)<-delta.nut.clin$Row.names
delta.nut.clin.bact<-merge(delta.nut.clin[ ,-1],delta.bact,by="row.names")
delta.nut.clin.abact<-merge(delta.nut.clin[ ,-1],delta.abact,by="row.names")

#add in data of intesrt - for example, what was the actual Daily Score at Time 0?
delta.nut.bact$Daily_score_Vis0<-df.tp0$Daily_score[match(delta.nut.bact$Row.names,df.tp0$PatID.full)]
delta.nut.abact$Daily_score_Vis0<-df.tp0$Daily_score[match(delta.nut.abact$Row.names,df.tp0$PatID.full)]
cor.test(delta.nut.bact$Lachnospira,delta.nut.bact$Daily_score,method="spearman")
ggplot(delta.nut.abact,aes(x=Daily_score, y=Lachnospira,color=Daily_score_Vis0))+
  geom_point(size=3)+
  scale_color_gradient(low="red",high="green")+
  labs(y="Delta Lachnospira (Vis1-Vis0)",x="Delta Daily_score (Vis1-Vis0)")
  

ggplot(delta.nut.clin,aes(x=Yogurt, y=HDL,color=Daily_score))+
  geom_point(size=3)+
  labs(y="delta Cholesterol (Vis1-Vis0)" ,x= "delta Yogurt (Vis1-Vis0")
  #geom_text(aes(label=Row.names),hjust=-0.5)

ggplot(delta.nut.abact,aes(x=Sweets,y=Blautia))+
  geom_point()
rownames(delta.bact)
