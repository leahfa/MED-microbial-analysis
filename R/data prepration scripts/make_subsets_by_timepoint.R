############################################################################################################
#Run AFTER running data preporation and variable definition script 
# (make sure you sourced functions_for_downstream_analysis and installed lmic package)
############################################################################################################
library(lmic)



####make subsets for each timepoint:####
df.tp0<-df[which(df$Time=="Vis0"),]
dall.tp0<-merge(df.tp0,dat.short,by.x="ID",by.y="row.names") 
#taxa=0 not removed here (originating from dat.s)
dat0<-subset.dat(df.tp0,dat)
dat0.short<-subset.dat(df.tp0,dat.short)
abdat0<-subset.dat(df.tp0,abdat.short)
dall.tp0<-merge(df.tp0,dat.short,by.x="ID",by.y="row.names") 
abdall.tp0<-merge(df.tp0,abdat.short,by.x="ID",by.y="row.names")  

df.tp1<-df[which(df$Time=="Vis1"),]
dall.tp1<-merge(df.tp1,dat.short,by.x="ID",by.y="row.names")
dat1<-subset.dat(df.tp1,dat)
dat1.short<-subset.dat(df.tp1,dat.short)
abdat1<-subset.dat(df.tp1 ,abdat.short) 
abdall.tp1<-merge(df.tp1,abdat.short,by.x="ID",by.y="row.names") 

df.tp6<-df[which(df$Time=="Vis6"),]
dall.tp6<-merge(df.tp6,dat.short,by.x="ID",by.y="row.names")
dat6<-subset.dat(df.tp6,dat)





