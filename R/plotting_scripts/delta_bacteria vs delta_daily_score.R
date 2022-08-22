#####################################################################################################
#Make sure you have package 'lmic' loaded and also sourced the file functions_for_downstream_analysis
# Script calculate_deltas_Vis1_Vis0 must be done first!
#####################################################################################################
library(ggplot2)
library(reshape2)
library(lmic)

dim(delta.nut.abact) #verify you have this object ready
#make correlation table:
res<-cor.1dim(clean.deltas(delta.abact,10,0.005),delta.nut$Daily_score,4)
#Tag taxa which pass significance thresholds:
tax<-res$Taxon[which(res$p<0.05 & res$FDR<0.2 & res$Rho>0.4)]
inds<-which(colnames(delta.nut.abact) %in% tax)
# Subset correlation results to significant intercations only:
temp<-delta.nut.abact[, c(1,grep("Daily_score",colnames(delta.nut.abact)), inds)]
# Melt dataframe for easy plotting:
dm<-melt(temp,id.vars=c("Daily_score","Row.names") )


ggplot(dm, aes(x=Daily_score,y=value ))+
  geom_point()+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE)+
  labs(x="delta MED score (Tend-T0)",
       y=paste0("delta Bacterial Absolute Abundances \n  (\u00b5g/mg)"))+
  theme(axis.title.y = element_text(face="italic"),
        strip.text = element_text(face="italic",size=8))+
  #geom_text(aes(label=Row.names),size=2)+
  facet_wrap(~variable,ncol=4,scales="free_y")

#make and write out correlation table:
identical(rownames(delta.nut),rownames(delta.abact))
res<-cor.1dim(clean.deltas(delta.abact,10,0),delta.nut$Daily_score,4)
write.csv(res, paste0(path,"delta_abact_daily_score.csv"))


