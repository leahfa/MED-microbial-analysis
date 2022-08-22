##########################################################################################################################
 #Extracting  PCoA coordinates and correlating to data of interest
##########################################################################################################################
library(vegan)
library(lmic)

#assign to keyx the metadata key of samples to work on. 
#Example1:
keyx<-df  #for all the data  (all timepoints)
#Example2:
keyx<-df.tp1 #for timepoint1 only (already prepared in script make_subsets_by_timepoints). Can also be recreated:
#Example 3:
keyx<-df[which(df$Time=="Vis1"),]
#Example 2 and 3 produce the same keyx!

#match the microbiome data to fit keyx exactly (same samples in same order):
datx<-subset.dat(keyx,dat,kind)
#remove features which are 0 across all samples:
datx<-filtCols(datx,0)
#set number of coordinates to extract
k=3 
#prepare the PCoA object:
pc<-cmdscale(vegdist(datx,"bray"),k)
colnames(pc)<-paste0(rep("PCoA",),seq(1:k))
#Verify is properly matched to metadata:
identical(rownames(pc),keyx$ID)

#Do simple correlation:
cor.test(pc[ ,2],keyx$Daily_score,method="spearman")
cor.test(pc[ ,1],keyx$WBC,method="spearman")

#Simple plot
plot(pc[ ,2],keyx$Daily_score,pch=16)
plot(pc[ ,1],keyx$WBC,pch=16)

#correlate all coordinates against metadata variables of interest:
res4<-cor.2dim( pc, keyx[ ,c(diet.cols,clin.cols,"BMI","age")],4,"no")





