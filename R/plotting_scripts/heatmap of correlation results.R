##########################################################################################################################
library(pheatmap)
library(RColorBrewer)
extra.cols.num.f<-c("age","BMI","steps","exercise","stress","pulse","sleep.hours")
library(VennDiagram)

#########################################################################################################################
#First prepare a dataframe of correlation results (`res4`) made with function cor.2dim (from my 'lmic' package)
##########################################################################################################################

library(pheatmap)

#### example for making a correlation object ####
identical(rownames(dat.short),df$ID) #verify df and dat are perfectly matched
res4<-cor.2dim(x.cors=df[ ,diet.cols],y.cors=dat.short,dig=4, cut.level="no")

#### make a vector of  taxa which have at least one correlation with q<0.2: ####
tapply(res4$q,res4$Var1,min)->z
sort(z)
keep<-names(z[which(z<0.2)])

#### make a vector of  taxa which have at least one correlation with R > t: ####
t=0.35
tapply(res4$Rho,res4$Var1,max)->z
sort(z)
keep2<-names(z[which(z>t)])
keep2

#### make a vector of  taxa which have at least one correlation with R < t-: ####
tapply(res4$Rho,res4$Var1,min)->z
sort(z)
keep3<-names(z[which(z<(-t))])
keep3

# keepall is the final vector of which taxa to show in heatmap. Here I chose it by which taxa have q<0.2,and abs(R)>0.35, 
#But can modified as needed
keepall=intersect(union(keep3,keep2),keep)
keepall

#subset correlation table to contain only taxa which passed teh thresholds:
all.cors<-res4[which(res4$Var1 %in% keepall),]
dim(all.cors)
####Turn tabular output to a correlation matrix, which an be used for clustering,heatmap, etc.:####

matR=matrix(NA,nrow=length(unique(all.cors$Var1  )),ncol=length(unique(all.cors$Var2 )))
colnames(matR)<-unique(all.cors$Var2)
rownames(matR)<-unique(all.cors$Var1)

#Make a matrix of Rs:
for ( i in 1:length(unique(all.cors$Var1))){
  name1=unique(all.cors$Var1)[i]
  for (j in 1:length(unique(all.cors$Var2))) {
    name2<-unique(all.cors$Var2)[j]
  ind=which(all.cors$Var1==name1 & all.cors$Var2==name2)
  matR[ i,j]<-all.cors$Rho[ind]
  }
  }

  

#And a matrix of q's:

matQ=matrix(NA,nrow=length(unique(all.cors$Var1  )),ncol=length(unique(all.cors$Var2 )))
colnames(matQ)<-unique(all.cors$Var2)
rownames(matQ)<-unique(all.cors$Var1)
#Make a matQix of Rs:
for ( i in 1:length(unique(all.cors$Var1))){
  name1=unique(all.cors$Var1)[i]
  for (j in 1:length(unique(all.cors$Var2))) {
    name2<-unique(all.cors$Var2)[j]
    ind=which(all.cors$Var1==name1 & all.cors$Var2==name2)
    matQ[ i,j]<-round(all.cors$q[ind],2)
  }
}

#### for easy visibility on heatmap, replace q qith asterix (if <0.2) or  blank ####

matQ[which(matQ<=0.2)]<-"*"
matQ[which(matQ>0.2)]<-""
dim(matQ)

#### make a vactor if italized taxa names, so will appear on plot in italics ####
newnames <- lapply(
  rownames(matR),
  function(x) bquote(italic(.(x))))


pheatmap(matR,display_numbers = matQ,fontsize_number = 12,number_color = "black",
         fontsize_col= 10,fontsize_row=9,
         #breaks=seq(-0.8, 0.8,0.001),
         color=colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(1611),
        # cellwidth=25,cellheight=10,
        clustering_callback = callback,
        legend = TRUE,
       labels_row = as.expression(newnames),
         cluster_cols = TRUE)


#### to save to disk ####
dev.copy(pdf,paste0(pathr,"Figure_3C.pdf"),width=4.5,height=4)
dev.off()


####if need a matrix of p values: #####


matP=matrix(NA,nrow=length(unique(all.cors$Var1  )),ncol=length(unique(all.cors$Var2 )))
colnames(matP)<-unique(all.cors$Var2)
rownames(matP)<-unique(all.cors$Var1)
for ( i in 1:length(unique(all.cors$Var1))){
  name1=unique(all.cors$Var1)[i]
  for (j in 1:length(unique(all.cors$Var2))) {
    name2<-unique(all.cors$Var2)[j]
    ind=which(all.cors$Var1==name1 & all.cors$Var2==name2)
    matP[ i,j]<-round(all.cors$p[ind],2)
  }
}
matP[which(matP>0.05)]<-""
