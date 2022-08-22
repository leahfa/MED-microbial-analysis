#########################################################################################################################
#Script data preparation and variable definition must be run first!
##########################################################################################################################



vars=c(diet.cols,clin.cols,extra.cols)
identical(df$ID,rownames(dat))
dim(df)
# function for looping through adonis models for logintudinal data:
make.adonis.sum.tab<-function(df,dat,vars) {
adonis.sum.tab<-as.data.frame(matrix(NA, nrow=length(vars),ncol=3))
colnames(adonis.sum.tab)<-c("Variable","R2","p")
for (i in 1:length(vars)){
  print(i)
  ind<-which(colnames(df)==vars[i])
  df$Var<-df[ ,ind]
  if (length(table(is.na(df$Var)))==1) {
x<-adonis(dat ~ Var   , data=df,strata=df$PatID.full,method="bray")
adonis.sum.tab[i,1]<-vars[i]
adonis.sum.tab[i,2]<-x[[1]]$R2[1]
adonis.sum.tab[i,3]<-x[[1]]$`Pr(>F)`[1]
  } else{
    df.s<-df[which(is.na(df$Var)==FALSE),]
    dat.s<-filtCols(subset.dat(df.s,dat,kind),0)
    x<-adonis(dat.s~ Var , data=df.s,strata=df.s$PatID.full,method="bray")
    adonis.sum.tab[i,1]<-vars[i]
    adonis.sum.tab[i,2]<-x[[1]]$R2[1]
    adonis.sum.tab[i,3]<-x[[1]]$`Pr(>F)`[1]
 
}
}
adonis.sum.tab$q<-p.adjust(adonis.sum.tab$p,"fdr")

return(adonis.sum.tab)
}
# function for looping through adonis models for cross-sectional data:
make.adonis.sum.tab2<-function(df,dat,vars) {
  adonis.sum.tab<-as.data.frame(matrix(NA, nrow=length(vars),ncol=3))
  colnames(adonis.sum.tab)<-c("Variable","R2","p")
  for (i in 1:length(vars)){
    print(i)
    ind<-which(colnames(df)==vars[i])
    df$Var<-df[ ,ind]
    if (length(table(is.na(df$Var)))==1) {
      x<-adonis(dat ~ Var   , data=df,method="bray")
      adonis.sum.tab[i,1]<-vars[i]
      adonis.sum.tab[i,2]<-x[[1]]$R2[1]
      adonis.sum.tab[i,3]<-x[[1]]$`Pr(>F)`[1]
    } else{
      df.s<-df[which(is.na(df$Var)==FALSE),]
      dat.s<-filtCols(subset.dat(df.s,dat,kind),0)
      x<-adonis(dat.s~ Var , data=df.s,method="bray")
      adonis.sum.tab[i,1]<-vars[i]
      adonis.sum.tab[i,2]<-x[[1]]$R2[1]
      adonis.sum.tab[i,3]<-x[[1]]$`Pr(>F)`[1]
      
    }
  }
  adonis.sum.tab$q<-p.adjust(adonis.sum.tab$p,"fdr")
  
  return(adonis.sum.tab)
}




adonis.bact.clin<-make.adonis.sum.tab(df,dat,clin.cols)
adonis.bact.nut<-make.adonis.sum.tab(df,dat,diet.cols)
adonis.bact.extras<-make.adonis.sum.tab(df,dat,extra.cols)



# can also be done on subsets (if created already, in script make_subsets_by_timepoint)

adonis.all.tp1<-make.adonis.sum.tab2(df.tp1,dat1,vars)
adonis.all.tp0<-make.adonis.sum.tab2(df.tp0,dat0,vars)
