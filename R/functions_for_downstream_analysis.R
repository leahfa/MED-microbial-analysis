#### transform RA to absolute abundances####
#total bacterial density per sample is in column "dens.id" of teh metadata key
make.abdat<-function(data,key,samp.id,dens.id){
  abdat<-matrix(NA,ncol=ncol(data),nrow=nrow(data))
  colnames(abdat)<-colnames(data)
  rownames(abdat)<-rownames(data)
  
  for (i in 1:nrow(data)){
    z<-which(colnames(key)==dens.id)
    z1<-which(colnames(key)==samp.id)
    dens<-key[which(key[ ,z1]==rownames(data)[i]),z]
    abdat[i, ]<-dens*data[i,]
  }
  return(abdat)
  print(rowSums(abdat)==key[ ,z])
}

#### calculate 'deltas' between two timeoints ####
#mat2 and mat1: subsets of 'dat' for each timepoint.  should be non-filtered,  but with same 
#patients in same order

make.delta.bact<-function(mat2,mat1,df){
  keep<-union(colnames(filtCols(mat2,0)),colnames(filtCols(mat1,0)))
  h1<-df$PatID.full[match(rownames(mat1),df$ID)]
  h2<-df$PatID.full[match(rownames(mat2),df$ID)]
  print(identical(h1,h2))
  delta.bact<-as.data.frame(matrix(NA,ncol=length(keep),nrow=nrow(mat1)))
  for ( i in 1:length(keep)){
    delta.bact[ ,i]<-mat2[ ,which(colnames(mat2)==keep[i])]-mat1[ ,which(colnames(mat1)==keep[i])]
  }
  colnames(delta.bact)<-keep
  rownames(delta.bact)<-h1
  return(delta.bact)
}

# Delta normlaized by RA at 1st timepoint:
make.delta.bact.norm<-function(mat2,mat1,df){
  #to avid div/0, run only ontaxa alreadypresent in Vis0
  keep<-intersect(colnames(filtCols(mat2,0)),colnames(filtCols(mat1,0)))
  h1<-df$PatID.full[match(rownames(mat1),df$ID)]
  h2<-df$PatID.full[match(rownames(mat2),df$ID)]
  print(identical(h1,h2))
  delta.bact<-as.data.frame(matrix(NA,ncol=length(keep),nrow=nrow(mat1)))
  for ( i in 1:length(keep)){
    delta.bact[ ,i]<-(mat2[ ,which(colnames(mat2)==keep[i])]-
                        mat1[ ,which(colnames(mat1)==keep[i])])/mat1[ ,which(colnames(mat1)==keep[i])]
  }
  colnames(delta.bact)<-keep
  rownames(delta.bact)<-h1
  is.na(delta.bact) <- sapply(delta.bact, is.infinite)
  return(delta.bact)
}

#### Same for fungal: ####

make.delta.fung<-function(mat2,mat1,df){
  keep<-union(colnames(filtCols(mat2,1)),colnames(filtCols(mat1,2)))
  h1<-df$PatID.full[match(rownames(mat1),df$fecal.sample)]
  h2<-df$PatID.full[match(rownames(mat2),df$fecal.sample)]
  print(identical(h1,h2))
  delta.fung<-as.data.frame(matrix(NA,ncol=length(keep),nrow=nrow(mat1)))
  for ( i in 1:length(keep)){
    delta.fung[ ,i]<-mat2[ ,which(colnames(mat2)==keep[i])]-mat1[ ,which(colnames(mat1)==keep[i])]
  }
  colnames(delta.fung)<-strip.names.silva.fungal.fixed(keep,7)
  rownames(delta.fung)<-h1
  return(delta.fung)
}


####remove features taht didnt change in at least n people by at least t : ####
clean.deltas<-function(delta.bact,n,t){
  keep=c()
  counter=1
  for( i  in 1:ncol(delta.bact)){
    x<-delta.bact[ ,i]
    y<-abs(x)
    numy<-length(y[which(y>0)])
    maxy<-max(y,na.rm=TRUE)
    if (numy>n & maxy>t) {
      keep[counter]<-colnames(delta.bact)[i]
      counter<-counter+1
    }
    
  }
  return(delta.bact[ ,which(colnames(delta.bact) %in% keep)])
}

###run KW test in loop across all (filtered) taxa against all categorical metadata ####
kruskal.2dim<-function(datobj,facobj,dig) {
  reslist<-list()
  for ( i in 1:ncol(datobj)) {
    resdf<-as.data.frame(matrix(NA,ncol=3,nrow=ncol(facobj)))
    colnames(resdf)<-c("Var","Factor","p")
    resdf$Var<-rep(colnames(datobj)[i],nrow(resdf))
    
    y=datobj[ ,i]
    for ( j in 1:ncol(facobj)){
      resdf[j,2]<-colnames(facobj)[j]
      resdf[j,3]<-round(kruskal.test(y~as.factor(facobj[,j]))$p.value,dig)
    }
    reslist[[i]]<-resdf
  }
  resfin<-do.call("rbind",reslist)
  resfin$q<-p.adjust(resfin$p,"fdr")
  return(resfin)
}