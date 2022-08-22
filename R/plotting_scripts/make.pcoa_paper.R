library(ggplot2)
library(vegan)
library(reshape2)
library(lmic)


key.s<-df
dat.s<-subset.dat(key.s,dat,kind) #order dat according to key
#### Make or define distance matrix####
#calculate a distnace matrix from a taxa table (relative or abosulte abundances) using bray-curtis
#or jaccard index. Result will be stored in `dista`. 
#if dat,s is a unifrac matrix (which is already a distance object), it will be directly stored to `dista`

dis="bray" # set required index - bray or jaccard:
if (kind!="unifrac"){
  if (dis=="bray") {X=F} else {X=T}
  dista=vegdist(dat.s,method=dis,binary=X)
} else {dista=dat.s}
identical(rownames(as.matrix(dista)),key.s$ID) # verifies it is still ordered same as in key

####plot PCOA on subset####
co=cmdscale(dista,k=5,eig=T) #calculates coordinates to display k dimensions on pCOA
eigs=co$eig
eigs2=eigenvals(co)
exp=round(100*(abs(eigs))/(sum(abs(eigs))),1) #gives same % explained as "PAST" does
dfw=as.data.frame(co[[1]]) #turns coordinates to format comfortable for downstream
rownames(dfw)=gsub("X","",rownames(dfw)) #just to repair sample names if tehy begain with a digit


forplot<-merge(dfw, key.s,by.x="row.names",by.y="ID") #merge corrdinates with the metadat info we have in key
forplot<-forplot[match(key.s$ID,forplot$Row.names),] #set order to be the same as in key (merge can change row order)
print(all.equal(as.character(forplot$Row.names),as.character(key.s$ID))) #verify its same order!

var1<-as.factor(key.s$Sweets)#set the variable you want to do ANOSIM on

anosim(dista,as.factor(var1),999)->a1
p1=a1$signif
r1=round(a1$statistic,2)
print(paste("anosim", "p=", p1,"| R=", r1))
i=1
j=2
name1="Sweets"

ggplot(forplot, aes(V1,V2,color=Daily_score))+
  geom_point(size=3)+
  scale_color_gradient(low="lightpink",high="blue4")+
 # geom_text(aes(label=Row.names),vjust=-0.9)+
  theme(axis.title = element_text(size=13),
        legend.title = element_text(size=12))+
#  # geom_point(aes(size=varT))+
  labs(x=paste0("Coordinate ",i," ",exp[i],"%"),y=paste0("Coordinate ",j," ",exp[j],"%"),
       color="MED score")

ggsave(paste0(pathr,"Figure3A.pdf"))
 


