############################################################################################################
# Source functions_for_downstream_analysis and install lmic package before starting!
############################################################################################################

library(lmic)


#### set paths for input metadata ####
path=paste0("C:/Users/",user,"/Dropbox/IBD/Healthy/")
user="user"
#or
path=paste0("/home/leah/Dropbox/IBD/Healthy/")

#### read in data ####
#keyh is a wide-format key, with all subject information (one line per patient)
keyh<-read.csv(paste0(path,"Pilot_table_for_analyis_110120_TF_noVitDsupp.csv"),stringsAsFactors = FALSE)
keyh[ ,-c(1:8)]<-apply(keyh[ ,-c(1:8)],2,as.numeric)
colnames(keyh)
cols.num <- c(9:ncol(keyh))
keyh[cols.num] <- sapply(keyh[cols.num],as.numeric)
sapply(keyh, class)
keyh$GI.symp<-"same"
keyh$GI.symp[which(keyh$GI.symptoms.improvement==1)]<-"better"
keyh$GI.symp[which(keyh$GI.symptoms.worsening==1)]<-"worse"
keyh$GI.symp[11]<-"same" # manual fix - has 1 for both cases

#df is a long-format key (one line per sample), lacking some subject information
df<-read.csv(paste0(path,"key_rearranged.csv"),stringsAsFactors = FALSE)
colnames(df)

#### match metadata (df) to microbial data (dat). Same samples, same order!####

keep<-intersect(rownames(dat),df$ID)
df<-df[which(df$ID %in% keep),]
dat<-subset.dat(df,dat)

#### add data from keyh to df ####
df$age<-keyh$age[match(df$PatID.full,keyh$participant)]
df$gender<-as.factor(keyh$gender[match(df$PatID.full,keyh$participant)])
df$stress<-keyh$Stress[match(df$PatID.full,keyh$participant)]
df$GI.symp<-as.factor(keyh$GI.symp[match(df$PatID.full,keyh$participant)])
df$General.feeling<-as.factor(keyh$General.Feeling[match(df$PatID.full,keyh$participant)])
df$sleep.hours<-keyh$Sleep.Hours[match(df$PatID.full,keyh$participant)]
df$steps<-keyh$Steps[match(df$PatID.full,keyh$participant)]
df$stool<-as.factor(keyh$Stool[match(df$PatID.full,keyh$participant)])
df$pulse<-keyh$pulse[match(df$PatID.full,keyh$participant)]
df$height<-keyh$height[match(df$PatID.full,keyh$participant)]
df$BMI<-df$weight/(df$height^2) 

#### define clinical and dietary variables: ####
diet.cols<-colnames(df)[18:33]
diet.cols
#remove from diet.cols variables which are close to 0):
rm=c("Soft.Drinks")
table(df$Butter)
df$Diet.Products<-as.factor(df$Diet.Products)
diet.cols<-diet.cols[which(! diet.cols %in% rm)]
diet.cols
clin.cols<-colnames(df)[3:17]
clin.cols
extra.cols<-c("exercise","age", "gender" ,"Bacterial.density","Fungal.density","stress","General.feeling",
              "GI.symp","sleep.hours","stool","steps","BMI","pulse")
extra.cols.num<-c("exercise","age", "Bacterial.density","Fungal.density","sleep.hours",
                  "steps","BMI","stress","pulse")





#### merge dat and df for easy plotting #####
# First make a version of dat with column names shortened to genus level:
dat.short<-dat
rowSums(dat)
colnames(dat.short)<-strip.names.silva.dada(colnames(dat),5)
colnames(dall)<-gsub("-","_",colnames(dall))
dall<-merge(df,dat.short,by.x="ID",by.y="row.names")

#### make a matrix of microbiome ABSOLUTE abundances ####
abdat<-make.abdat(dat,df,samp.id = "ID",dens.id = "Bacterial.density")
#### merge abdat (similar to dat, but absolute abundances) and df for easy plotting ####
temp<-subset.dat(df, abdat)
abdat.short<-temp
colnames(abdat.short)<-strip.names.silva.dada(colnames(temp),5)
abdall<-merge(df,abdat.short,by.x="ID",by.y="row.names")


#### some plotting examples ####
ggplot(dall,aes(x=Daily_score,y=Faecalibacterium))+
  #geom_line(linetype="dashed")+
  geom_point(aes(color=Time),size=3)+
  theme(axis.title.y = element_text(size=12,face="bold.italic"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.text=element_text(size=10))

####some simple correlation examples ####
cor.test(abdall$Lachnospira,abdall$Daily_score,method="spearman")
cor.test(abdall$Faecalibacterium,abdall$Daily_score,method="spearman")



