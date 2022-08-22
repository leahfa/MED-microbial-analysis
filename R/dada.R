
library(dada2)

path1="/scratch200/leah/MDT1/Order_431844/files/" #path to files cleaned of primers, reverse primers, and LSU/SSU sequences
path="/scratch200/leah/MDT1/Order_431844/" #folder am working in
taxa.ref<-"/urigo/leah/resources/silva_nr_v132_train_set.fa.gz"
species.ref<-"/urigo/leah/resources/silva_species_assignment_v132.fa.gz"
dir.create(path,"files_filt_by_dada/")
names_F=list.files(path1,".R1")
names_R=list.files(path1,".R2")
fnFs=paste0(path1,names_F)
fnRs=paste0(path1,names_R)

temp<-sapply(strsplit(basename(names_F), "_S"), `[`, 1)
#sample.names<-sapply(strsplit(temp,"-"),"[",2)
sample.names<-temp
filtFs <- paste0(path,"files_filt_by_dada/",sample.names,"_R1.fastq.gz") #make vector of names for filrtered files. add ".gz" if need to compress. 
filtRs <- paste0(path,"files_filt_by_dada/",sample.names,"_R2.fastq.gz") #make vector of names for filrtered files. add ".gz" if need to compress.



out <- filterAndTrim(fnFs, filtFs, fnRs,filtRs, maxN = 0, maxEE=c(2,2),truncQ = 2,truncLen=c(210,210),rm.phix = TRUE, compress =TRUE, multithread = TRUE)
out=as.data.frame(out) #format out for easy viewing
out$per.survived<-out$reads.out/out$reads.in
out$SampID=gsub(".fastq.gz","",rownames(out))
out
write.table(out,paste0(path,"summary_dadafilt.txt"), row.names=F,sep="\t",quote=F)


keep<-which(names_F %in% rownames(out))

names_F<-names_F[keep]
names_R<-names_R[keep]
filtFs<-filtFs[keep]
filtRs<-filtRs[keep]




errF=learnErrors(filtFs,multithread=TRUE)
#png(paste0(path,"plotErrorsF.png"),width=700,height=600)
#plotErrors(errF, nominalQ=TRUE)
#dev.off()

errR=learnErrors(filtRs,multithread=TRUE)
#png(paste0(path,"plotErrorsR.png"),width=700,height=600)
#plotErrors(errR, nominalQ=TRUE)
#dev.off()

print(dada2:::checkConvergence(errF))
print(dada2:::checkConvergence(errR))


derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names

names(derepFs) <- names_F
names(derepRs) <- names_R


dadaFs=dada(derepFs,err=errF, multithread = TRUE)
dadaRs=dada(derepRs,err=errR, multithread = TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
print(dim(seqtab))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
print(sort(rowSums(seqtab.nochim)/rowSums(seqtab)))

saveRDS(seqtab.nochim,paste0(path,"st.alltax.rds"))

#### assign taxonomy ####
st.filt<-seqtab.nochim

taxa <- assignTaxonomy(st.filt,taxa.ref , multithread = 6, tryRC = FALSE)
taxa<-addSpecies(taxa,species.ref,tryRC=FALSE,allowMultiple=TRUE)
#### assign asvIDs to sequence table conames ####
#for easier downstream analysis, use asvID instead of seqeunecs in sequence table. Add corresponding asvIDs to the taxa table, 
#which retains the sequence information and thus serves as a complete asvID/sequence/taxonomy lookup key
identical(rownames(taxa),colnames(st.filt))
taxa<-as.data.frame(taxa)
taxa$asvID<-paste0("asv",seq(1:nrow(taxa))) #give each raw sequence a asv identifer
colnames(st.filt)<-taxa$asvID 
taxa$seq<-rownames(taxa)
rownames(taxa)<-NULL
write.table(st.filt,paste0(path,"st.alltax.txt"),sep="\t")
write.csv(taxa, paste0(path,"taxa.csv"))

