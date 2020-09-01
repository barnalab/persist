library(tidyverse)
library(reshape2)
library(data.table)
library(viridis)
library(scales)
library(limma)
library(edgeR)
library(zoo)
library(hues)
library(ggrepel)
library(deSolve)
library(Biostrings)

setwd("~/vac/233x/analysis/")
options(stringsAsFactors=F)

#Construct annotations
constructs=data.frame(fread("./Collated.tsv",sep="\t",header=T))
constructs.use=subset(constructs,Barcode.Amplicon.Sequencing.2%in%c("TRUE","spike-in"))
rownames(constructs.use)=constructs.use$SequenceID
constructs.use.spikein=subset(constructs.use,Barcode.Amplicon.Sequencing.2=="spike-in")
constructs.use.others=subset(constructs.use,Barcode.Amplicon.Sequencing.2!="spike-in")

#Rename some annoying names
constructs.use[,1]=make.names(constructs.use[,1],unique=T)
constructs.use[grep("CoV(.*)Nluc_hHBB",constructs.use[,1]),3]=c("vary5UTR_Nluc_HBB_CoV2") #CoV2 5'UTR mutants
constructs.use$Experiment[constructs.use$Experiment=="HBB_Nluc_varyCDS"]="HBB_varyCDS-NLuc_HBB_F30Pepper"
constructs.use$Experiment[constructs.use$Experiment=="vary5UTR_NlucP_HBA1"]="NlucP"
constructs.use$Experiment[constructs.use$Experiment=="vary5UTR_NlucP_HBB"]="NlucP"

constructs.use$CDSlen=nchar(constructs.use$CDS.sequence)

constructs.use$Experiment=factor(constructs.use$Experiment,levels=unique(constructs.use$Experiment)[c(3,2,1,4,11,9,10,7,5,6,8)]) #Group orders for plotting

#Counts table
dat=data.frame(fread("./joined.txt",header=T))
rownames(dat)=dat[,1]
colnames(dat)=c("id",unlist(lapply(strsplit(colnames(dat)[-1],"\\."),"[[",4)))
dat=dat[-nrow(dat),]
dat.mat=dat[,-1]
rownames(dat.mat)=dat[,1]

#All raw counts
write.table(dat.mat,file="./tables/all_counts.tsv",sep="\t",row.names=T,col.names=T,quote=F)

#Sample annotations
samples=data.frame(fread("./samples.csv",header=T))
rownames(samples)=samples$index
spikeins=subset(constructs.use,Barcode.Amplicon.Sequencing.2=="spike-in")$SequenceID
usespikeins=spikeins[1:3]

#Remove fraction 1, up to fraction 11 for sum
fraction_start=2
fraction_end=16
weight_fraction_end=11

#In-cell stability, long
icstability_long_spikein=subset(samples,(Experiment=="stability")&(Spikein==1)&(Group=="In-cell")&(Type=="Long"))$index
icstability_long_spikein_labels=samples[icstability_long_spikein,"Time"]

dat.mat.icstability_long_spikein=dat.mat[,icstability_long_spikein]
colnames(dat.mat.icstability_long_spikein)=icstability_long_spikein_labels

#Flag no primer binding site constructs
long_used=constructs.use[unique(c(grep("acatttgcttctgacacaactgtgttcac",tolower(constructs.use$Full.sequence)),grep("ATGGCCGTTTACCCATACGATGTTCCTGAC",toupper(constructs.use$Full.sequence)))),]$SequenceID #remove no HBB or Nluc since they are low count PCR artifacts
dat.mat.icstability_long_spikein=dat.mat.icstability_long_spikein[long_used,]
dat.mat.icstability_long_spikein.log2=log2(dat.mat.icstability_long_spikein)

icstability_long.design=model.matrix(~0+factor(icstability_long_spikein_labels))
colnames(icstability_long.design)=unlist(lapply(strsplit(colnames(icstability_long.design),")"),"[[",2))

#Spike-in log linear fit
'
ID	Length	Conc	uM	dilution
120002B1	897	1057	3.46	1000
120010B1	897	1071	3.51	100
220023B1	882	996	3.32	10
310333T3	953	1058	3.26	1
'
spikein_concs=c(3460,351,33.2)#,3.26)
dat.mat.icstability_long_spikein.log2.scaled=dat.mat.icstability_long_spikein.log2
for (i in 1:ncol(dat.mat.icstability_long_spikein.log2))
{
  spikeinfit=lm(dat.mat.icstability_long_spikein.log2[usespikeins,c(1,i)])
  dat.mat.icstability_long_spikein.log2.scaled[,i]=predict.lm(spikeinfit,newdata=dat.mat.icstability_long_spikein.log2[,c(1,i)])
  
  #Using actual amounts
  #spikeinfit=lm(data.frame(std=log2(spikein_concs),dat=dat.mat.icstability_long_spikein.log2[usespikeins,i]))
  #dat.mat.icstability_long_spikein.log2.scaled[,i]=predict.lm(spikeinfit,newdata=data.frame(dat=dat.mat.icstability_long_spikein.log2[,c(i)]))
}

rownames(dat.mat.icstability_long_spikein.log2.scaled)=rownames(dat.mat.icstability_long_spikein.log2)
v.icstability_long=voom(2^(dat.mat.icstability_long_spikein.log2.scaled[!rownames(dat.mat.icstability_long_spikein.log2.scaled)%in%spikeins,]),icstability_long.design,lib.size=rep(1,15))

icstability_long.fit=lmFit(v.icstability_long,icstability_long.design)
icstability_long.fit.coef=icstability_long.fit$coefficients
icstability_long.fit.se=icstability_long.fit$stdev.unscaled*icstability_long.fit$sigma
icstability_long.fit.upper=icstability_long.fit.coef+icstability_long.fit.se
icstability_long.fit.lower=icstability_long.fit.coef-icstability_long.fit.se


icstability_long.fit.coef.scale=icstability_long.fit.coef-matrix(apply(icstability_long.fit.coef,1,max),nrow=nrow(icstability_long.fit.coef),ncol=ncol(icstability_long.fit.coef),byrow=F)
icstability_long.fit.upper.scale=icstability_long.fit.upper-matrix(apply(icstability_long.fit.upper,1,max),nrow=nrow(icstability_long.fit.upper),ncol=ncol(icstability_long.fit.upper),byrow=F)
icstability_long.fit.lower.scale=icstability_long.fit.lower-matrix(apply(icstability_long.fit.lower,1,max),nrow=nrow(icstability_long.fit.lower),ncol=ncol(icstability_long.fit.lower),byrow=F)

icstability_long.fit.coef.scale.melt=melt(icstability_long.fit.coef.scale)
colnames(icstability_long.fit.coef.scale.melt)=c("id","time","value")
icstability_long.fit.coef.scale.melt$upper=melt(icstability_long.fit.upper.scale)$value
icstability_long.fit.coef.scale.melt$lower=melt(icstability_long.fit.lower.scale)$value
icstability_long.fit.coef.scale.melt$time=factor(icstability_long.fit.coef.scale.melt$time)
icstability_long.fit.coef.scale.melt$totalcount=log2(rowSums(2^icstability_long.fit.coef)[icstability_long.fit.coef.scale.melt$id])

#Try to fit decay rate (of sorts) for ordering
icstability_long.fit.coef.dr=list() 
#Let's just take average since only 2 timepoints to fit
for(i in rownames(icstability_long.fit.coef))
{
  #icstability_long.fit.coef.dr[i]=lm(icstability_long.fit.coef[i,1:3]~c(1,7,12))$coefficients[2] #Actual fitting (base 2)
  #ln N/No = -decay_rate*t
  decay_rate_6=-log(2^(icstability_long.fit.coef.scale[i,2]))/6
  decay_rate_11=-log(2^(icstability_long.fit.coef.scale[i,3]))/11
  decay_rate=mean(c(decay_rate_6,decay_rate_11))
  icstability_long.fit.coef.dr[i]=decay_rate
}
icstability_long.fit.coef.dr=unlist(icstability_long.fit.coef.dr)


#Fraction of mRNA remainingrelative to time zero
write.table(2^icstability_long.fit.coef.scale,file="./tables/long_spikein_incellstability_percent.tsv",row.names=T,col.names=T,quote=F,sep="\t")

icstability_long.fit.coef.scale.melt$id=factor(icstability_long.fit.coef.scale.melt$id,levels=rownames(icstability_long.fit.coef.scale)[order(icstability_long.fit.coef.dr[rownames(icstability_long.fit.coef.scale)])])

fig.icstability.long=ggplot(subset(icstability_long.fit.coef.scale.melt),aes(x=time,y=id,fill=value,colour=value))+theme_classic()+geom_tile()+scale_fill_viridis(option="magma",oob=squish,na.value="darkgrey")+scale_colour_viridis(option="magma",oob=squish,na.value="darkgrey")+xlab("")+ylab("")+scale_y_discrete(labels=constructs.use[levels(icstability_long.fit.coef.scale.melt$id),1])#,limits=c(0,0.5))

pdf("./plots/long_spikein_incellstability.pdf",width=8,height=20)
print(fig.icstability.long)
dev.off()

fig.icstability.long.line=ggplot(subset(icstability_long.fit.coef.scale.melt),aes(x=time,y=value,group=id,colour=totalcount))+theme_classic()+geom_line(alpha=0.33)+scale_colour_viridis(name="Total abundance",option="magma")+xlab("Time")+ylab("Log abundance")

pdf("./plots/long_spikein_incellstability_line.pdf",width=6,height=6)
print(fig.icstability.long.line)
dev.off()

#Fractions
fraction_spikein=subset(samples,(Experiment=="polysomes")&(Spikein==1))$index
fraction_spikein_labels=samples[fraction_spikein,"Fraction"]

dat.mat.fraction_spikein=dat.mat[,fraction_spikein]
colnames(dat.mat.fraction_spikein)=fraction_spikein_labels
dat.mat.fraction_spikein.log2=log2(dat.mat.fraction_spikein)

#Log linear fit
dat.mat.fraction_spikein.log2.scaled=dat.mat.fraction_spikein.log2
for (i in 1:ncol(dat.mat.fraction_spikein.log2))
{
  spikeinfit=lm(dat.mat.fraction_spikein.log2[usespikeins,c(1,i)])
  dat.mat.fraction_spikein.log2.scaled[,i]=predict.lm(spikeinfit,newdata=dat.mat.fraction_spikein.log2[,c(1,i)])
}
rownames(dat.mat.fraction_spikein.log2.scaled)=rownames(dat.mat.fraction_spikein.log2)
v.fraction=voom(2^dat.mat.fraction_spikein.log2.scaled[!(rownames(dat.mat.fraction_spikein.log2.scaled)%in%spikeins),],lib.size=rep(1,16))[,fraction_start:fraction_end]
v.fraction.scale=2^v.fraction$E/rowSums(2^v.fraction$E)
rownames(v.fraction.scale)=rownames(v.fraction)
colnames(v.fraction.scale)=colnames(v.fraction)

#v.fraction.weights=c(0,0,0,1,1.75,2.75,3.5,4.5,5.75,7,8.5,10,12,14,17,20)[fraction_start:fraction_end]
v.fraction.weights=c(0,0,0,1,1.5,2.5,3.5,4.5,5.5,7,8.5,10,12,14,17,20)[fraction_start:fraction_end]
v.fraction.scale.weight=v.fraction.scale*matrix(v.fraction.weights,nrow=nrow(v.fraction.scale),ncol=ncol(v.fraction.scale),byrow=T)
rownames(v.fraction.scale.weight)=rownames(v.fraction)
colnames(v.fraction.scale.weight)=colnames(v.fraction)

#v.fraction.group.zscale
counts.scale.list=list()
for(cdslen in constructs.use[rownames(v.fraction.scale),"CDSlen"])
{
  cdslenids=subset(constructs.use[rownames(v.fraction.scale),],CDSlen==cdslen)$SequenceID
  counts=dat.mat.fraction_spikein[cdslenids,fraction_start:weight_fraction_end]
  counts.scale=t(apply(counts/matrix(colSums(counts),nrow=nrow(counts),ncol=ncol(counts),byrow=T),1,scale,center=T,scale=T))
  colnames(counts.scale)=colnames(counts)
  counts.scale.list[[as.character(cdslen)]]=counts.scale
}
counts.scale=do.call(rbind,counts.scale.list)

v.fraction.scale.melt=melt(v.fraction.scale)
v.fraction.scale.weight.melt=melt(v.fraction.scale.weight)
v.fraction.zscale.melt=melt(counts.scale)
colnames(v.fraction.scale.melt)=c("id","fraction","value")
colnames(v.fraction.scale.weight.melt)=c("id","fraction","value")
colnames(v.fraction.zscale.melt)=c("id","fraction","value")

v.fraction.scale.melt$fraction=factor(v.fraction.scale.melt$fraction)
v.fraction.scale.weight.melt$fraction=factor(v.fraction.scale.weight.melt$fraction)
v.fraction.scale.melt$totalcount=rowSums(2^v.fraction$E)[v.fraction.scale.melt$id]
v.fraction.zscale.melt$fraction=factor(v.fraction.zscale.melt$fraction)

#Order by row sum of weighted %mRNA = "mean ribosome load"; using up to fraction 12
v.fraction.scale.melt$id=factor(v.fraction.scale.melt$id,levels=rownames(v.fraction.scale)[order(constructs.use[rownames(v.fraction.scale),"Experiment"],rowSums(v.fraction.scale.weight[,1:(weight_fraction_end-fraction_start)]))])
v.fraction.scale.weight.melt$id=factor(v.fraction.scale.weight.melt$id,levels=rownames(v.fraction.scale.weight)[order(constructs.use[rownames(v.fraction.scale),"Experiment"],rowSums(v.fraction.scale.weight[,1:(weight_fraction_end-fraction_start)]))])

v.fraction.scale.melt$labels=factor(constructs.use[as.character(v.fraction.scale.melt$id),1],levels=constructs.use[levels(v.fraction.scale.melt$id),1])
v.fraction.scale.weight.melt$labels=factor(constructs.use[as.character(v.fraction.scale.weight.melt$id),1],levels=constructs.use[levels(v.fraction.scale.weight.melt$id),1])

#%mRNA
write.table(v.fraction.scale,file="./tables/fractions_percentmrna.tsv",row.names=T,col.names=T,quote=F,sep="\t")

#Weighted %mRNA
write.table(v.fraction.scale.weight,file="./tables/fractions_weightedpercentmrna.tsv",row.names=T,col.names=T,quote=F,sep="\t")

#z-scaled
write.table(counts.scale,file="./tables/fractions_zscaled.tsv",row.names=T,col.names=T,quote=F,sep="\t")

fig.fractions.line=ggplot(v.fraction.scale.melt,aes(x=fraction,y=value,group=id,colour=log2(totalcount)))+theme_classic()+geom_line(alpha=0.33)+scale_colour_viridis()+theme(legend.position="none")
pdf("./plots/spikein_fractions_line.pdf",width=6,height=6)
print(fig.fractions.line)
dev.off()

fig.fractions=ggplot(v.fraction.scale.melt,aes(x=fraction,y=labels,fill=(value),colour=(value)))+theme_classic()+geom_tile()+scale_fill_viridis(name="% mRNA",option="magma",oob=squish,na.value="darkgrey")+scale_colour_viridis(name="% mRNA",option="magma",oob=squish,na.value="darkgrey")+xlab("")+ylab("")
pdf("./plots/spikein_fractions.pdf",width=10,height=24)
print(fig.fractions)
dev.off()

fig.fractions.weight=ggplot(v.fraction.scale.weight.melt,aes(x=fraction,y=labels,fill=(value),colour=(value)))+theme_classic()+geom_tile()+scale_fill_viridis(name="Weighted % mRNA",option="magma",oob=squish,na.value="darkgrey")+scale_colour_viridis(name="Weighted % mRNA",option="magma",oob=squish,na.value="darkgrey")+xlab("")+ylab("")
pdf("./plots/spikein_fractions_weighted.pdf",width=10,height=24)
print(fig.fractions.weight)
dev.off()

#Comparison with previous 64x data
previds1=subset(constructs,Barcode.Amplicon.Sequencing.1=="TRUE"&Designer=="Barna")$SequenceID
previds2=subset(constructs,Barcode.Amplicon.Sequencing.2=="TRUE"&Designer=="Barna")$SequenceID
previds2=previds2[!previds2%in%spikeins]
previds.names=constructs.use[previds1,1]
names(previds.names)=previds1
previds.names[is.na(previds.names)]=subset(constructs,SequenceID%in%previds1[is.na(previds.names)])[,1]

#Previous 64x data
prevdat=read.csv("./Relative_abundances.csv")
colnames(prevdat)[1]="id"
prevdat[prevdat$id%in%prevdat$id[duplicated(prevdat$id)],][1:3,-1]=prevdat[prevdat$id%in%prevdat$id[duplicated(prevdat$id)],][1:3,-1]+prevdat[prevdat$id%in%prevdat$id[duplicated(prevdat$id)],][4:6,-1]
prevdat=prevdat[!duplicated(prevdat$id),]
rownames(prevdat)=prevdat[,1]
prevdat=prevdat[,-1]

prevdat.use=as.matrix(prevdat[previds1,])
rownames(prevdat.use)=previds1
prevdat.use=(prevdat.use/matrix(colSums(prevdat.use,na.rm=T),nrow=nrow(prevdat.use),ncol=ncol(prevdat.use),byrow=T))[,2:14]
colnames(prevdat.use)=2:14
prevdat.use.scale=t(apply(prevdat.use,1,scale,center=T,scale=T))
colnames(prevdat.use.scale)=2:14
prevdat.use.orders=order(apply(t(apply(prevdat.use.scale[,1:13],1,rollmean,3)),1,which.max))

counts.old=dat.mat.fraction_spikein[previds2,2:14]
counts.old.scale=t(apply(counts.old/matrix(colSums(dat.mat.fraction_spikein[previds2,2:14]),nrow=length(previds2),ncol=13,byrow=T),1,scale,center=T,scale=T))
colnames(counts.old.scale)=colnames(counts.old)
counts.old.orders=order(apply(t(apply(counts.old.scale[,4:13],1,rollmean,3)),1,which.max))

counts.old.scale.melt=melt(counts.old.scale)
colnames(counts.old.scale.melt)=c("id","fraction","value")

prevdat.use.scale.melt=melt(prevdat.use.scale)
colnames(prevdat.use.scale.melt)=c("id","fraction","value")
prevdat.use.scale.melt$labels=factor(previds.names[prevdat.use.scale.melt$id],levels=previds.names[prevdat.use.orders])
prevdat.use.scale.melt$fraction=factor(prevdat.use.scale.melt$fraction,levels=2:14)
counts.old.scale.melt$labels=factor(constructs.use[as.character(counts.old.scale.melt$id),1])
levels(counts.old.scale.melt$labels)[60]="hHBB_Nluc_hHBB" #Not the same barcode construct, but to plot it anyway
levels(counts.old.scale.melt$labels)[104]="TEV_Nluc_hHBA1" #Not the same barcode construct, but to plot it anyway
counts.old.scale.melt$labels=factor(as.character(counts.old.scale.melt$labels),levels=previds.names[prevdat.use.orders])
counts.old.scale.melt$fraction=factor(counts.old.scale.melt$fraction,levels=2:14)

counts.old.scale.melt=subset(counts.old.scale.melt,!is.na(as.character(labels)))
prevdat.use.scale.melt.empty=subset(prevdat.use.scale.melt,as.character(labels)%in%levels(counts.old.scale.melt$labels)[!levels(counts.old.scale.melt$labels)%in%counts.old.scale.melt$labels])
prevdat.use.scale.melt.empty$value=NA
counts.old.scale.melt=rbind(counts.old.scale.melt,prevdat.use.scale.melt.empty)

#233x fractions, %mRNA    
v.fraction.scale.melt.old=subset(v.fraction.scale.melt,constructs.use[as.character(v.fraction.scale.melt$id),"Designer"]=="Barna")
levels(v.fraction.scale.melt.old$labels)[165]="hHBB_Nluc_hHBB" #Not the same barcode construct, but to plot it anyway
levels(v.fraction.scale.melt.old$labels)[224]="TEV_Nluc_hHBA1" #Not the same barcode construct, but to plot it anyway
v.fraction.scale.melt.old$labels=factor(as.character(v.fraction.scale.melt.old$labels),levels=previds.names[prevdat.use.orders])
v.fraction.scale.melt.old=subset(v.fraction.scale.melt.old,!is.na(as.character(labels)))
v.fraction.scale.melt.old=rbind(v.fraction.scale.melt.old[,c("id","fraction","value","labels")],subset(prevdat.use.scale.melt.empty,as.numeric(fraction)%in%c(2:12)))

#233x fractions, %mRNA weighted
v.fraction.scale.weight.melt.old=subset(v.fraction.scale.weight.melt,constructs.use[as.character(v.fraction.scale.weight.melt$id),"Designer"]=="Barna")
levels(v.fraction.scale.weight.melt.old$labels)[165]="hHBB_Nluc_hHBB" #Not the same barcode construct, but to plot it anyway
levels(v.fraction.scale.weight.melt.old$labels)[224]="TEV_Nluc_hHBA1" #Not the same barcode construct, but to plot it anyway
v.fraction.scale.weight.melt.old$labels=factor(as.character(v.fraction.scale.weight.melt.old$labels),levels=previds.names[prevdat.use.orders])
v.fraction.scale.weight.melt.old=subset(v.fraction.scale.weight.melt.old,!is.na(as.character(labels)))
v.fraction.scale.weight.melt.old=rbind(v.fraction.scale.weight.melt.old[,c("id","fraction","value","labels")],subset(prevdat.use.scale.melt.empty,as.numeric(fraction)%in%c(2:12)))

fig.fractions.old=ggplot(counts.old.scale.melt,aes(x=fraction,y=labels,fill=(value),colour=(value)))+theme_classic()+geom_tile()+scale_fill_viridis(name="% mRNA",option="magma",oob=squish,na.value="darkgrey",limits=c(-2,2))+scale_colour_viridis(name="% mRNA",option="magma",oob=squish,na.value="darkgrey",limits=c(-2,2))+xlab("")+ylab("")
fig.fractions.prev=ggplot(prevdat.use.scale.melt,aes(x=fraction,y=labels,fill=(value),colour=(value)))+theme_classic()+geom_tile()+scale_fill_viridis(name="% mRNA",option="magma",oob=squish,na.value="darkgrey",limits=c(-2,2))+scale_colour_viridis(name="% mRNA",option="magma",oob=squish,na.value="darkgrey",limits=c(-2,2))+xlab("")+ylab("")
fig.vfractions.old=ggplot(v.fraction.scale.melt.old,aes(x=fraction,y=labels,fill=(value),colour=(value)))+theme_classic()+geom_tile()+scale_fill_viridis(name="% mRNA",option="magma",oob=squish,na.value="darkgrey")+scale_colour_viridis(name="% mRNA",option="magma",oob=squish,na.value="darkgrey")+xlab("")+ylab("")
fig.vfractions.weight.old=ggplot(v.fraction.scale.weight.melt.old,aes(x=fraction,y=labels,fill=(value),colour=(value)))+theme_classic()+geom_tile()+scale_fill_viridis(name="% mRNA",option="magma",oob=squish,na.value="darkgrey")+scale_colour_viridis(name="% mRNA",option="magma",oob=squish,na.value="darkgrey")+xlab("")+ylab("")

fig.compare=cowplot::plot_grid(fig.fractions.prev,fig.fractions.old,fig.vfractions.old,fig.vfractions.weight.old,ncol=4,align="hv",axis="tbr")
pdf("./plots/64x_comparison.pdf",width=20,height=8)
print(fig.compare)
dev.off()

#Short in cell stability
icstability_short_spikein=subset(samples,(Experiment=="stability")&(Spikein==1)&(Group=="In-cell")&(Type=="Short"))$index
icstability_short_spikein_labels=samples[icstability_short_spikein,"Time"]

dat.mat.icstability_short_spikein=dat.mat[,icstability_short_spikein]
colnames(dat.mat.icstability_short_spikein)=icstability_short_spikein_labels
dat.mat.icstability_short_spikein.log2=log2(dat.mat.icstability_short_spikein)

icstability_short.design=model.matrix(~0+factor(icstability_short_spikein_labels))
colnames(icstability_short.design)=unlist(lapply(strsplit(colnames(icstability_short.design),")"),"[[",2))

dat.mat.icstability_short_spikein.log2.scaled=dat.mat.icstability_short_spikein.log2
for (i in 1:ncol(dat.mat.icstability_short_spikein.log2))
{
  spikeinfit=lm(dat.mat.icstability_short_spikein.log2[usespikeins,c(1,i)])
  dat.mat.icstability_short_spikein.log2.scaled[,i]=predict.lm(spikeinfit,newdata=dat.mat.icstability_short_spikein.log2[,c(1,i)])
}
rownames(dat.mat.icstability_short_spikein.log2.scaled)=rownames(dat.mat.icstability_short_spikein.log2)
v.icstability_short=voom(2^(dat.mat.icstability_short_spikein.log2.scaled[!rownames(dat.mat.icstability_short_spikein.log2.scaled)%in%spikeins,]),icstability_short.design,lib.size=rep(1,15))

icstability_short.fit=lmFit(v.icstability_short,icstability_short.design)

icstability_short.fit.coef=icstability_short.fit$coefficients
icstability_short.fit.se=icstability_short.fit$stdev.unscaled*icstability_short.fit$sigma
icstability_short.fit.upper=icstability_short.fit.coef+icstability_short.fit.se
icstability_short.fit.lower=icstability_short.fit.coef-icstability_short.fit.se

icstability_short.fit.coef.scale=icstability_short.fit.coef-matrix(apply(icstability_short.fit.coef,1,max),nrow=nrow(icstability_short.fit.coef),ncol=ncol(icstability_short.fit.coef),byrow=F)
icstability_short.fit.upper.scale=icstability_short.fit.upper-matrix(apply(icstability_short.fit.upper,1,max),nrow=nrow(icstability_short.fit.upper),ncol=ncol(icstability_short.fit.upper),byrow=F)
icstability_short.fit.lower.scale=icstability_short.fit.lower-matrix(apply(icstability_short.fit.lower,1,max),nrow=nrow(icstability_short.fit.lower),ncol=ncol(icstability_short.fit.lower),byrow=F)

icstability_short.fit.coef.scale.melt=melt(icstability_short.fit.coef.scale)
colnames(icstability_short.fit.coef.scale.melt)=c("id","time","value")
icstability_short.fit.coef.scale.melt$upper=melt(icstability_short.fit.upper.scale)$value
icstability_short.fit.coef.scale.melt$lower=melt(icstability_short.fit.lower.scale)$value
icstability_short.fit.coef.scale.melt$time=factor(icstability_short.fit.coef.scale.melt$time)
icstability_short.fit.coef.scale.melt$totalcount=log2(rowSums(2^icstability_short.fit.coef)[icstability_short.fit.coef.scale.melt$id])

#Fraction mRNA remamining relative to time zero
write.table(2^icstability_short.fit.coef.scale,file="./tables/short_spikein_incellstability_percent.tsv",row.names=T,col.names=T,quote=F,sep="\t")

icstability_short.fit.coef.scale.melt$id=factor(icstability_short.fit.coef.scale.melt$id,levels=rownames(icstability_short.fit.coef.scale)[order(icstability_short.fit.coef.scale[,3])])
fig.icstability.short=ggplot(subset(icstability_short.fit.coef.scale.melt),aes(x=time,y=id,fill=value,colour=value))+theme_classic()+geom_tile()+scale_fill_viridis(option="magma",oob=squish,na.value="darkgrey")+scale_colour_viridis(option="magma",oob=squish,na.value="darkgrey")+xlab("")+ylab("")+scale_y_discrete(labels=constructs.use[levels(icstability_short.fit.coef.scale.melt$id),1])

pdf("./plots/short_spikein_incellstability.pdf",width=9,height=22)
print(fig.icstability.short)
dev.off()

fig.icstability.short.line=ggplot(subset(icstability_short.fit.coef.scale.melt),aes(x=time,y=2^value,group=id,colour=totalcount))+theme_classic()+geom_line(alpha=0.33)+scale_colour_viridis(name="Total abundance",option="magma")+xlab("Time")+ylab("Log abundance")

pdf("./plots/short_spikein_incellstability_line.pdf",width=9,height=22)
print(fig.icstability.short)
dev.off()

#In-vitro stability, long
ivstability_long_spikein=subset(samples,(Experiment=="stability")&(Spikein==1)&(Group=="In-vitro")&(Type=="Long"))$index
ivstability_long_spikein_labels=samples[ivstability_long_spikein,"Sample_Name"]

dat.mat.ivstability_long_spikein=dat.mat[,ivstability_long_spikein]
colnames(dat.mat.ivstability_long_spikein)=ivstability_long_spikein_labels
long_used=constructs.use[unique(c(grep("acatttgcttctgacacaactgtgttcac",tolower(constructs.use$Full.sequence)),grep("ATGGCCGTTTACCCATACGATGTTCCTGAC",toupper(constructs.use$Full.sequence)))),]$SequenceID #remove no HBB or Nluc
dat.mat.ivstability_long_spikein=dat.mat.ivstability_long_spikein[long_used,]

dat.mat.ivstability_long_spikein.log2=log2(dat.mat.ivstability_long_spikein)

ivstability_long.design=model.matrix(~0+factor(ivstability_long_spikein_labels))
colnames(ivstability_long.design)=unlist(lapply(strsplit(colnames(ivstability_long.design),")"),"[[",2))
ivstability_long.design=ivstability_long.design[,ivstability_long_spikein_labels]

dat.mat.ivstability_long_spikein.log2.scaled=dat.mat.ivstability_long_spikein.log2
for (i in 1:ncol(dat.mat.ivstability_long_spikein.log2))
{
  spikeinfit=lm(dat.mat.ivstability_long_spikein.log2[usespikeins,c(1,i)])
  dat.mat.ivstability_long_spikein.log2.scaled[,i]=predict.lm(spikeinfit,newdata=dat.mat.ivstability_long_spikein.log2[,c(1,i)])
}
rownames(dat.mat.ivstability_long_spikein.log2.scaled)=rownames(dat.mat.ivstability_long_spikein.log2)
v.ivstability_long=voom(2^(dat.mat.ivstability_long_spikein.log2.scaled[!rownames(dat.mat.ivstability_long_spikein.log2.scaled)%in%spikeins,]),ivstability_long.design,lib.size=rep(1,10))

ivstability_long.fit=lmFit(v.ivstability_long,ivstability_long.design)
ivstability_long.fit.coef=ivstability_long.fit$coefficients
ivstability_long.fit.se=ivstability_long.fit$stdev.unscaled*ivstability_long.fit$sigma
ivstability_long.fit.upper=ivstability_long.fit.coef+ivstability_long.fit.se
ivstability_long.fit.lower=ivstability_long.fit.coef-ivstability_long.fit.se

ivstability_long.fit.coef.dr=list()
ivstability_long.fit.coef.init=list()
for(i in rownames(ivstability_long.fit.coef))
{
  ivstability_long.fit.coef.dr.fit=lm(log(2^ivstability_long.fit.coef[i,3:6])~c(0,0.5,1,2,3,4,5,6,16,24)[3:6])
  ivstability_long.fit.coef.dr[i]=(-(ivstability_long.fit.coef.dr.fit$coefficients[2]))
  ivstability_long.fit.coef.init[i]=(ivstability_long.fit.coef.dr.fit$coefficients[1])
}
ivstability_long.fit.coef.dr=unlist(ivstability_long.fit.coef.dr)
ivstability_long.fit.coef.init=unlist(ivstability_long.fit.coef.init)

ivstability_long.fit.coef.scale=ivstability_long.fit.coef-matrix(apply(ivstability_long.fit.coef,1,max),nrow=nrow(ivstability_long.fit.coef),ncol=ncol(ivstability_long.fit.coef),byrow=F)
ivstability_long.fit.upper.scale=NA
ivstability_long.fit.lower.scale=NA

#plot(c(25,2.5,0.25,0.025)/c(953,897,897,882),dat.mat.ivstability_long_spikein.log2[spikeins[c(4,1,2,3)],1])

ivstability_long.fit.coef.scale.melt=melt(ivstability_long.fit.coef.scale)
colnames(ivstability_long.fit.coef.scale.melt)=c("id","time","value")
ivstability_long.fit.coef.scale.melt$upper=melt(ivstability_long.fit.upper.scale)$value
ivstability_long.fit.coef.scale.melt$lower=melt(ivstability_long.fit.lower.scale)$value
ivstability_long.fit.coef.scale.melt$time=factor(ivstability_long.fit.coef.scale.melt$time,levels=ivstability_long_spikein_labels)
ivstability_long.fit.coef.scale.melt$totalcount=ivstability_long.fit.coef[ivstability_long.fit.coef.scale.melt$id,1]

#Fraction of mRNA relative to time zero
write.table(2^ivstability_long.fit.coef.scale,file="./tables/long_spikein_invitrostability_percent.tsv",row.names=T,col.names=T,quote=F,sep="\t")

ivstability_long.fit.coef.scale.melt$id=factor(ivstability_long.fit.coef.scale.melt$id,levels=rownames(ivstability_long.fit.coef.scale)[order(ivstability_long.fit.coef.dr)])

fig.ivstability.long=ggplot(subset(ivstability_long.fit.coef.scale.melt),aes(x=time,y=id,fill=value,colour=value))+theme_classic()+geom_tile()+scale_fill_viridis(option="magma",oob=squish,na.value="darkgrey")+scale_colour_viridis(option="magma",oob=squish,na.value="darkgrey")+xlab("")+ylab("")+scale_y_discrete(labels=constructs.use[levels(ivstability_long.fit.coef.scale.melt$id),1])
pdf("./plots/long_spikein_invitrostability.pdf",width=10,height=20)
print(fig.ivstability.long)
dev.off()

fig.ivstability.long.line=ggplot(subset(ivstability_long.fit.coef.scale.melt),aes(x=time,y=value,group=id,colour=totalcount))+theme_classic()+geom_line(alpha=0.33)+scale_colour_viridis(name="Total abundance",option="magma")+xlab("Time")+ylab("Log abundance")
pdf("./plots/long_spikein_invitrostability_line.pdf",width=6,height=6)
print(fig.ivstability.long.line)
dev.off()

#short stability, in vitro
ivstability_short_spikein=subset(samples,(Experiment=="stability")&(Spikein==1)&(Group=="In-vitro")&(Type=="Short"))$index
ivstability_short_spikein_labels=samples[ivstability_short_spikein,"Sample_Name"]

dat.mat.ivstability_short_spikein=dat.mat[,ivstability_short_spikein]
colnames(dat.mat.ivstability_short_spikein)=ivstability_short_spikein_labels
dat.mat.ivstability_short_spikein.log2=log2(dat.mat.ivstability_short_spikein)

ivstability_short.design=model.matrix(~0+factor(ivstability_short_spikein_labels))
colnames(ivstability_short.design)=unlist(lapply(strsplit(colnames(ivstability_short.design),")"),"[[",2))
ivstability_short.design=ivstability_short.design[,ivstability_short_spikein_labels]

dat.mat.ivstability_short_spikein.log2.scaled=dat.mat.ivstability_short_spikein.log2
for (i in 1:ncol(dat.mat.ivstability_short_spikein.log2))
{
  spikeinfit=lm(dat.mat.ivstability_short_spikein.log2[usespikeins,c(1,i)])
  dat.mat.ivstability_short_spikein.log2.scaled[,i]=predict.lm(spikeinfit,newdata=dat.mat.ivstability_short_spikein.log2[,c(1,i)])
}
rownames(dat.mat.ivstability_short_spikein.log2.scaled)=rownames(dat.mat.ivstability_short_spikein.log2)
v.ivstability_short=voom(2^(dat.mat.ivstability_short_spikein.log2.scaled[!rownames(dat.mat.ivstability_short_spikein.log2.scaled)%in%spikeins,]),ivstability_short.design,lib.size=rep(1,10))

ivstability_short.fit=lmFit(v.ivstability_short,ivstability_short.design)
ivstability_short.fit.coef=ivstability_short.fit$coefficients
ivstability_short.fit.se=ivstability_short.fit$stdev.unscaled*ivstability_short.fit$sigma
ivstability_short.fit.upper=ivstability_short.fit.coef+ivstability_short.fit.se
ivstability_short.fit.lower=ivstability_short.fit.coef-ivstability_short.fit.se

ivstability_short.fit.coef.scale=ivstability_short.fit.coef-matrix(apply(ivstability_short.fit.coef,1,max),nrow=nrow(ivstability_short.fit.coef),ncol=ncol(ivstability_short.fit.coef),byrow=F)
ivstability_short.fit.upper.scale=NA
ivstability_short.fit.lower.scale=NA

ivstability_short.fit.coef.scale.melt=melt(ivstability_short.fit.coef.scale)
colnames(ivstability_short.fit.coef.scale.melt)=c("id","time","value")
ivstability_short.fit.coef.scale.melt$upper=melt(ivstability_short.fit.upper.scale)$value
ivstability_short.fit.coef.scale.melt$lower=melt(ivstability_short.fit.lower.scale)$value
ivstability_short.fit.coef.scale.melt$time=factor(ivstability_short.fit.coef.scale.melt$time,levels=ivstability_short_spikein_labels)
ivstability_short.fit.coef.scale.melt$totalcount=log2(rowSums(2^ivstability_short.fit.coef)[ivstability_short.fit.coef.scale.melt$id])

write.table(2^ivstability_short.fit.coef.scale,file="./tables/short_spikein_invitrostability_percent.tsv",row.names=T,col.names=T,quote=F,sep="\t")

ivstability_short.fit.coef.scale.melt$id=factor(ivstability_short.fit.coef.scale.melt$id,levels=rownames(ivstability_short.fit.coef.scale)[order(rowMeans(apply(ivstability_short.fit.coef.scale,2,rank)[,-1]))])

fig.ivstability.short=ggplot(subset(ivstability_short.fit.coef.scale.melt),aes(x=time,y=id,fill=value,colour=value))+theme_classic()+geom_tile()+scale_fill_viridis(option="magma",oob=squish,na.value="darkgrey")+scale_colour_viridis(option="magma",oob=squish,na.value="darkgrey")+xlab("")+ylab("")+scale_y_discrete(labels=constructs.use[levels(ivstability_short.fit.coef.scale.melt$id),1])
pdf("./plots/short_spikein_invitrostability.pdf",width=6,height=6)
print(fig.ivstability.long)
dev.off()

fig.ivstability.short.line=ggplot(subset(ivstability_short.fit.coef.scale.melt),aes(x=time,y=2^value,group=id,colour=totalcount))+theme_classic()+geom_line(alpha=0.33)+scale_colour_viridis(name="Total abundance",option="magma")+xlab("Time")+ylab("Log abundance")
pdf("./plots/short_spikein_invitrostability_line.pdf",width=6,height=6)
print(fig.ivstability.short.line)
dev.off()


#Load individual luciferase data and format them
#Experiment 8
luc_exp8_pest=read.csv("./exp8_nluc_pest.csv",header=T,sep="\t")
luc_exp8_pest=data.frame(luc_exp8_pest %>% group_by(construct,time) %>% summarise(luc.mean_pest=mean(log2(Nluc/Fluc),na.rm=T),luc.sd_pest=sd(log2(Nluc/Fluc),na.rm=T), rna.mean_pest=mean(log2(rna),na.rm=T), rna.sd_pest=sd(log2(rna),na.rm=T)))
luc_exp8_pest.constructs=subset(constructs.use,Experiment=="NlucP")[,c(1,4)]
rownames(luc_exp8_pest.constructs)=luc_exp8_pest.constructs[,1]
luc_exp8_pest.constructs$constructs.np=c("320047B1","120019B1","120018B1","120020B1","220027B1","220025B1","220022B1","120009B1","320053B1","520079B1","120003B1","220030B1")
luc_exp8_pest$construct.id=luc_exp8_pest.constructs[luc_exp8_pest$construct,2]
luc_exp8_pest$construct.id.np=luc_exp8_pest.constructs[luc_exp8_pest$construct,3]

#Experiment 7
luc_exp7=read.csv("./exp7_nluc.csv",header=T,sep="\t")
luc_exp7=data.frame(luc_exp7 %>% group_by(construct,time) %>% summarise(luc.mean_exp7=mean(log2(Nluc/Fluc),na.rm=T),luc.sd_exp7=sd(log2(Nluc/Fluc),na.rm=T), rna.mean_exp7=mean(log2(rna),na.rm=T), rna.sd_exp7=sd(log2(rna),na.rm=T)))
luc_exp7.constructs=c("320047B1","120019B1","120018B1","120020B1","220027B1","220025B1","220022B1","120009B1","320053B1","520079B1","120003B1","220030B1")
names(luc_exp7.constructs)=unlist(lapply(strsplit(constructs.use[luc_exp7.constructs,1],"\\."),"[[",1))
luc_exp7$construct.id=luc_exp7.constructs[luc_exp7$construct]
luc_exp7$construct.id.np=luc_exp7.constructs[luc_exp7$construct]

#Experiment 10
luc_exp10=read.csv("./exp10_nluc.csv",header=T,sep="\t")
luc_exp10=data.frame(luc_exp10 %>% group_by(construct,time) %>% summarise(luc.mean_exp10=mean(log2(Nluc/Fluc),na.rm=T),luc.sd_exp10=sd(log2(Nluc/Fluc),na.rm=T)))
luc_exp10.constructs=c("320047B1","120019B1","120018B1","120020B1","220027B1","220025B1","220022B1","120009B1","320053B1","520079B1","120003B1","220030B1")
names(luc_exp10.constructs)=unlist(lapply(strsplit(constructs.use[luc_exp10.constructs,1],"\\."),"[[",1))
luc_exp10$construct.id=luc_exp10.constructs[luc_exp10$construct]
luc_exp10$construct.id.np=luc_exp10.constructs[luc_exp10$construct]

#Joined table
luc_joined=Reduce(function(x,y) plyr::join(x,y, by=c("construct.id.np","time"),type="full",match="all"), list(luc_exp7,luc_exp10,luc_exp8_pest))

#6hr data
luc_joined.6=subset(luc_joined,time==6)
luc_joined.6$decay_rate=icstability_long.fit.coef.dr[luc_joined.6$construct.id.np]
luc_joined.6$weighted_ribosomes=rowSums(v.fraction.scale.weight[luc_joined.6$construct.id,])
luc_joined.6$length=constructs.use[luc_joined.6$construct.id,"Length.to.order"]

luc_joined.6.poly=melt(v.fraction.scale[luc_joined.6$construct.id,])
luc_joined.6.poly$labels=constructs.use[as.character(luc_joined.6.poly$Var1),1]
fig.luc.poly=ggplot(luc_joined.6.poly,aes(x=Var2,y=value,group=labels,colour=labels))+theme_classic()+geom_line()+scale_colour_iwanthue()

#differential equation system
diff_func=function(t,state,params) {with(as.list(c(state,params)), {
  dmrna=(-deg_rate_mrna)*mrna
  dprot=translation*mrna-deg_rate_prot*prot
  list(c(dmrna,dprot))
})}

#18+-11min NlucP; Nluc >6h (no degradation seen after 6h...)
deg_rate_prot=log(2)/100000
#deg_rate_prot=log(2)/0.33

cz_prot_noribo=list()
cz_prot=list()
cz_mrna=list()
times=seq(0,64,by=0.01)
for(i in rownames(icstability_long.fit.coef))
{
  inv_mrna_length=1/as.numeric(constructs.use[i,"Length.to.order"])
  num_ribosome=as.numeric(rowSums(v.fraction.scale.weight)[i])
  
  decay_rate_mrna=as.numeric(icstability_long.fit.coef.dr[i])
  decay_rate_prot=deg_rate_prot
  
  params=c(deg_rate_mrna=decay_rate_mrna,translation=num_ribosome,deg_rate_prot=decay_rate_prot)
  params_noribo=c(deg_rate_mrna=decay_rate_mrna,translation=1,deg_rate_prot=decay_rate_prot)
  state=c(mrna=inv_mrna_length,prot=0)
  
  solved=ode(y=state,times=times,func=diff_func,parms=params)
  solved_noribo=ode(y=state,times=times,func=diff_func,parms=params_noribo)
  
  cz_prot[[i]]=solved[,3]
  cz_prot_noribo[[i]]=solved_noribo[,3]
  cz_mrna[[i]]=solved[,2]
}  

cz_prot.mat=do.call(cbind,cz_prot)[,luc_joined.6$construct.id.np]
cz_prot.all.mat=t(do.call(cbind,cz_prot))
colnames(cz_prot.all.mat)=times
cz_prot_noribo.all.mat=t(do.call(cbind,cz_prot_noribo))
colnames(cz_prot_noribo.all.mat)=times
cz_mrna.mat=do.call(cbind,cz_mrna)[,luc_joined.6$construct.id.np]
cz.melt=melt(cz_prot.mat)
cz.melt$mrna=melt(cz_mrna.mat)$value
cz.melt$Var2=constructs.use[as.character(cz.melt$Var2),1]

#Predicted expression vs time
pdf("./plots/predicted_prot.pdf",width=8,height=6)
ggplot(cz.melt,aes(x=times[Var1],y=value,group=Var2,colour=Var2))+theme_classic()+geom_line()+scale_colour_iwanthue()+ylab("Predicted expression")+xlab("Time")
dev.off()

#Predicted mRNA vs time
pdf("./plots/predicted_mrna.pdf",width=8,height=6)
ggplot(cz.melt,aes(x=times[Var1],y=mrna,group=Var2,colour=Var2))+theme_classic()+geom_line()+scale_colour_iwanthue()+ylab("Predicted mRNA")+xlab("Time")
#dev.off()

#Predicted vs measured expression
#luc_joined.6$pred=subset(cz.melt,Var1==which.max(apply(cz_prot.mat,1,cor,2^luc_joined.6$luc.mean_pest)))$value
luc_joined.6$pred=subset(cz.melt,Var1==which.max(apply(cz_prot.mat,1,cor,2^luc_joined.6$luc.mean_exp7)))$value
#times[which.max(apply(cz_prot.mat,1,cor,2^luc_joined.6$luc.mean_pest))]
#times[which.max(apply(cz_prot.mat,1,cor,2^luc_joined.6$luc.mean_exp7))]
pdf("./plots/predicted_measured.pdf",width=8,height=6)
ggplot(subset(luc_joined.6),aes(y=2^luc.mean_exp7,x=(pred)))+theme_classic()+geom_point()+geom_text_repel(aes(label=construct))+xlab("Predicted")+ylab("Exp10 Nluc/Fluc")
#ggplot(subset(luc_joined.6),aes(y=2^luc.mean_pest,x=(pred)))+theme_classic()+geom_point()+geom_text_repel(aes(label=construct))+xlab("Predicted")+ylab("Exp10 Nluc/Fluc")
dev.off()

luc_joined.24=subset(luc_joined,time==24)
luc_joined.24$pred=subset(cz.melt,Var1==which.max(apply(cz_prot.mat,1,cor,2^luc_joined.24$luc.mean_exp7)))$value
#luc_joined.24$pred=subset(cz.melt,Var1==which.max(apply(cz_prot.mat,1,cor,2^luc_joined.24$luc.mean_pest)))$value
times[which.max(apply(cz_prot.mat,1,cor,2^luc_joined.24$luc.mean_exp7))]
#times[which.max(apply(cz_prot.mat,1,cor,2^luc_joined.24$luc.mean_pest))]
#ggplot(subset(luc_joined.24),aes(y=2^luc.mean_exp7,x=(pred)))+theme_classic()+geom_point()+geom_text_repel(aes(label=construct))+xlab("Predicted")+ylab("Exp10 Nluc/Fluc")
#ggplot(subset(luc_joined.24),aes(y=2^luc.mean_pest,x=(pred)))+theme_classic()+geom_point()+geom_text_repel(aes(label=construct))+xlab("Predicted")+ylab("Exp10 Nluc/Fluc")

luc_joined.12=subset(luc_joined,time==12)
luc_joined.12$pred=subset(cz.melt,Var1==which.max(apply(cz_prot.mat,1,cor,2^luc_joined.12$luc.mean_exp7)))$value
#times[which.max(apply(cz_prot.mat,1,cor,2^luc_joined.12$luc.mean_exp7))]
#ggplot(subset(luc_joined.12),aes(y=2^luc.mean_exp7,x=(pred)))+theme_classic()+geom_point()+geom_text_repel(aes(label=construct))+xlab("Predicted")+ylab("Exp10 Nluc/Fluc")


#IDs ordred by design type and fraction
plotids=rownames(v.fraction.scale)

#Predicted prots
cz.plot=data.frame(ids=names(cz_prot),pred=cz_prot.all.mat[,ncol(cz_prot.all.mat)],noribo_pred=cz_prot_noribo.all.mat[,ncol(cz_prot_noribo.all.mat)])
cz.plot$pred.use=cz.plot$pred
cz.plot[intersect(subset(constructs.use,!Designer=="Barna")$SequenceID,cz.plot$ids),"pred.use"]=cz.plot[intersect(subset(constructs.use,!Designer=="Barna")$SequenceID,cz.plot$ids),"noribo_pred"]

cz.plot.missing=plotids[!plotids%in%rownames(cz.plot)]
cz.plot=rbind(cz.plot,data.frame(ids=cz.plot.missing,pred=NA,noribo_pred=NA,pred.use=NA))

plotids.ordered=NULL
plotids.groups=levels(constructs.use$Experiment)

#Ordering by groups
for (group in plotids.groups)
{
  group_seqid=subset(constructs.use[plotids,],Experiment==group)$SequenceID
  
  #Order by polysome hclust
  #group_order=hclust(dist(v.fraction.scale[group_seqid,]),"average")$order
  
  #Order by in vitro stability
  #group_seqid_use=intersect(rownames(ivstability_long.fit.coef.scale),group_seqid)
  #group_order=rep(NA,length(group_seqid))
  #names(group_order)=group_seqid
  #group_order[group_seqid_use]=order(rowMeans(apply(ivstability_long.fit.coef.scale[,2:6],2,rank))[group_seqid_use])
  #group_order[group_seqid_use]=order(ivstability_long.fit.coef.dr[group_seqid_use])
  
  #Order by predicted
  group_order=order(cz.plot[group_seqid,"pred.use"],decreasing=T)
  
  group_seqid=rev(group_seqid[group_order])
  plotids.ordered=c(plotids.ordered,group_seqid)
}
plotids=plotids.ordered
plotids=plotids[!is.na(plotids)]

icstability_long.fit.dr.plot=data.frame(id=plotids,labels=constructs.use[plotids,1],value=icstability_long.fit.coef.dr[plotids])
icstability_long.fit.dr.plot$id=factor(as.character(icstability_long.fit.dr.plot$id),levels=plotids)
icstability_long.fit.dr.plot$labels=factor(icstability_long.fit.dr.plot$labels,levels=constructs.use[plotids,1])

fig.icstability.long.plot=ggplot(icstability_long.fit.dr.plot,aes(y=labels,x=value))+theme_classic()+geom_bar(stat="identity",fill="black")+scale_x_continuous(expand=c(0,0))+ylab("")+xlab("In-cell\nRNA decay")+theme(legend.position="bottom",axis.text.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank())

ivstability_long.fit.dr.plot=data.frame(id=plotids,labels=constructs.use[plotids,1],value=ivstability_long.fit.coef.dr[plotids])
ivstability_long.fit.dr.plot$id=factor(as.character(ivstability_long.fit.dr.plot$id),levels=plotids)
ivstability_long.fit.dr.plot$labels=factor(ivstability_long.fit.dr.plot$labels,levels=constructs.use[plotids,1])

fig.ivstability.long.plot=ggplot(ivstability_long.fit.dr.plot,aes(y=labels,x=value))+theme_classic()+geom_bar(stat="identity",fill="black")+scale_x_continuous(expand=c(0,0))+ylab("")+xlab("In-vitro\nRNA decay")+theme(legend.position="bottom",axis.text.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank())


v.fraction.scale.weight.tr.plot=data.frame(id=plotids,labels=constructs.use[plotids,1],value=rowSums(v.fraction.scale.weight)[plotids])
v.fraction.scale.weight.tr.plot$value=v.fraction.scale.weight.tr.plot$value-min(v.fraction.scale.weight.tr.plot$value,na.rm=T)
v.fraction.scale.weight.tr.plot$id=factor(as.character(v.fraction.scale.weight.tr.plot$id),levels=plotids)
v.fraction.scale.weight.tr.plot$labels=factor(v.fraction.scale.weight.tr.plot$labels,levels=constructs.use[plotids,1])

fig.fraction.tr.plot=ggplot(v.fraction.scale.weight.tr.plot,aes(y=labels,x=value))+theme_classic()+geom_bar(stat="identity",fill="black")+scale_x_continuous(expand=c(0,0))+ylab("")+xlab("Ribo number\n(-min)")+theme(legend.position="bottom",axis.text.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank())

cz.plot$labels=factor(constructs.use[cz.plot$ids,1],levels=constructs.use[plotids,1])

fig.cz.plot=ggplot(cz.plot,aes(y=labels,x=pred.use))+theme_classic()+geom_bar(stat="identity",fill="black")+scale_x_continuous(expand=c(0,0))+ylab("")+xlab("Predicted")+theme(legend.position="bottom",axis.text.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank())

ivstability_long.fit.coef.scale.melt.plot=subset(ivstability_long.fit.coef.scale.melt,as.character(id)%in%plotids)
ivstability_long.fit.coef.scale.melt.plot$id=factor(as.character(ivstability_long.fit.coef.scale.melt.plot$id),levels=plotids)
ivstability_long.fit.coef.scale.melt.plot=ivstability_long.fit.coef.scale.melt.plot[!is.na(ivstability_long.fit.coef.scale.melt.plot$id),]

ivstability_long.fit.coef.scale.melt.plot.missing=plotids[!plotids%in%unique(as.character(ivstability_long.fit.coef.scale.melt.plot$id))]
ivstability_long.fit.coef.scale.melt.plot=rbind(ivstability_long.fit.coef.scale.melt.plot,data.frame(id=rep(ivstability_long.fit.coef.scale.melt.plot.missing,each=length(levels(ivstability_long.fit.coef.scale.melt.plot$time))),time=rep(levels(ivstability_long.fit.coef.scale.melt.plot$time),length(ivstability_long.fit.coef.scale.melt.plot.missing)),value=NA,upper=NA,lower=NA))

ivstability_long.fit.coef.scale.melt.plot$labels=factor(constructs.use[as.character(ivstability_long.fit.coef.scale.melt.plot$id),1],levels=constructs.use[plotids,1])
levels(ivstability_long.fit.coef.scale.melt.plot$time)=as.character(c(0,0.5,1,2,3,4,5,6,16,24))

fig.ivstability.long.scale.plot=ggplot(subset(ivstability_long.fit.coef.scale.melt.plot),aes(x=time,y=labels,fill=(value),colour=(value)))+theme_classic()+geom_tile()+scale_fill_viridis(option="magma",oob=squish,na.value="darkgrey",name="% mRNA")+scale_colour_viridis(option="magma",oob=squish,na.value="darkgrey",name="% mRNA")+xlab("Time")+ylab("")+theme(axis.text.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.position="bottom")

v.fraction.scale.melt.plot=subset(v.fraction.scale.melt,as.character(id)%in%plotids)
v.fraction.scale.melt.plot$id=factor(as.character(v.fraction.scale.melt.plot$id),levels=plotids)
v.fraction.scale.melt.plot=v.fraction.scale.melt.plot[!is.na(v.fraction.scale.melt.plot$id),]
v.fraction.scale.melt.plot$labels=factor(constructs.use[as.character(v.fraction.scale.melt.plot$id),1],levels=constructs.use[plotids,1])

fig.fraction.plot=ggplot(subset(v.fraction.scale.melt.plot,fraction%in%c(2:11)),aes(x=fraction,y=labels,fill=(value),colour=(value)))+theme_classic()+geom_tile()+scale_fill_viridis(option="magma",oob=squish,na.value="darkgrey",name="% mRNA")+scale_colour_viridis(option="magma",oob=squish,na.value="darkgrey",name="% mRNA")+xlab("Fractions 2-11")+ylab("")+theme(axis.text.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.position="bottom")

v.fraction.scale.weight.melt.plot=subset(v.fraction.scale.weight.melt,as.character(id)%in%plotids)
v.fraction.scale.weight.melt.plot$id=factor(as.character(v.fraction.scale.weight.melt.plot$id),levels=levels(v.fraction.scale.weight.melt$id)[levels(v.fraction.scale.weight.melt$id)%in%plotids])
v.fraction.scale.weight.melt.plot=v.fraction.scale.weight.melt.plot[!is.na(v.fraction.scale.weight.melt.plot$id),]
v.fraction.scale.weight.melt.plot$labels=factor(constructs.use[as.character(v.fraction.scale.weight.melt.plot$id),1],levels=constructs.use[plotids,1])

fig.fraction.weight.plot=ggplot(subset(v.fraction.scale.weight.melt.plot,fraction%in%c(2:11)),aes(x=fraction,y=labels,fill=(value),colour=(value)))+theme_classic()+geom_tile()+scale_fill_viridis(option="magma",oob=squish,na.value="darkgrey",name="Weighted % mRNA")+scale_colour_viridis(option="magma",oob=squish,na.value="darkgrey",name="Weighted % mRNA")+xlab("Fractions 2-11")+ylab("")+theme(axis.text.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.position="bottom")


v.fraction.zscale.melt.plot=subset(v.fraction.zscale.melt,as.character(id)%in%plotids)
v.fraction.zscale.melt.plot$id=factor(as.character(v.fraction.zscale.melt.plot$id),levels=levels(v.fraction.zscale.melt$id)[levels(v.fraction.zscale.melt$id)%in%plotids])
v.fraction.zscale.melt.plot=v.fraction.zscale.melt.plot[!is.na(v.fraction.zscale.melt.plot$id),]
v.fraction.zscale.melt.plot$labels=factor(constructs.use[as.character(v.fraction.zscale.melt.plot$id),1],levels=constructs.use[plotids,1])

fig.fraction.zscale.plot=ggplot(subset(v.fraction.zscale.melt.plot,fraction%in%c(2:11)),aes(x=fraction,y=labels,fill=(value),colour=(value)))+theme_classic()+geom_tile()+scale_fill_viridis(option="magma",oob=squish,na.value="darkgrey",name="Weighted % mRNA")+scale_colour_viridis(option="magma",oob=squish,na.value="darkgrey",name="Weighted % mRNA")+xlab("Fractions 2-11")+ylab("")+theme(axis.text.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.position="bottom")

#For labeling groups
groupdf=constructs.use[plotids,]
groupdf$SequenceID=factor(groupdf$SequenceID,levels=plotids)
groupdf$labels=factor(constructs.use[plotids,1],levels=constructs.use[plotids,1])
fig.groupbar=ggplot(groupdf,aes(x=1,y=labels,fill=Experiment,colour=Experiment))+theme_classic()+geom_tile()+scale_fill_iwanthue()+scale_colour_iwanthue()+scale_x_discrete(expand=c(0,0))+xlab("")+ylab("")+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.position="bottom")

fig.panel=cowplot::plot_grid(
  fig.groupbar,
  fig.fraction.plot,
  fig.fraction.weight.plot,
  fig.fraction.zscale.plot,
  fig.ivstability.long.scale.plot,
  fig.ivstability.long.plot,
  fig.icstability.long.plot,
  fig.fraction.tr.plot,
  fig.cz.plot,
  ncol=9,nrow=1,align="hv",axis="tbr",rel_widths=c(0.45,0.3,0.3,0.3,0.3,0.1,0.1,0.1,0.1))

pdf("./plots/data_panel_pred.pdf",width=36,height=24)
print(fig.panel)
dev.off()

