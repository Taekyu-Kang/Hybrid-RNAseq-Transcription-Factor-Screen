options(stringsAsFactors=F)
library(DESeq2)
library(pasilla)
library(Biobase)
library(GSVA)
library(ggplot2)
library(ggpubr)
library(EnvStats)
library(gplots)
library(psych)
library(corrplot)
library(rtracklayer)

#perform ANOVA with TPM data to determine those showing condition specific ASE
#read in TPM data, remove duplicates
f1_TPM_org <- read.table("/home/tkang/primary_ASE/all_f1_reps/python_test/2022-06-22_F1hybridTPM.txt",header=T,quote="",sep="\t")
f1_TPM_org <- f1_TPM_org[!is.na(f1_TPM_org$Genes),]
f1_TPM_org <- f1_TPM_org[!duplicated(f1_TPM_org$Genes),]
rownames(f1_TPM_org) <- f1_TPM_org$Genes
f1_TPM_org1 <- f1_TPM_org[,-1]

#remove those showing 0 in more than half the samples
f1dat <- f1_TPM_org1
f1dat[f1dat == 0] <- NA
f1filt <- apply(f1dat,1,function(x) sum(is.na(x)) >= 6)
f1filt1 <- unname(f1filt)
f1dat1 <- f1dat[!f1filt1,]
f1dat1[is.na(f1dat1)] <- 0

#remove those with high coefficient of variation of ASE across replicates
f1dat1$F1.1_clog2fc <- log2(f1dat1$PWK.CONT.1 / f1dat1$STF.CONT.1)
f1dat1$F1.1_xlog2fc <- log2(f1dat1$PWK.SEN.1 / f1dat1$STF.SEN.1)

f1dat1$F1.2_clog2fc <- log2(f1dat1$PWK.CONT.2 / f1dat1$STF.CONT.2)
f1dat1$F1.2_xlog2fc <- log2(f1dat1$PWK.SEN.2 / f1dat1$STF.SEN.2)

f1dat1$F1.4_clog2fc <- log2(f1dat1$PWK.CONT.3 / f1dat1$STF.CONT.3)
f1dat1$F1.4_xlog2fc <- log2(f1dat1$PWK.SEN.3 / f1dat1$STF.SEN.3)

f1dat2 <- f1dat1[,c(13:18)]

#remove those with inf or NA
f1dat2.1 <- f1dat2
f1dat2.2 <- do.call(data.frame,                      
                    lapply(f1dat2.1,
                           function(x) replace(x, is.infinite(x), NA)))
f1dat2.3 <- do.call(data.frame,                      
                    lapply(f1dat2.2,
                           function(x) replace(x, is.nan(x), NA)))
rownames(f1dat2.3) <- rownames(f1dat2.1)

f1dat2.3$cont_mean <- rowMeans(f1dat2.3[,c(1,3,5)],na.rm=T)
f1dat2.3$cont_sd <- rowSds(as.matrix(f1dat2.3[,c(1,3,5)]),na.rm=T)

f1dat2.3$sen_mean <- rowMeans(f1dat2.3[,c(2,4,6)],na.rm=T)
f1dat2.3$sen_sd <- rowSds(as.matrix(f1dat2.3[,c(2,4,6)]),na.rm=T)

f1dat3 <- f1dat2.3[,c(7:10)]
f1dat3[is.na(f1dat3)] <- 0

f1dat3$cont_cov <- f1dat3$cont_sd / f1dat3$cont_mean
f1dat3$sen_cov <- f1dat3$sen_sd / f1dat3$sen_mean

quantile(f1dat3$cont_cov,probs=seq(0,1,0.05),na.rm=T)
quantile(f1dat3$sen_cov,probs=seq(0,1,0.1),na.rm=T)

#take everything with abs(CoV) < 5
f1dat6 <- f1dat3[abs(f1dat3$cont_cov) < 5,]
f1dat7 <- f1dat6[abs(f1dat6$sen_cov) < 5,]
f1dat7 <- f1dat7[!is.na(f1dat7$cont_mean),]

keep <- rownames(f1dat7)

f1dat.new <- f1_TPM_org[!f1filt1,]

f1dat.final <- f1dat.new[f1dat.new$Genes %in% keep,]

f1dat.final <- f1dat.final[!is.na(f1dat.final$PWK.CONT.1),]
f1dat.final <- f1dat.final[,-1]

#now run ANOVA
final <- data.frame(Gene=rownames(f1dat.final),anova_F=c(rep("",nrow(f1dat.final))),anova_pval=c(rep("",nrow(f1dat.final))))
rownames(final) <- final$Gene

genelist <- final$Gene

for (i in genelist){
  dat <- f1dat.final[i,]
  
  anova_df <- data.frame(TPM=t(dat))
  colnames(anova_df) <- "TPM"
  anova_df$Allele <- sapply(strsplit(rownames(anova_df),".",fixed=T),"[",1)
  anova_df$Condition <- sapply(strsplit(rownames(anova_df),".",fixed=T),"[",2)
  
  aout <- aov(TPM ~ Allele * Condition,data=anova_df)
  aout_s <- summary(aout)
  a_pval <- aout_s[[1]][["Pr(>F)"]][3]
  a_F <- aout_s[[1]][["F value"]][3]
  
  final[i,2] <- a_F
  final[i,3] <- a_pval
}

final <- final[order(final$anova_pval),]
final$BH <- p.adjust(final$anova_pval,method="BH")

write.table(final[,c(1,2)],file="/Volumes/tkang/primary_ASE/CoV_filt/allF1_anovaF.txt",sep="\t",eol="\n",quote=F,col.names=T,row.names=F)

anova_F <- sapply(final$anova_F,as.numeric)
vq <- quantile(anova_F)
top25 <- vq[[4]]
ase <- rownames(final[as.numeric(final$anova_F) > top25,])
notase <- rownames(final[as.numeric(final$anova_F) < top25,])

write.table(ase,file="/Volumes/tkang/primary_ASE/CoV_filt/top25_anova_ASE.txt",sep="\t",eol="\n",quote=F,col.names=F,row.names=F)
write.table(notase,file="/Volumes/tkang/primary_ASE/CoV_filt/top25_anova_notASE.txt",sep="\t",eol="\n",quote=F,col.names=F,row.names=F)

##################################################################################################################################
##################################################################################################################################
#analyze the table after allnested_call_var_gtrd_sanger_bychrom_clean.py
#do CoV cutoff
anova_NA <- read.table("/Volumes/tkang/primary_ASE/CoV_filt/NA_top25anova_DESeq_norm_GTRDallTF_chisq_forR.txt",header=F,quote="",sep="\t")
anova_NA1 <- anova_NA[anova_NA$V2 > 250,]
anova_NA2 <- anova_NA1[anova_NA1$V3 > 250,]
anova_NA3 <- anova_NA2[anova_NA2$V4 > 250,]
anova_NA4 <- anova_NA3[anova_NA3$V5 > 250,]

#fisher's
anova_out <- data.frame(Genes=c(anova_NA4[,1]),fisher_pval=c(rep("",nrow(anova_NA4))))
rownames(anova_out) <- anova_out$Genes

for (i in 1:nrow(anova_NA4)){
  tempdf <- data.frame(ASE=c(anova_NA4[i,2],anova_NA4[i,4]),noASE=c(anova_NA4[i,3],anova_NA4[i,5]))
  rownames(tempdf) <- c("var","novar")
  gene <- anova_NA4[i,1]
  
  fisher_out <- fisher.test(tempdf,alternative="greater")
  #or <- as.numeric(chisq_out[[1]])
  pval <- fisher_out$p.value
  #chistat <- chisq_out[[1]]
  anova_out[gene,2] <- pval
  #anova_out[gene,3] <- chistat
  #anova_out[gene,3] <- or
}

anova_out$BH <- p.adjust(anova_out$fisher_pval,method="BH")
anova_out <- anova_out[order(anova_out$fisher_pval),]

