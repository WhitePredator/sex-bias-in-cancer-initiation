library(data.table)
library(dplyr)
library(purrr)
library(readr)
library(ggplot2)
library(forcats)
library(rstatix)
library(limma)
library(grid)
library(ggrepel)
library(biomaRt)
library(ggpubr)
library('org.Hs.eg.db')
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
memory.limit(50000000)

clinical.tcga=readRDS("~/clinical.joint.rds")

express.tcga=fread("~/tcga_RSEM_Hugo_norm_count.gz")
express.tcga=data.frame(express.tcga)
express.tcga=express.tcga[!duplicated(express.tcga$sample),]
rownames(express.tcga)=express.tcga$sample
express.tcga=express.tcga[,-1]
colnames(express.tcga)=gsub("\\.","-",colnames(express.tcga))
express.tcga=na.omit(express.tcga)
express.tcga.cpm=2^express.tcga-1
express.tcga.cpm= t(t(express.tcga.cpm)/colSums(express.tcga.cpm) * 1000000)
keep=rowSums(express.tcga.cpm >= 1 ) >= 1
express.tcga=express.tcga[keep,]
express.tcga=express.tcga[,colnames(express.tcga) %in% c(paste(clinical.joint$patient,"-01",sep = ""),melanoma.sample.06)]
sum.row=data.frame(colSums(express.tcga))
patient.outlier=rownames(sum.row)[sum.row$"colSums.express.tcga."<10000]
express.tcga=express.tcga[,!colnames(express.tcga) %in% patient.outlier]

msi.3=fread("~/TCGASubtype.20170308.tsv.gz")
msi.3=left_join(msi.3,clinical.tcga[,c(9,6,8)],by=c("sampleID"="sample"))
msi.3$patient=substr(msi.3$sampleID,1,12)

cancer.focus=c("STAD","COAD","UCEC","LGG","GBM","ACC","THCA","SKCM","PCPG","HNSC","LUAD","PRAD","BRCA")

STAD.sample=msi.3[msi.3$disease=="STAD",]
STAD.sample=STAD.sample[,c(1,8,10,11,12)]
STAD.sample=na.omit(STAD.sample)
STAD.sample$group=ifelse(STAD.sample$Subtype_other %in% c("HM-indel","HM-SNV"),"MSI+","MSI-")

COAD.sample=msi.3[msi.3$disease=="COAD",]
COAD.sample=COAD.sample[,c(1,8,10,11,12)]
COAD.sample=na.omit(COAD.sample)
COAD.sample$group=ifelse(COAD.sample$Subtype_other %in% c("HM-indel","HM-SNV"),"MSI+","MSI-")

UCEC.sample=msi.3[msi.3$disease=="UCEC",]
UCEC.sample=UCEC.sample[UCEC.sample$Subtype_Selected %in% c("UCEC.CN_LOW","UCEC.CN_HIGH","UCEC.POLE","UCEC.MSI"),]
UCEC.sample=UCEC.sample[,c(1,9,10,11,12)]
UCEC.sample=na.omit(UCEC.sample)
UCEC.sample$group=ifelse(UCEC.sample$Subtype_Selected %in% c("UCEC.POLE","UCEC.MSI"),"MSI+","MSI-")


LGG.sample=msi.3[msi.3$disease=="LGG",]
LGG.sample=LGG.sample[LGG.sample$Subtype_other %in% c("G-CIMP-high","G-CIMP-low"),]
LGG.sample=LGG.sample[,c(1,8,10,11,12)]
LGG.sample=na.omit(LGG.sample)
LGG.sample$group=ifelse(LGG.sample$Subtype_other=="G-CIMP-high","CIMP+","CIMP-")

GBM.sample=msi.3[msi.3$disease=="GBM",]
GBM.sample=GBM.sample[GBM.sample$Subtype_other %in% c("G-CIMP-high","G-CIMP-low"),]
GBM.sample=GBM.sample[,c(1,8,10,11,12)]
GBM.sample=na.omit(GBM.sample)
GBM.sample$group=ifelse(GBM.sample$Subtype_other=="G-CIMP-high","CIMP+","CIMP-")

ACC.sample=msi.3[msi.3$disease=="ACC",]
ACC.sample=ACC.sample[ACC.sample$Subtype_Selected %in% c("ACC.CIMP-high","ACC.CIMP-low"),]
ACC.sample=ACC.sample[,c(1,9,10,11,12)]
ACC.sample=na.omit(ACC.sample)
ACC.sample$group=ifelse(ACC.sample$Subtype_Selected=="ACC.CIMP-high","CIMP+","CIMP-")

THCA.sample=msi.3[msi.3$disease=="THCA",]
THCA.sample=THCA.sample[THCA.sample$Subtype_DNAmeth %in% c("follicular","CpG island methylated"),]
THCA.sample=THCA.sample[,c(1,3,10,11,12)]
THCA.sample=na.omit(THCA.sample)
THCA.sample$group=ifelse(THCA.sample$Subtype_DNAmeth=="CpG island methylated","CIMP+","Follicular")

SKCM.sample=msi.3[msi.3$disease=="SKCM",]
SKCM.sample=SKCM.sample[SKCM.sample$Subtype_DNAmeth %in% c("CpG island-methylated","normal-like","hyper-methylated","hypo-methylated"),]
SKCM.sample=SKCM.sample[,c(1,3,10,11,12)]
SKCM.sample=na.omit(SKCM.sample)
SKCM.sample$group=ifelse(SKCM.sample$Subtype_DNAmeth %in% c("CpG island-methylated","hyper-methylated"),"Hyper-meth","Hypo-meth")

PCPG.sample=msi.3[msi.3$disease=="PCPG",]
PCPG.sample=PCPG.sample[PCPG.sample$Subtype_DNAmeth %in% c("low-methylated","hyper-methylated"),]
PCPG.sample=PCPG.sample[,c(1,3,10,11,12)]
PCPG.sample=na.omit(PCPG.sample)
PCPG.sample$group=ifelse(PCPG.sample$Subtype_DNAmeth=="hyper-methylated","Hyper-meth","Hypo-meth")

HNSC.sample=msi.3[msi.3$disease=="HNSC",]
HNSC.sample=HNSC.sample[HNSC.sample$Subtype_DNAmeth %in% c("hyper-methylated","hypo-methylated","normal-like","CpG island hyper-methylated"),]
HNSC.sample=HNSC.sample[,c(1,3,10,11,12)]
HNSC.sample=na.omit(HNSC.sample)
HNSC.sample$group=ifelse(HNSC.sample$Subtype_DNAmeth %in% c("CpG island hyper-methylated","hyper-methylated"),"Hyper-meth","Hypo-meth")

LUAD.sample=msi.3[msi.3$disease=="LUAD",]
LUAD.sample=LUAD.sample[LUAD.sample$Subtype_DNAmeth %in% c("high","low"),]
LUAD.sample=LUAD.sample[,c(1,3,10,11,12)]
LUAD.sample=na.omit(LUAD.sample)
LUAD.sample$group=ifelse(LUAD.sample$Subtype_DNAmeth=="high","Hyper-meth","Hypo-meth")

PRAD.sample=msi.3[msi.3$disease=="PRAD",]
PRAD.sample=PRAD.sample[PRAD.sample$Subtype_DNAmeth %in% c("1","4"),]
PRAD.sample=PRAD.sample[,c(1,3,10,11,12)]
PRAD.sample=na.omit(PRAD.sample)
PRAD.sample$group=ifelse(PRAD.sample$Subtype_DNAmeth =="1","Hyper-meth","Hypo-meth")

BRCA.sample=msi.3[msi.3$disease=="BRCA",]
BRCA.sample=BRCA.sample[BRCA.sample$Subtype_mRNA %in% c("LumA","LumB","Her2","Basal"),]
BRCA.sample=BRCA.sample[,c(1,2,10,11,12)]
BRCA.sample=na.omit(BRCA.sample)
BRCA.sample$group=ifelse(BRCA.sample$Subtype_mRNA %in% c("LumA","LumB"),"ER+","ER-")
BRCA.sample=BRCA.sample[BRCA.sample$gender=="female",]

group.list=list(STAD.sample,COAD.sample,UCEC.sample,LGG.sample,GBM.sample,ACC.sample,THCA.sample,SKCM.sample,PCPG.sample,
                HNSC.sample,LUAD.sample,PRAD.sample,BRCA.sample)
names(group.list)=cancer.focus

for(i in 1:length(group.list)){
  if(i<4){
    group.list[[i]]$group=relevel(factor(group.list[[i]]$group),ref= "MSI-")
  }else if(i<7){
    group.list[[i]]$group=relevel(factor(group.list[[i]]$group),ref= "CIMP-")
  }else if(i==7){
    group.list[[i]]$group=relevel(factor(group.list[[i]]$group),ref= "Follicular")
  }else if(i<13){
    group.list[[i]]$group=relevel(factor(group.list[[i]]$group),ref= "Hypo-meth")
  }else{
    group.list[[i]]$group=relevel(factor(group.list[[i]]$group),ref= "ER-")
  }
}


DEGs.pvalue=list()
for(j in 1:length(cancer.focus)){
  cancer.sample=group.list[[j]]$sampleID
  cancer.express=express.tcga[colnames(express.tcga) %in% cancer.sample]
  cancer.express=na.omit(cancer.express)
  thr.express=cancer.express["THRA",]+cancer.express["THRB",]
  hif.express=cancer.express["HIF1A",]-cancer.express["HIF3A",]
  combine.gene=rbind(thr.express,hif.express)
  rownames(combine.gene)=c("THR","HIF")
  cancer.express=rbind(cancer.express,combine.gene)
  design.tmp=data.frame(sample=colnames(cancer.express))
  design.tmp=left_join(design.tmp,group.list[[j]][,c(1,6)],by=c("sample"="sampleID"))
  design.tmp=left_join(design.tmp,clinical.tcga,by="sample")
  design.tmp=design.tmp[,-3]
  design.tmp=na.omit(design.tmp)
  cancer.express=cancer.express[,design.tmp$sample]
  for(k in ncol(design.tmp):2){
    if(length(unique(design.tmp[,k]))==1){
      design.tmp=design.tmp[,-k]
    }
  }
  form.limma <- as.formula(paste("~ ",paste(colnames(design.tmp)[-1],collapse="+"),sep=""))
  design <- model.matrix(form.limma,design.tmp)
  vfit <- lmFit(cancer.express,design)
  efit <- eBayes(vfit,robust=TRUE,trend=TRUE)
  DEG.cimp = topTable(efit, coef=2, n=Inf)
  DEG.cimp$gene=rownames(DEG.cimp)
  DEGs.pvalue[[cancer.focus[j]]]=DEG.cimp
}



library(openxlsx)
write.xlsx(DEGs.pvalue, "~/subtype.grouped.DEGs.xlsx",rowNames = TRUE) 

