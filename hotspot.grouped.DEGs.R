library(ggplot2)
library(data.table)
library(dplyr)
library(grid)
library(ggrepel)
library(biomaRt)
library('org.Hs.eg.db')
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(limma)

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

probe=fread("~/gencode.v22.annotation.gene.probeMap")
express.ucec=data.frame(fread("E:/my_study/Project/TCGA.DATA/cancertype/TCGA-UCEC.htseq_counts.tsv.gz"))
express.ucec=left_join(express.ucec,probe[,c(1,2)],by=c("Ensembl_ID"="id"))
express.ucec=express.ucec[,c(1,ncol(express.ucec),2:(ncol(express.ucec)-1))]
express.ucec=na.omit(express.ucec)
express.ucec=express.ucec[!duplicated(express.ucec$gene),]
rownames(express.ucec)=express.ucec$gene
express.ucec=express.ucec[,-c(1,2)]
colnames(express.ucec)=substr(colnames(express.ucec),1,15)
express.ucec=express.ucec[,grepl(".01$",colnames(express.ucec))]
colnames(express.ucec)=gsub("\\.","-",colnames(express.ucec))
group=rep("UCEC",ncol(express.ucec))
express.ucec.cpm=2^express.ucec-1
express.ucec.cpm <- normalizeQuantiles(express.ucec.cpm)
y <- DGEList(express.ucec.cpm, group=group, genes=rownames(express.ucec))
keep <- filterByExpr(y)
y <- y[keep,, keep.lib.sizes=FALSE]
y=log2(y$counts+1)
save.y=data.frame(y)

# saveRDS(save.y,file = "E:/my_study/Project/TCGA.DATA/UCEC.rds")
# express.tcga=readRDS("e:/my_study/Project/TCGA.DATA/pancancer/tcga.express.filter.1.1.rds")
# express.ucec=readRDS("E:/my_study/Project/TCGA.DATA/UCEC.rds")
colnames(express.ucec)=gsub("\\.","-",colnames(express.ucec))


maf=fread("E:/my_study/Project/TCGA.DATA/MAF/mc3.v0.2.8.PUBLIC.maf")
maf.keep=readRDS(file = "E:/my_study/Project/cancer_risk/code2020.9.22/maf.keep.rds")
maf.keep$join=paste(maf.keep$gene,maf.keep$patient,maf.keep$HGVSc,sep = ".")
maf2paint=maf[,c("Hugo_Symbol","Tumor_Sample_Barcode","Transcript_ID","Chromosome","Start_Position","HGVSc","HGVSp","Variant_Classification")]
maf2paint$join=paste(maf2paint$Hugo_Symbol,substr(maf2paint$Tumor_Sample_Barcode,1,12),maf2paint$HGVSc,sep = ".")
maf.use=left_join(maf.keep,maf2paint,by="join")
maf.use=maf.use[!duplicated(maf.use$join),]
maf.use=maf.use[,c(8:15)]
colnames(maf.use)=c("Hugo_Symbol","Tumor_Sample_Barcode","Transcript_ID","Chromosome","Start_Position","HGVSc","HGVSp","Variant_Classification")
maf.use$submitter_id=substr(maf.use$Tumor_Sample_Barcode,1,12)
clinical.biolinks=readRDS("E:/my_study/Project/TCGA.DATA/clinical/clinical_tcgabiolinks.rds")
clinical.joint=clinical.biolinks[,c(1,5,6,7,25,28,31,33,42)]
clinical.joint$submitter_id=paste(clinical.joint$submitter_id,"-01",sep = "")
maf.use=left_join(maf.use,clinical.joint[,c(1,7,8,9)],by="submitter_id")
maf.use$Tumor_Sample_Barcode=substr(maf.use$Tumor_Sample_Barcode,1,12)
maf.use=maf.use[!duplicated(maf.use),]
maf.stad=maf.use[maf.use$disease=="STAD",]  
maf.coad=maf.use[maf.use$disease=="COAD",]  
maf.ucec=maf.use[maf.use$disease=="UCEC",]
maf.skcm=maf.use[maf.use$disease=="SKCM",]

maf.table=maf[,c(1,16)]
maf.table$Tumor_Sample_Barcode=substr(maf.table$Tumor_Sample_Barcode,1,12)
maf.table=left_join(maf.table,clinical.joint[,c(1,7,8,9)],by=c("Tumor_Sample_Barcode"="submitter_id"))
maf.table=maf.table[!duplicated(maf.table$Tumor_Sample_Barcode),]

####stad####
stad.mutation.table=data.frame(patient=maf.table$Tumor_Sample_Barcode[maf.table$disease=="STAD"])
stad.mutation.table=na.omit(stad.mutation.table)
stad.mutation.table$PTEN=ifelse(stad.mutation.table$patient %in% maf.stad$Tumor_Sample_Barcode[maf.stad$Hugo_Symbol=="PTEN"],1,0)
stad.mutation.table$PTEN.hot=ifelse(stad.mutation.table$patient %in% maf.stad$Tumor_Sample_Barcode[maf.stad$Hugo_Symbol=="PTEN" & maf.stad$HGVSp=="p.Lys267ArgfsTer9"],1,0)
stad.mutation.table$PGM5=ifelse(stad.mutation.table$patient %in% maf.stad$Tumor_Sample_Barcode[maf.stad$Hugo_Symbol=="PGM5"],1,0)
stad.mutation.table$PGM5.hot=ifelse(stad.mutation.table$patient %in% maf.stad$Tumor_Sample_Barcode[maf.stad$Hugo_Symbol=="PGM5" & maf.stad$HGVSp=="p.Ile98Val"],1,0)
stad.mutation.table$LARP4B=ifelse(stad.mutation.table$patient %in% maf.stad$Tumor_Sample_Barcode[maf.stad$Hugo_Symbol=="LARP4B"],1,0)
stad.mutation.table$LARP4B.hot=ifelse(stad.mutation.table$patient %in% maf.stad$Tumor_Sample_Barcode[maf.stad$Hugo_Symbol=="LARP4B" & maf.stad$HGVSp=="p.Thr163HisfsTer47"],1,0)
stad.mutation.table$HOTGENE=ifelse(stad.mutation.table$PTEN==1 | stad.mutation.table$PGM5==1 | stad.mutation.table$LARP4B==1,1,0)
stad.mutation.table$HOTPOT=ifelse(stad.mutation.table$PTEN.hot==1 | stad.mutation.table$PGM5.hot==1 | stad.mutation.table$LARP4B.hot==1,1,0)
stad.mutation.table=left_join(stad.mutation.table,clinical.joint[,c(1,7)],by=c("patient"="submitter_id"))
stad.mutation.table$patient=paste(stad.mutation.table$patient,"-01",sep = "")
stad.mutation.table=stad.mutation.table[stad.mutation.table$patient %in% colnames(express.tcga),]
####UCEC####
ucec.mutation.table=data.frame(patient=maf.table$Tumor_Sample_Barcode[maf.table$disease=="UCEC"])
ucec.mutation.table=na.omit(ucec.mutation.table)
ucec.mutation.table$LARP4B=ifelse(ucec.mutation.table$patient %in% maf.ucec$Tumor_Sample_Barcode[maf.ucec$Hugo_Symbol=="LARP4B"],1,0)
ucec.mutation.table$LARP4B.hot=ifelse(ucec.mutation.table$patient %in% maf.ucec$Tumor_Sample_Barcode[maf.ucec$Hugo_Symbol=="LARP4B" & maf.ucec$HGVSp=="p.Thr163HisfsTer47"],1,0)
ucec.mutation.table$PTEN=ifelse(ucec.mutation.table$patient %in% maf.ucec$Tumor_Sample_Barcode[maf.ucec$Hugo_Symbol=="PTEN"],1,0)
ucec.mutation.table$PTEN.hot=ifelse(ucec.mutation.table$patient %in% maf.ucec$Tumor_Sample_Barcode[maf.ucec$Hugo_Symbol=="PTEN" & maf.ucec$HGVSp=="p.Lys267ArgfsTer9"],1,0)
ucec.mutation.table$PGM5=ifelse(ucec.mutation.table$patient %in% maf.ucec$Tumor_Sample_Barcode[maf.ucec$Hugo_Symbol=="PGM5"],1,0)
ucec.mutation.table$PGM5.hot=ifelse(ucec.mutation.table$patient %in% maf.ucec$Tumor_Sample_Barcode[maf.ucec$Hugo_Symbol=="PGM5" & maf.ucec$HGVSp=="p.Ile98Val"],1,0)
ucec.mutation.table$HOTGENE=ifelse((ucec.mutation.table$LARP4B==1 | ucec.mutation.table$PTEN==1 | ucec.mutation.table$PGM5==1), 1, 0)
ucec.mutation.table$HOTPOT=ifelse((ucec.mutation.table$LARP4B.hot==1 | ucec.mutation.table$PTEN.hot==1 | ucec.mutation.table$PGM5.hot==1), 1, 0)
ucec.mutation.table=left_join(ucec.mutation.table,clinical.joint[,c(1,7)],by=c("patient"="submitter_id"))
ucec.mutation.table$patient=paste(ucec.mutation.table$patient,"-01",sep = "")
ucec.mutation.table=ucec.mutation.table[!duplicated(ucec.mutation.table),]
ucec.mutation.table=ucec.mutation.table[ucec.mutation.table$patient %in% colnames(express.ucec),]

####COAD####
coad.mutation.table=data.frame(patient=maf.table$Tumor_Sample_Barcode[maf.table$disease=="COAD"])
coad.mutation.table$PGM5=ifelse(coad.mutation.table$patient %in% maf.coad$Tumor_Sample_Barcode[maf.coad$Hugo_Symbol=="PGM5"],1,0)
coad.mutation.table$PGM5.hot=ifelse(coad.mutation.table$patient %in% maf.coad$Tumor_Sample_Barcode[maf.coad$Hugo_Symbol=="PGM5" & maf.coad$HGVSp=="p.Ile98Val"],1,0)
coad.mutation.table$PTEN=ifelse(coad.mutation.table$patient %in% maf.coad$Tumor_Sample_Barcode[maf.coad$Hugo_Symbol=="PTEN"],1,0)
coad.mutation.table$PTEN.hot=ifelse(coad.mutation.table$patient %in% maf.coad$Tumor_Sample_Barcode[maf.coad$Hugo_Symbol=="PTEN" & maf.coad$HGVSp=="p.Lys267ArgfsTer9"],1,0)
coad.mutation.table$LARP4B=ifelse(coad.mutation.table$patient %in% maf.coad$Tumor_Sample_Barcode[maf.coad$Hugo_Symbol=="LARP4B"],1,0)
coad.mutation.table$LARP4B.hot=ifelse(coad.mutation.table$patient %in% maf.coad$Tumor_Sample_Barcode[maf.coad$Hugo_Symbol=="LARP4B" & maf.coad$HGVSp=="p.Thr163HisfsTer47"],1,0)
coad.mutation.table$HOTGENE=ifelse(coad.mutation.table$LARP4B==1 | coad.mutation.table$PTEN==1 | coad.mutation.table$PGM5==1, 1, 0)
coad.mutation.table$HOTPOT=ifelse(coad.mutation.table$LARP4B.hot==1 | coad.mutation.table$PTEN.hot==1 | coad.mutation.table$PGM5.hot==1, 1, 0)
coad.mutation.table=left_join(coad.mutation.table,clinical.joint[,c(1,7)],by=c("patient"="submitter_id"))
coad.mutation.table$patient=paste(coad.mutation.table$patient,"-01",sep = "")
coad.mutation.table=coad.mutation.table[coad.mutation.table$patient %in% colnames(express.tcga),]

mutation.table.list=list(STAD=stad.mutation.table,COAD=coad.mutation.table,UCEC=ucec.mutation.table)

cancer=c("STAD","COAD","UCEC")

gene.hotpot.pvalue=list()
for(j in 1:length(cancer)){
  if(cancer[j]!="UCEC"){
    cancer.sample=clinical.joint$submitter_id[clinical.joint$disease==cancer[j]]
    cancer.express=express.tcga[colnames(express.tcga) %in% cancer.sample]
  }else{
    cancer.express=express.ucec
  }
  design.tmp=data.frame(sample=colnames(cancer.express))
  design.tmp=left_join(design.tmp,mutation.table.list[[cancer[j]]][,c(1,9,10)],by=c("sample"="patient"))
  colnames(design.tmp)[3]="gender"
  colnames(design.tmp)[2]="hotpot"
  design.tmp=na.omit(design.tmp)
  cancer.express=cancer.express[,design.tmp$sample]
  gender=design.tmp$gender
  hotpot=design.tmp$hotpot
  if(j!=3){
    design <- model.matrix(~gender*hotpot)
  }else{
    design <- model.matrix(~hotpot)
  }
  vfit <- lmFit(cancer.express,design)
  efit <- eBayes(vfit,robust=TRUE,trend=TRUE)
  if(j!=3){
    DEG.hotpot = topTable(efit, coef=3, n=Inf)
  }else{
    DEG.hotpot = topTable(efit, coef=2, n=Inf)
  }
  gene.hotpot.pvalue[[cancer[j]]]=DEG.hotpot
}


library(openxlsx)
write.xlsx(gene.hotpot.pvalue, "~/hotspot.grouped.DEGs.xlsx",rowNames = TRUE) 
