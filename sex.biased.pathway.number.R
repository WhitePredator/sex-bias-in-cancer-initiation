
library(dplyr)
library(data.table)
library(purrr)
library(KEGGREST)
library(org.Hs.eg.db)
library(tidyverse)
library(limma)

express.gtex=fread("~/gtex_RSEM_Hugo_norm_count.gz")
express.gtex=data.frame(express.gtex)
express.gtex=express.gtex[!duplicated(express.gtex$sample),]
rownames(express.gtex)=express.gtex$sample
express.gtex=express.gtex[,-1]
colnames(express.gtex)=gsub("\\.","-",colnames(express.gtex))
express.gtex=express.gtex[,colnames(express.gtex) %in% clinic.gtex$Sample]
express.gtex=na.omit(express.gtex)
express.gtex.cpm=2^express.gtex-1
express.gtex.cpm= t(t(express.gtex.cpm)/colSums(express.gtex.cpm) * 1000000)
keep=rowSums(express.gtex.cpm >= 1 ) >= 1
express.gtex=express.gtex[keep,]
sum.row=data.frame(colSums(express.gtex))
patient.outlier=rownames(sum.row)[sum.row$"colSums.express.gtex."<10000]
express.gtex=express.gtex[,!colnames(express.gtex) %in% patient.outlier]

clinic.gtex=fread("~/gtex/GTEX_phenotype.gz")
clinic.gtex=clinic.gtex[!grepl("Cells",clinic.gtex$`body_site_detail (SMTSD)`),]
clinic.gtex=clinic.gtex[clinic.gtex$Sample %in% colnames(express.gtex),]

hsa_path_eg  <- keggLink("pathway", "hsa") %>% 
  tibble(pathway = ., eg = sub("hsa:", "", names(.)))
pathnames=data.frame(pathway=names(keggList("pathway", "hsa")),name=keggList("pathway", "hsa"))
hsa_path_eg=left_join(hsa_path_eg,pathnames,by="pathway")
hsa_kegg_anno <- hsa_path_eg %>%
  mutate(
    symbol = mapIds(org.Hs.eg.db, eg, "SYMBOL", "ENTREZID"),
    ensembl = mapIds(org.Hs.eg.db, eg, "ENSEMBL", "ENTREZID")
  )
pathway.hsa=unique(hsa_kegg_anno$pathway)

cancer.type.gtex=c("Bladder","Kidney","Liver","Skin","Colon","Stomach","Pancreas","Lung","Thyroid","Brain","Esophagus","Adrenal Gland")

clinic.gtex=clinic.gtex[clinic.gtex$"_primary_site" %in% cancer.type.gtex,]
clinic.gtex.list=split(clinic.gtex,clinic.gtex$"_primary_site")
express.gtex=express.gtex[,colnames(express.gtex) %in% clinic.gtex$Sample]

group=data.frame(Sample=colnames(express.gtex))
group=left_join(group,clinic.gtex[,c(1,3)],by="Sample")
express.gtex=t(express.gtex)
express.gtex=data.frame(express.gtex)
express.gtex$group=group$"_primary_site"
express.list=split(express.gtex,express.gtex$group)

for(i in 1:12){
  express.list[[i]]=data.frame(t(express.list[[i]]))
  rowname=rownames(express.list[[i]])
  express.list[[i]]=data.frame(apply(express.list[[i]],2,as.numeric))
  rownames(express.list[[i]])=rowname
  colnames(express.list[[i]])=gsub("\\.","-",colnames(express.list[[i]]))
}
sample.size=100
repeat.time=1000
sig.repeat.matrix=data.frame(matrix(ncol = 12,nrow = repeat.time))
colnames(sig.repeat.matrix)=names(clinic.gtex.list)
set.seed(31415926)

kegg.func=function(kegg.term,hsa_kegg_anno,cancer.exp){
  kegg.gene=hsa_kegg_anno$symbol[hsa_kegg_anno$pathway==kegg.term]
  kegg.gene.express=cancer.exp[kegg.gene,]
  kegg.gene.express=na.omit(kegg.gene.express)
  kegg.score=apply(kegg.gene.express,2,mean)
  kegg.score=data.frame(kegg.score)
  colnames(kegg.score)=kegg.term
  return(kegg.score)
}

for(i in 1:12){
  print(i)
  male.sample=clinic.gtex.list[[i]]$Sample[clinic.gtex.list[[i]]$"_gender"=="male"]
  female.sample=clinic.gtex.list[[i]]$Sample[clinic.gtex.list[[i]]$"_gender"=="female"]
  male.size=ceiling(length(male.sample)*sample.size/nrow(clinic.gtex.list[[i]]))
  female.size=ceiling(length(female.sample)*sample.size/nrow(clinic.gtex.list[[i]]))
  if(nrow(clinic.gtex.list[[i]])>sample.size){
    for(j in 1:repeat.time){
      male.tmp=sample(length(male.sample),male.size)
      female.tmp=sample(length(female.sample),female.size)
      express.tmp=express.list[[i]][,c(male.sample[male.tmp],female.sample[female.tmp])]
      cancer.kegg=sapply(pathway.hsa, kegg.func,hsa_kegg_anno, express.tmp, simplify = FALSE)
      cancer.kegg=purrr::reduce(cancer.kegg,bind_cols)
      pathway.score=t(cancer.kegg)
      design.tmp=data.frame(Sample=colnames(pathway.score))
      design.tmp=left_join(design.tmp,clinic.gtex[,c(1,2,4)],by="Sample")
      colnames(design.tmp)[3]="gender"
      colnames(design.tmp)[2]="primary"
      gender=design.tmp$gender
      primary=design.tmp$primary
      if(length(unique(primary))!=1){
        design <- model.matrix(~gender+primary)
      }else{
        design <- model.matrix(~gender)
      }
      score.limma <- voom(pathway.score, design)
      vfit <- lmFit(pathway.score,design)
      efit <- eBayes(vfit,robust=TRUE,trend=TRUE)
      tempOutput.mutation = topTable(efit, coef=2, n=Inf,sort.by = "none")
      sig.repeat.matrix[j,i]=sum(tempOutput.mutation$P.Value<0.05)
    }
  }
  
  else{
    express.tmp=express.list[[i]]
    cancer.kegg=sapply(pathway.hsa, kegg.func,hsa_kegg_anno, express.tmp, simplify = FALSE)
    cancer.kegg=purrr::reduce(cancer.kegg,bind_cols)
    pathway.score=t(cancer.kegg)
    design.tmp=data.frame(Sample=colnames(pathway.score))
    design.tmp=left_join(design.tmp,clinic.gtex.list[[i]][,c(1,2,4)],by="Sample")
    colnames(design.tmp)[3]="gender"
    colnames(design.tmp)[2]="primary"
    gender=design.tmp$gender
    primary=design.tmp$primary
    if(length(unique(primary))!=1){
      design <- model.matrix(~gender+primary)
    }else{
      design <- model.matrix(~gender)
    }
    score.limma <- voom(pathway.score, design)
    vfit <- lmFit(pathway.score,design)
    efit <- eBayes(vfit,robust=TRUE,trend=TRUE)
    tempOutput.mutation = topTable(efit, coef=2, n=Inf,sort.by = "none")
    sig.repeat.matrix[,i]=sum(tempOutput.mutation$P.Value<0.05)
  }
}

library(forcats)
sig.repeat.melt=melt(sig.repeat.matrix,value.name = "number",variable.name = "tissue")
sig.repeat.melt %>%
  mutate(tissue = fct_reorder(tissue, dplyr::desc(number))) %>%
  ggplot( aes(x = tissue,y=number,fill = tissue))+
  geom_boxplot(outlier.shape=NA) +
  theme_classic()+
  theme(legend.position="none")+
  scale_y_continuous(breaks=seq(-10,170,10), expand = c(0, 0.3),limits = c(-10,170)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


