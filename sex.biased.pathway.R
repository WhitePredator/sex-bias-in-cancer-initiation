
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(limma)
library(ggpubr)
library(rstatix)
library(purrr)
library(grid)


clinic.gtex=fread("e:/my_study/Project/GTEX/GTEX_phenotype.gz")
clinic.gtex=clinic.gtex[!grepl("Cells",clinic.gtex$`body_site_detail (SMTSD)`),]


express.gtex=fread("E:/my_study/Project/GTEX/gtex_RSEM_Hugo_norm_count.gz")
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

clinic.gtex=clinic.gtex[clinic.gtex$Sample %in% colnames(express.gtex),]
colnames(clinic.gtex)[1]="sample"

library(KEGGREST)
library(org.Hs.eg.db)
library(tidyverse)

hsa_path_eg  <- keggLink("pathway", "hsa") %>% 
  tibble(pathway = ., eg = sub("hsa:", "", names(.)))
pathnames=data.frame(pathway=names(keggList("pathway", "hsa")),name=keggList("pathway", "hsa"))
hsa_path_eg$pathway=substr(hsa_path_eg$pathway,6,100)
hsa_path_eg=left_join(hsa_path_eg,pathnames,by="pathway")
hsa_kegg_anno <- hsa_path_eg %>%
  mutate(
    symbol = mapIds(org.Hs.eg.db, eg, "SYMBOL", "ENTREZID"),
    ensembl = mapIds(org.Hs.eg.db, eg, "ENSEMBL", "ENTREZID")
  )


pathway.hsa=unique(hsa_kegg_anno$pathway)

kegg.func=function(kegg.term,hsa_kegg_anno,cancer.exp){
  kegg.gene=hsa_kegg_anno$symbol[hsa_kegg_anno$pathway==kegg.term]
  kegg.gene.express=cancer.exp[kegg.gene,]
  kegg.gene.express=na.omit(kegg.gene.express)
  kegg.score=apply(kegg.gene.express,2,mean)
  kegg.score=data.frame(kegg.score)
  colnames(kegg.score)=kegg.term
  return(kegg.score)
}


cancer.func=function(cancer.type, clinic, express,pathway.hsa,TCGA){
  if(TCGA){
    cancer.p=clinic$sample[clinic$disease==cancer.type]
  }else{
    cancer.p=clinic$sample[clinic$"_primary_site"==cancer.type]
  }
  cancer.exp=express[,colnames(express) %in% cancer.p]
  cancer.kegg=sapply(pathway.hsa, kegg.func,hsa_kegg_anno, cancer.exp, simplify = FALSE)
  cancer.kegg=purrr::reduce(cancer.kegg,bind_cols)
  return(cancer.kegg)
}



cancer.type.gtex=c("Kidney","Liver","Skin","Colon","Stomach","Pancreas","Lung","Thyroid","Brain","Esophagus","Adrenal Gland")

cancer.kegg.gtex=sapply(cancer.type.gtex, cancer.func, clinic.gtex, express.gtex,pathway.hsa,0,simplify = FALSE)




####GTEX####
pathway.test.gtex=list()
for(i in 1:length(cancer.kegg.gtex)){
  pathway.score=t(cancer.kegg.gtex[[i]])
  design.tmp=data.frame(sample=colnames(pathway.score))
  design.tmp=left_join(design.tmp,clinic.gtex[,c(1,2,4)],by="sample")
  colnames(design.tmp)[3]="gender"
  colnames(design.tmp)[2]="primary"
  gender=design.tmp$gender
  primary=design.tmp$primary
  if(length(unique(primary))!=1){
    design <- model.matrix(~gender+primary)
  }else{
    design <- model.matrix(~gender)
  }
  vfit <- lmFit(pathway.score,design)
  efit <- eBayes(vfit,robust=TRUE,trend=TRUE)
  tempOutput.mutation = topTable(efit, coef=2, n=Inf,sort.by = "none")
  tempOutput.mutation$pathway=rownames(tempOutput.mutation)
  tempOutput.mutation=left_join(tempOutput.mutation,hsa_kegg_anno[,c(1,3)],by="pathway")
  tempOutput.mutation=tempOutput.mutation[!duplicated(tempOutput.mutation),]
  pathway.test.gtex[[names(cancer.kegg.gtex)[i]]]=tempOutput.mutation
  
}

library(openxlsx)
write.xlsx(pathway.test.gtex, "~/gtex.sexbiased.pathways.xlsx") 


####GnRH secretion####
library(forcats)

GnRH_secretion=data.frame(tissue=cancer.type.gtex,p.value=NA,fc=NA)
for(i in 1:length(cancer.type.gtex)){
  GnRH_secretion$p.value[i]=pathway.test.gtex[[cancer.type.gtex[i]]]$P.Value[pathway.test.gtex[[cancer.type.gtex[i]]]$name=="GnRH secretion - Homo sapiens (human)"]
  GnRH_secretion$fc[i]=pathway.test.gtex[[cancer.type.gtex[i]]]$logFC[pathway.test.gtex[[cancer.type.gtex[i]]]$name=="GnRH secretion - Homo sapiens (human)"]
}
GnRH_secretion$lg.p=-log10(GnRH_secretion$p.value)
GnRH_secretion$fc=-GnRH_secretion$fc

GnRH_secretion %>%
  mutate(tissue = fct_reorder(tissue, dplyr::desc(fc))) %>%
  ggplot( aes(x = tissue,y=fc,fill = tissue))+
  geom_col() +
  theme_classic()+
  theme(legend.position="none")+
  scale_y_continuous(breaks=seq(-0.04,0.11,0.02), expand = c(0, 0),limits = c(-0.04,0.11)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

####GnRH signaling####
GnRH_signaling=data.frame(tissue=cancer.type.gtex,p.value=NA,fc=NA)
for(i in 1:length(cancer.type.gtex)){
  GnRH_signaling$p.value[i]=pathway.test.gtex[[cancer.type.gtex[i]]]$P.Value[pathway.test.gtex[[cancer.type.gtex[i]]]$name=="GnRH signaling pathway - Homo sapiens (human)"]
  GnRH_signaling$fc[i]=pathway.test.gtex[[cancer.type.gtex[i]]]$logFC[pathway.test.gtex[[cancer.type.gtex[i]]]$name=="GnRH signaling pathway - Homo sapiens (human)"]
}
GnRH_signaling$lg.p=-log10(GnRH_signaling$p.value)
GnRH_signaling$fc=-GnRH_signaling$fc

GnRH_signaling %>%
  mutate(tissue = fct_reorder(tissue, dplyr::desc(fc))) %>%
  ggplot( aes(x = tissue,y=fc,fill = tissue))+
  geom_col() +
  theme_classic()+
  theme(legend.position="none")+
  scale_y_continuous(breaks=seq(-0.04,0.11,0.02), expand = c(0, 0),limits = c(-0.04,0.11)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


