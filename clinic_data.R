library(TCGAbiolinks)
library(data.table)
library(dplyr)

clinical=fread("~/Survival_SupplementalTable_S1_20171025_xena_sp.gz")
cancer.type=unique(clinical$`cancer type abbreviation`)
cancer.type.tcga=paste("TCGA-",cancer.type,sep = "")

clinical.biolinks=GDCquery_clinic(project = "TCGA-ACC", type = "clinical")

for (i in 2:length(cancer.type.tcga)) {
  clinical.temp=GDCquery_clinic(project = cancer.type.tcga[i], type = "clinical")
  clinical.biolinks=rbind(clinical.biolinks,clinical.temp)
}

saveRDS(clinical.biolinks,"~/clinical_tcgabiolinks.rds")

clinical.joint=clinical.biolinks[,c(1,5,6,7,25,28,31,33,42)]
clinical.joint$sample=paste(clinical.joint$submitter_id,"-01",sep = "")
clinical.joint$alcohol_history[is.na(clinical.joint$alcohol_history)]="not.reported"
clinical.joint$years_smoked=ifelse(is.na(clinical.joint$years_smoked),"not.reported","yes")
colnames(clinical.joint)=c("patient","primary.diagnosis","tumor.stage","age","alcohol","smoke","gender","race","disease","sample")
clinical.joint.list=split(clinical.joint,clinical.joint$disease)
for(i in 1:33){
  clinical.joint.list[[i]][is.na(clinical.joint.list[[i]]$age),"age"]=mean(na.omit(clinical.joint.list[[i]]$age))
}
clinical.joint.1=purrr::reduce(clinical.joint.list,bind_rows)

clinic <- fread("~/TcgaTargetGTEX_phenotype.txt.gz")
project <- fread("~/TCGA_GTEX_category.txt")
melanoma.sample=clinic$sample[clinic$detailed_category=="Skin Cutaneous Melanoma"]
melanoma.sample.06=melanoma.sample[grepl("-06$",melanoma.sample)]
clinic=clinic[clinic$"_sample_type" %in% c("Normal Tissue","Primary Tumor") | clinic$sample %in% melanoma.sample.06,]
gender.none=as.character(clinic$sample[clinic$`_gender`!="Male" & clinic$`_gender`!="Female"])
clinic=clinic[!clinic$sample %in% gender.none,]
clinic$patient=substr(clinic$sample,1,12)
clinical.joint=left_join(clinical.joint.1[1:9],clinic[,c(1:5,8)],by="patient")
clinical.joint=clinical.joint[,c(1,3:10)]
clinical.joint[is.na(clinical.joint$sample),"sample"]=paste(clinical.joint[is.na(clinical.joint$sample),"patient"],"-01",sep = "")


saveRDS(clinical.joint,"~/clinical.joint.rds")
