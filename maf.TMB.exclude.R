
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ape)
library(purrr)

maf.all = fread("~/mc3.v0.2.8.PUBLIC.maf")
maf.all = maf.all[, c(1, 9, 16, 35)]
maf.all$Tumor_Sample_Barcode = substr(maf.all$Tumor_Sample_Barcode, 1, 12)
colnames(maf.all)[3] = "patient"
clinical.joint = readRDS("~/clinical.joint.rds")

maf = left_join(maf.all, clinical.joint[,c(1,8)], by = "patient")
maf = na.omit(maf)
maf=maf[!duplicated(maf),]
cancer.type = na.omit(unique(maf$disease))

exclude.sample = NA
for (i in 1:length(cancer.type)) {
  cancer.maf = maf[maf$disease == cancer.type[i], ]
  mutation.count = na.omit(count(cancer.maf, patient))
  colnames(mutation.count)[2] = "mutations"
  q.exclude = quantile(mutation.count$mutations, probs = 0.95)
  exclude.temp = mutation.count$patient[mutation.count$mutations > q.exclude]
  exclude.sample = c(exclude.sample, exclude.temp)
  
}
exclude.sample=na.omit(exclude.sample)
saveRDS(exclude.sample,"~/exclude.sample.rds")

