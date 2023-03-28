rm(list = ls())
gc()
library(ape)
library(purrr)
library(dplyr)
library(data.table)
library(tidyr)


maf=fread("~/mc3.v0.2.8.PUBLIC.maf")
maf.sample=unique(substr(maf$Tumor_Sample_Barcode,1,15))
maf.sample.01=maf.sample[grepl("-01$",maf.sample)]
maf.sample.02356=maf.sample[grepl("-02$",maf.sample) | grepl("-03$",maf.sample) | grepl("-05$",maf.sample) | grepl("-06$",maf.sample)]
maf.duplicate=maf.sample.02356[substr(maf.sample.02356,1,12) %in% substr(maf.sample.01,1,12)]

gff3.37=read.gff("~/Homo_sapiens.GRCh37.87.gff3.gz")
exclude.sample <- readRDS("~/exclude.sample.rds")
clinical.biolinks=readRDS("~/clinical_tcgabiolinks.rds")

clinical.joint=clinical.biolinks[,c(1,5,6,7,25,28,31,33,42)]
clinical.joint$alcohol_history[is.na(clinical.joint$alcohol_history)]="not.reported"
clinical.joint$years_smoked=ifelse(is.na(clinical.joint$years_smoked),"not.reported","yes")
colnames(clinical.joint)[1]="patient"
cancer.type=unique(clinical.joint$disease)
nrac.type=c("ACC",  "BLCA", "CHOL", "COAD", "DLBC", "ESCA", "GBM",  "HNSC", "KICH", "KIRC", "KIRP",
            "LIHC", "LUAD", "LUSC", "MESO", "PAAD", "PCPG", "READ", "SARC", "SKCM",
            "STAD", "THCA", "THYM", "UVM" )

gff3.37.gene=separate_rows(gff3.37,attributes,sep = ";")
gff3.37.gene=gff3.37.gene[grep("Name",gff3.37.gene$attributes),]
gff3.37.gene$attributes=substr(gff3.37.gene$attributes,6,100)
gff3.37.gene=gff3.37.gene[,c(9,1,4,5,7,3)]
colnames(gff3.37.gene)=c("gene","chromosome_name","start_position","end_position","stand","type")
chr.length=gff3.37[gff3.37$type=="chromosome",]
autosome.length=sum(chr.length$end[1:22])
x.chr.length=chr.length$end[24]
y.chr.length=59373566  ###https://grch37.ensembl.org/info/genome/genebuild/assembly.html
male.length=2*autosome.length+x.chr.length+y.chr.length
female.length=2*(autosome.length+x.chr.length)


maf.tmp=maf[,c(1,9,16,35,72,73)]
maf.tmp$Tumor_Sample_Barcode=substr(maf.tmp$Tumor_Sample_Barcode,1,15)
maf.tmp=maf.tmp[!maf.tmp$Tumor_Sample_Barcode %in% maf.duplicate,]
colnames(maf.tmp)[c(1,3)]=c("gene","patient")

######
maf.tmp$gene[maf.tmp$gene=="GPR98"]="ADGRV1"
gff3.37.gene$gene[gff3.37.gene$gene=="GPR98"]="ADGRV1"
######
maf.tmp$patient=substr(maf.tmp$patient,1,12)
maf.tmp=maf.tmp[!maf.tmp$patient %in% exclude.sample,]


maf.fuc.lost=maf.tmp[maf.tmp$Variant_Classification %in% 
                       c("Nonsense_Mutation","Nonstop_Mutation","Translation_Start_Site",
                         "Frame_Shift_Ins","Frame_Shift_Del"),]

maf.sift.polyphen=maf.tmp[grepl("deleterious\\(",maf.tmp$SIFT) & grepl("probably_damaging",maf.tmp$PolyPhen),]

maf.repeat=maf.tmp
maf.repeat=maf.repeat[grepl(">",maf.repeat$HGVSc),]
maf.repeat$mutation.loci=NA
grepsub=function(x){
  x[7]=substr(x[4],1,grepRaw(">",x[4])-2)
}
maf.repeat$mutation.loci=apply(maf.repeat,1,grepsub)

maf.repeat$tmp=paste(maf.repeat$gene,maf.repeat$mutation.loci,sep="-")
tmp.count=data.frame(table(maf.repeat$tmp))
maf.high.freq=maf.repeat[maf.repeat$tmp %in% tmp.count$Var1[tmp.count$Freq>=5],]

maf.all=rbind(maf.fuc.lost,maf.sift.polyphen)
maf.all=rbind(maf.all,maf.high.freq[,-(7:8)])
maf.all=maf.all[!duplicated(maf.all),]
maf.all=maf.all[maf.all$Variant_Classification!="Silent",]

maf.all=maf.all[,c(1,2,3)]
maf.all=left_join(maf.all,clinical.joint[,c(1,7,8,9)],by="patient")
maf.all=maf.all[maf.all$disease %in% nrac.type,]


maf.genes=data.frame(gene=unique(maf.all$gene))
maf.genes=left_join(maf.genes,gff3.37.gene,by="gene")
maf.genes=maf.genes[!grepl("GL",maf.genes$chromosome_name),]
par1=c(10001,2781479)
par2=c(155701383,156030895)
maf.genes$par=ifelse(maf.genes$chromosome_name=="X" & 
                       ((maf.genes$start_position>=par1[1] &maf.genes$start_position<=par1[2]) | 
                          (maf.genes$start_position>=par2[1] & maf.genes$start_position<=par2[2])),1,0)

maf.genes$gene.copy[!maf.genes$chromosome_name %in% c("X","Y") | maf.genes$par==1 ]=2
maf.genes$gene.copy[maf.genes$chromosome_name=="Y" ]=1
maf.genes$gene.copy[is.na(maf.genes$gene.copy) ]=1.5

maf.all=left_join(maf.all,maf.genes[,c(1,2)],by="gene")


maf.all=maf.all[!duplicated(maf.all[,c(1,3)]),]

maf.all.list=split(maf.all,list(maf.all$disease))


get.gene.test=function(gene=gene,cancer.maf=cancer.maf,cancer.patient=cancer.patient,gender.table=gender.table,side=side){
  
  gene.mutation=cancer.patient
  gene.mutation$mutation=ifelse(gene.mutation$patient %in% unique(cancer.maf$patient[cancer.maf$gene==gene]),1,0)
  gene.copy=maf.genes$gene.copy[maf.genes$gene==gene]
  
  ####LRT#####
  gene.mutation.table=table(gene.mutation[,c("gender","mutation")])
  if(ncol(gene.mutation.table)!=1){
  
    a=gene.mutation.table[1,2]
    b=gene.mutation.table[1,1]
    c=gene.mutation.table[2,2]
    d=gene.mutation.table[2,1]
    gene.chr=maf.genes$chromosome_name[maf.genes$gene==gene]
    lamda=gender.table$Freq[gender.table$gender=="female" & gender.table$chromosome_name==gene.chr]/
      gender.table$Freq[gender.table$gender=="male" & gender.table$chromosome_name==gene.chr]
    lamda=ifelse(is.infinite(lamda),1,lamda)
    
    ####fisher exact test ####
    fisher.p=fisher.test(gene.mutation.table)$p.value
    
    if (maf.genes$gene.copy[maf.genes$gene==gene]!=1) {
      pxy=c/(c+d)
      pxx=a/(a+b)
      likelihood.univ=1e200*(1-pxx)^b*(pxx)^a*(1-pxy)^d*(pxy)^c
      
      e=lamda*(a+b+c+d)
      f=-(a*lamda+a+b*lamda+c+c*lamda+d)
      g=a+c
      x1=(-f-(f^2-4*e*g)^(1/2))/(2*e)
      x2=(-f+(f^2-4*e*g)^(1/2))/(2*e)
      if ((f^2-4*e*g)>0 & x1>=0 & x1<=1) {
        x.root=x1
      }else{
        stop("quadratic equation likelihood has some problems")
      }
      likelihood.rest=1e200*(lamda*x.root)^a*(1- lamda*x.root)^b*(x.root)^c*(1-x.root)^d
      
      statistic=2*(log(likelihood.univ)-log(likelihood.rest))
      pvalue.ng=1-pchisq(statistic,df = 1)
    }else{
      pvalue.ng=fisher.p
      lamda=0
    }
    
    
    expect.male=(c+d)*x.root
    expect.female=(a+b)*lamda*x.root
    if(a>=expect.female & c<=expect.male){
      direction="female"
    }else if(a<=expect.female & c>=expect.male){
      direction="male"
    }else{
      direction="none"
    }
    
    ####logistic regression####
    gene.mutation=data.frame(patient=unique(cancer.maf$patient))
    gene.mutation$mutation=ifelse(gene.mutation$patient %in% unique(cancer.maf$patient[cancer.maf$gene==gene]),1,0)
    gene.mutation=left_join(gene.mutation,clinical.joint,by="patient")
    colnames(gene.mutation)[5]="age"
    gene.mutation[is.na(gene.mutation$age),"age"]=mean(na.omit(gene.mutation$age))
    gene.mutation$gender=ifelse(gene.mutation$gender=="male",1,0)
    # dummy.feature <- setdiff(colnames(gene.mutation),c("patient","mutation","gender","age"))
    dummy.feature <- c("tumor_stage","alcohol_history","years_smoked","race")  
    lr.dummy <- dummy.data.frame(gene.mutation, names=dummy.feature)
    dummy.list <- attr(lr.dummy,"dummies")
    if (length(dummy.list)!=0) {
      rm.col <- c()
      for (k in 1:length(dummy.list)){
        rm.col <- c(rm.col, dummy.list[[k]][length(dummy.list[[k]])])
      }
      lr.dummy=lr.dummy[,-rm.col]
    }
    colnames(lr.dummy) <- gsub(" ", ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub("-", ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub('\\[', ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub('\\]', ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub('\\(', ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub('\\)', ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub(',', ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub(';', ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub('/', ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub('\\+', ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub('&', ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub("'", ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub('\\|', ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub('>', ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub('=', ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub('%', ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub('<', ".", colnames(lr.dummy))
    colnames(lr.dummy) <- gsub(':', ".", colnames(lr.dummy))
    
    exclude.col <- match(c("primary_diagnosis","disease","patient","mutation"), colnames(lr.dummy))
    form <- as.formula(paste("mutation ~ ",paste(colnames(lr.dummy)[-exclude.col],collapse="+"),sep=""))
    fit=glm(form, data = lr.dummy, family  = binomial(link ="logit"))
    fit.aic=step(fit,scope=list(lower=as.formula(mutation ~ gender), upper=as.formula(mutation ~ .)), trace = 0,steps = step_a,direction = side)
    lr.pvalue=summary(fit.aic)$coefficients["gender",4]
    
    
    pvalue.return=data.frame(gene=gene,  lrt.ng=pvalue.ng,  lr.p=lr.pvalue, 
                             male1=gene.mutation.table[2,2],male0=gene.mutation.table[2,1],
                             female1=gene.mutation.table[1,2],female0=gene.mutation.table[1,1],
                             expect.male=expect.male,expect.female=expect.female,lamda=lamda,
                             direction=direction)
    
    return=list(pvalue.return=pvalue.return,gene.mutation=gene.mutation)
    return(return)
  }
  
}

get.test.results=function(cancer.maf=cancer.maf,
                          clinical.joint=clinical.joint,maf.genes=maf.genes,
                          side=side){
  if(cancer.maf$disease[1]=="SKCM"){
    cancer.maf=cancer.maf[cancer.maf$race=="white",]
  }
  gene.freq=data.frame(table(cancer.maf$gene))
  genes.005=as.character(gene.freq$Var1[gene.freq$Freq>=0.05*length(unique(cancer.maf$patient))])
  cancer.patient=data.frame(patient=unique(cancer.maf$patient))
  cancer.patient=left_join(cancer.patient,clinical.joint[,c(1,7)],by="patient")
  
  gender.table=data.frame(table(cancer.maf[,c(4,7)]))
  gender.table=gender.table[!grepl("GL",gender.table$chromosome_name),]
  
  library(snowfall)
  library(parallel)
  sfInit(parallel = TRUE, cpus = detectCores()-2)
  sfExport("get.gene.test")
  sfExport("clinical.joint")
  sfExport("maf.genes")
  sfExport("cancer.maf")
  sfExport("cancer.patient")
  sfExport("gender.table")
  sfExport("step_a")
  sfExport("side")
  sfLibrary(dplyr)
  sfLibrary(purrr)
  sfLibrary(tidyr)
  sfLibrary(maxLik)
  sfLibrary(dummies)
  test.matrix=sfSapply(genes.005, get.gene.test, cancer.maf, cancer.patient, gender.table,side, simplify = FALSE)
  mutation.list=purrr::transpose(test.matrix)[[2]]
  test.matrix=reduce(purrr::transpose(test.matrix)[[1]],bind_rows)
  test.matrix=left_join(test.matrix,maf.genes[,c(1,2,6,7)],by="gene")
  test.matrix$cancer=cancer.maf$disease[1]
  test.matrix$lrt.ng.fdr=p.adjust(test.matrix$lrt.ng,method = "fdr")
  return=list(mutation.list=mutation.list,test.matrix=test.matrix)
  return(return)
}

library(snowfall)
library(parallel)
sfInit(parallel = TRUE, cpus = detectCores()-2)
sfExport("get.gene.test")
sfExport("get.test.results")
sfLibrary(dplyr)
sfLibrary(purrr)
step_a=100
side="both"
sfExport("step_a")
sfExport("side")
test.result=sfLapply(maf.all.list,get.test.results,clinical.joint,maf.genes,side)
mutation.list=purrr::transpose(test.result)[[1]]
gene.result=purrr::transpose(test.result)[[2]]
gene.result.reduce=purrr::reduce(gene.result,bind_rows)
test.significant=gene.result.reduce[gene.result.reduce$lrt.ng.fdr<=0.05,]
test.significant=test.significant[test.significant$lr.p<=0.05,]
sfStop()

write.csv(gene.result.reduce,"~/bias.gene.csv")

