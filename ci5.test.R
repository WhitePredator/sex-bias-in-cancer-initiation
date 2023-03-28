library(ggplot2)
setwd("~/CI5/")

name = list.files("CI5-Xd/")
dir = paste("~/CI5/CI5-Xd/",name,sep="")
n = length(dir)
csv_list=list()
for (i in 1:n) {
  csv_list[[i]]=read.csv(file=dir[i],header = F)
}

registry <- read.delim2("~/CI5/CI5-Xd/registry.txt", header=FALSE)

country_list=list()

canada_num=34
canada=csv_list[[canada_num]]
country_list[["canada"]]=canada

usa.npcr42.white.num=209
usa.npcr42.white=csv_list[[usa.npcr42.white.num]]
country_list[["usa_white"]]=usa.npcr42.white

usa.npcr42.black.num=211
usa.npcr42.black=csv_list[[usa.npcr42.black.num]]
country_list[["usa_black"]]=usa.npcr42.black

south_korea_num=253
south_korea=csv_list[[south_korea_num]]
country_list[["south_korea"]]=south_korea

england_num=392
england=csv_list[[england_num]]
country_list[["england"]]=england


three_col=csv_list[[1]][,1:3]


#china
china_num=215:228
china=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in china_num) {
  a=csv_list[[i]][,4:5]
  china=china+a
}
china=cbind(three_col,china)
country_list[["china"]]=china

#japan
japan_num=245:252
japan=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in japan_num) {
  a=csv_list[[i]][,4:5]
  japan=japan+a
}
japan=cbind(three_col,japan)
country_list[["japan"]]=japan

#france
france_num=299:309
france=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in france_num) {
  a=csv_list[[i]][,4:5]
  france=france+a
}
france=cbind(three_col,france)
country_list[["france"]]=france

#germany
germany_num=310:318
germany=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in germany_num) {
  a=csv_list[[i]][,4:5]
  germany=germany+a
}
germany=cbind(three_col,germany)
country_list[["germany"]]=germany

#italy
italy_num=321:353
italy=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in italy_num) {
  a=csv_list[[i]][,4:5]
  italy=italy+a
}
italy=cbind(three_col,italy)
country_list[["italy"]]=italy

#spain
spain_num=368:380
spain=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in spain_num) {
  a=csv_list[[i]][,4:5]
  spain=spain+a
}
spain=cbind(three_col,spain)
country_list[["spain"]]=spain


country_dict=matrix(data = NA,nrow =length(country_list) ,ncol = 1)
for (i in 1:length(country_list)) {
  country_dict[i,1]=names(country_list[i])
}


t=c(32,37,42,47,52,57,62,67)
country_cancer_list=list()
cancer_num=16   ###  关注的癌种数量
country_num=length(country_dict)
cancer_dict=matrix(data = 0, nrow = cancer_num, ncol = 3  )    
cancer_dict[,1]=1:cancer_num

row.names(cancer_dict)=c("bladder", "colon", "kidney", "liver", "luad", "lusc", "melanoma", "pancreas", "rectum", "stomach","thyroid","brain","oesophagus","Oral","Gallbladder","Larynx")
cancer_dict[,2]=c(163,42,160,59,79,78,100,70,49,35,200,194,23,3,68,75)
cancer_dict[,3]=c(163,42,160,59,79,78,100,70,49,35,200,194,23,3,68,75)


colnames(cancer_dict)=c("dict","ci5_start","ci5_end")   


for (j in 1:country_num) {
  country.matrix=matrix(data=0, ncol=length(t), nrow=cancer_num*2)
  country.matrix=data.frame(country.matrix)
  rownames(country.matrix)[seq(1,cancer_num*2,2)]=paste(rownames(cancer_dict),"_male",sep = "")
  rownames(country.matrix)[seq(2,cancer_num*2,2)]=paste(rownames(cancer_dict),"_female",sep = "")
  population_count=country_list[[j]]$V5[1:38]
  for (i in 1:cancer_num) {
    #####male
    incidence_count=matrix(data = 0,ncol = 1,nrow = length(t))
    for (k in (cancer_dict[i,2]:cancer_dict[i,3])) {
      a=country_list[[j]][((k-1)*38+7):((k-1)*38+14),4]
      incidence_count=incidence_count+a
    }
    country.matrix[2*i-1,]=incidence_count/population_count[7:14]
    #####female
    incidence_count=matrix(data = 0,ncol = 1,nrow = length(t))
    for (k in (cancer_dict[i,2]:cancer_dict[i,3])) {
      a=country_list[[j]][((k-1)*38+19+7):((k-1)*38+19+14),4]
      incidence_count=incidence_count+a
    }
    country.matrix[2*i,]=incidence_count/population_count[26:33]
  }
  country_cancer_list[[country_dict[j]]]=country.matrix
}


#### parallel  test

sex=c(rep(1,8),rep(0,8))
t_log=log10(rep(t,2))
country_cancer_test_list=list()


for (j in 1:country_num) {
  test_matrix=matrix(data=0,nrow = cancer_num ,ncol = 1)
  colnames(test_matrix)="p_value"
  rownames(test_matrix)=rownames(cancer_dict)
  for (i in 1:cancer_num) {
    cancer_log=log10(rbind(country_cancer_list[[j]][2*i-1,],country_cancer_list[[j]][2*i,])+10^-100)  
    cancer_log=as.vector(t(cancer_log))
    fit=summary(lm(cancer_log~t_log+sex+sex*t_log))
    test_matrix[i,1]=fit[[4]][4,4] 
    
  }
  eval(parse(text=paste(paste('country_cancer_test_list[["',country_dict[j],'"]]=',sep=''), 'test_matrix')))
}

cancer_country_test_matrix=matrix(data=0,nrow =cancer_num*country_num ,ncol = 1)

for (k in 1:(cancer_num*country_num)) {
  if((k%%country_num)!=0){
    cancer_country_test_matrix[k,]=country_cancer_test_list[[k%%country_num]][ceiling(k/(country_num)),]
  }
  else{
    cancer_country_test_matrix[k,]=country_cancer_test_list[[country_num]][ceiling(k/(country_num)),]
  }
}
rownames(cancer_country_test_matrix)=paste(rep(rownames(cancer_dict),each=country_num),'_',rep(country_dict,cancer_num),sep='')
colnames(cancer_country_test_matrix)=colnames(test_matrix)
cancer_country_test_matrix=data.frame(cancer_country_test_matrix)
cancer_country_test_matrix$p_value_adjusted=p.adjust(cancer_country_test_matrix$p_value,method = "bonferroni")


write.csv(cancer_country_test_matrix,"~/unparallel.test.csv")

