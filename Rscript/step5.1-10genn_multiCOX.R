rm(list = ls())
setwd("H:/210E盘文档/colon cancer/0DataAndScript")

load("1.primary data/train/ESP.RData")
load("1.primary data/train/hub_probe.RData")
{
  #GSE39582 drug related information
  GSE39582_clin_all <- read.table("4.ACT/0.data/GSE39582clinical.txt", sep="\t",header=T,na.strings = c(" ","NA"),stringsAsFactors = F)
  GSE39582_clin_mid <- GSE39582_clin_all[which(GSE39582_clin_all$GEOID%in%colnames(GSE39582_exp_matrix)),]
  drug_all_1 <- GSE39582_clin_mid[which(GSE39582_clin_mid$chemotherapy_adjuvant==" Y"),]
  ##使用FU的联合用药
  drug_all <- drug_all_1[complete.cases(drug_all_1$chemotherapy_adjuvant_type),]
  ##使用FU的联合用药
  drug_5fu_FUFOL_FOLFOX_FOLFIRI <- drug_all[which(drug_all$chemotherapy_adjuvant_type==" 5FU"|drug_all$chemotherapy_adjuvant_type==" FUFOL"|drug_all$chemotherapy_adjuvant_type==" FOLFOX"|drug_all$chemotherapy_adjuvant_type==" FOLFIRI"),]
  
  GSE39582.ACT.sample = drug_5fu_FUFOL_FOLFOX_FOLFIRI$GEOID
}

##re-construct RSF prognosis model according to the RFS time of GSE39582 with 5-FU ACT
library(randomForestSRC)
library(survival)
library(Hmisc)
#select the sample of training set
hub_genes <- as.character(hub_probe$gene[1:18])

hub_genes
###
hub_genes <- c("THY1","BUB1","DEPDC1","STON1","ATAD2","HSD17B2","TPX2","AURKB","CYR61","FAM84A")

#hub_genes <- c(hub_genes,"CYR61","FAM84A")
dim(drug_all_1) #203  21
dim(drug_all) #144  21
dim(drug_5fu_FUFOL_FOLFOX_FOLFIRI) #142  21

sample_5FU <- drug_5fu_FUFOL_FOLFOX_FOLFIRI$GEOID
GSE39582.RFS <- GSE39582_exp_matrix[hub_genes,sample_5FU]
#prognosis model
time <- clinical_infor_df[sample_5FU,]$time
event <- clinical_infor_df[sample_5FU,]$event
train.matrix <- na.omit(data.frame(time = time, event = event, t(GSE39582.RFS), stringsAsFactors = F))
#select the optimal seed?
set.seed(1244746)
grow.RFS=rfsrc(Surv(time,event) ~.,data =train.matrix,ntree = 1000,nsplit = 2)
grow.RFS

library(survminer)
#######use the model.grow$predicted, which is in-bag predicted value, 
#######to apply maxstat, then use the cutoff value on testing data;
GSE39582_rfs_score = grow.RFS$predicted
summary(GSE39582_rfs_score)
GSE39582_rfs_score.surv = cbind(train.matrix[,1:2],GSE39582_rfs_score)
head(GSE39582_rfs_score.surv)
#Determining the optimal cutpoint for GSE39582_rfs_score
GSE39582_rfs_score.cut = surv_cutpoint(
  data = GSE39582_rfs_score.surv,
  time = "time",
  event = "event",
  variables = "GSE39582_rfs_score",
  progressbar = F
)
summary(GSE39582_rfs_score.cut)
cut_off = GSE39582_rfs_score.cut$cutpoint$cutpoint
#Plot the cutpoint for GSE39582_rfs_score
plot(GSE39582_rfs_score.cut, "GSE39582_rfs_score", palette = "npg")
GSE39582_rfs_score.cat = surv_categorize(GSE39582_rfs_score.cut)
GSE39582_rfs_score.cat$GSE39582_rfs_score = plyr::revalue(GSE39582_rfs_score.cat$GSE39582_rfs_score, 
                                                          c("low"="ACT-benefit", "high"="ACT-nonbenefit"))


all(rownames(GSE39582_rfs_score.surv) == rownames(GSE39582_rfs_score.cat))
GSE39582.info.ACT = cbind(GSE39582_rfs_score.cat, GSE39582_rfs_score.surv$GSE39582_rfs_score)
colnames(GSE39582.info.ACT)[3:4] = c("group","score")
#save(GSE39582.info.ACT, file = "Script/2022.11.03/GSE39582.10gene.RData")



##################GSE17536 testing cohort
load("1.primary data/test/GSE17536.RData")
test_matrix <- test_matrix[probe2gene_loc,]
test_matrix <- test_matrix[sig_probe_loc,]
rownames(test_matrix) <- all.gene #assign gene symbols to matrix
#drug related information##########提取GSE17536中期使用过药物的样本
GSE17536_mid_drug<-read.table("4.ACT/0.data/GSE17536_mid_drug.txt",sep="\t",header=T,stringsAsFactors = F)

GSE17536_sample_ACT <- GSE17536_mid_drug$geoid
test_cli <- test_cli[GSE17536_sample_ACT,]
#predict testing set
test <- data.frame(t(test_matrix[hub_genes,GSE17536_sample_ACT]),stringsAsFactors = F)
pred <- predict(grow.RFS, newdata=test)
score <- pred$predicted
group <- rep(0, length(score))
for (i in 1:length(group)) {
  if( score[i]<=cut_off) {
    group[i]<-"low"
  } else{
    group[i]<-"high"
  }
}

{
  ##5-years
  test_cli[which(test_cli[,1]>60),2] = 0
  test_cli[which(test_cli[,1]>60),1] = 60
}

surv.df <- data.frame(test_cli[,1], test_cli[,2], score,stringsAsFactors = F)
colnames(surv.df)[1:2] <- c("surtime", "censor")
test.survdiff <- survdiff(Surv(surtime,censor)~group,data=surv.df)
HR_temp <- (test.survdiff$obs[2]/test.survdiff$exp[2])/(test.survdiff$obs[1]/test.survdiff$exp[1])
HR_temp <- round(HR_temp,4)
p.val_temp <- round(1 - pchisq(test.survdiff$chisq, 1),4)
up95_temp <- exp(log(as.numeric(HR_temp)) + qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
up95_temp <- round(up95_temp,4)
low95_temp <- exp(log(as.numeric(HR_temp)) - qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
low95_temp <- round(low95_temp,4)
c_index_temp <- rcorr.cens(-surv.df$score, Surv(surv.df$surtime, surv.df$censor))[1]
c_index_temp <- round(c_index_temp,4)
plot(survfit(Surv(surtime,censor)~group,data=surv.df),xlab="Survival in months", ylab="Survival probability",
     col=c("red","darkblue"),bty="l",lwd=4)
title(main=paste("GSE17536", "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
legend("topright", legend=c(paste0("high (N=",table(group)["high"],")"),paste0("low (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                             paste0("p valu:",p.val_temp),
                             paste0("HR:",HR_temp),
                             paste0("up95:",up95_temp),
                             paste0("low95:",low95_temp)))


GSE17536_clinical <- read.table(paste0("1.primary data/test/","GSE17536",".frma.sample.tsv"),
                                header = T,sep = "\t",stringsAsFactors = F)
rownames(GSE17536_clinical) <- GSE17536_clinical$Sample
Age <- GSE17536_clinical[rownames(test_cli),"Age"]
Sex <- GSE17536_clinical[rownames(test_cli),"Sex"]
Stage <- factor(GSE17536_clinical[rownames(test_cli),"Stage"], levels = c(2,3))
Group <- factor(group,levels=c("low","high"))

cox.info.gse17536 <- cbind(test_cli,Age) %>% cbind(Sex) %>% cbind(Group) %>% cbind(Stage)
## group
uni.cox.group <- coxph(Surv(time, event) ~ Group, data = cox.info.gse17536)
summary(uni.cox.group)
## age
uni.cox.age <- coxph(Surv(time, event) ~ Age, data = cox.info.gse17536)
summary(uni.cox.age)
## sex
uni.cox.sex <- coxph(Surv(time, event) ~ Sex, data = cox.info.gse17536)
summary(uni.cox.sex)
####
uni.cox.stage <- coxph(Surv(time, event) ~ Stage, data = cox.info.gse17536)
summary(uni.cox.stage)
## age+sex+group+stage
mul.cox <- coxph(Surv(time, event) ~ Group+Age+Sex+Stage, data = cox.info.gse17536)
summary(mul.cox)
#基础森林图
ggforest(mul.cox,  #coxph得到的Cox回归结果
         data = cox.info.gse17536,  #数据集
         main = 'Hazard ratio of GSE17536',  #标题
         cpositions = c(0.05, 0.15, 0.35),  #前三列距离
         fontsize = 1, #字体大小
         refLabel = 'reference', #相对变量的数值标签，也可改为1
         noDigits = 3 #保留HR值以及95%CI的小数位数
)


##################GSE29621
library(dplyr)
GSE29621.sample <- read.table("1.primary data/test/GSE29621.sample.txt",
                              header = T,sep = "\t",stringsAsFactors = F)
GSE29621.sample.23 <- filter(GSE29621.sample,stage==2|stage==3)
GSE29621.sample.23 <- na.omit(GSE29621.sample.23)
GSE29621.sample.drug <- intersect(GSE29621.sample.23$gse17536,GSE17536_sample_ACT)
GSE29621.sample.drug21 <- GSE29621.sample.23$sample[which(GSE29621.sample.23$gse17536 %in% GSE29621.sample.drug)]
load("1.primary data/test/GSE29621.RData")
test_matrix <- test_matrix[probe2gene_loc,]
test_matrix <- test_matrix[sig_probe_loc,]
rownames(test_matrix) <- all.gene #assign gene symbols to matrix
test_cli <- test_cli[GSE29621.sample.drug21,]
#predict testing set
test <- data.frame(t(test_matrix[hub_genes,GSE29621.sample.drug21]),stringsAsFactors = F)
pred <- predict(grow.RFS, newdata=test)
score <- pred$predicted
group <- rep(0, length(score))
for (i in 1:length(group)) {
  if( score[i]<=cut_off) {
    group[i]<-"low"
  } else{
    group[i]<-"high"
  }
}

{
  ##5-years
  test_cli[which(test_cli[,1]>60),2] = 0
  test_cli[which(test_cli[,1]>60),1] = 60
}

surv.df <- data.frame(test_cli[,1], test_cli[,2], score,stringsAsFactors = F)
colnames(surv.df)[1:2] <- c("surtime", "censor")
test.survdiff <- survdiff(Surv(surtime,censor)~group,data=surv.df)
HR_temp <- (test.survdiff$obs[2]/test.survdiff$exp[2])/(test.survdiff$obs[1]/test.survdiff$exp[1])
HR_temp <- round(HR_temp,4)
p.val_temp <- round(1 - pchisq(test.survdiff$chisq, 1),4)
up95_temp <- exp(log(as.numeric(HR_temp)) + qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
up95_temp <- round(up95_temp,4)
low95_temp <- exp(log(as.numeric(HR_temp)) - qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
low95_temp <- round(low95_temp,4)
c_index_temp <- rcorr.cens(-surv.df$score, Surv(surv.df$surtime, surv.df$censor))[1]
c_index_temp <- round(c_index_temp,4)
plot(survfit(Surv(surtime,censor)~group,data=surv.df),xlab="Survival in months", ylab="Survival probability",
     col=c("red","darkblue"),bty="l",lwd=4)
title(main=paste("GSE29621", "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
legend("topright", legend=c(paste0("high (N=",table(group)["high"],")"),paste0("low (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                             paste0("p valu:",p.val_temp),
                             paste0("HR:",HR_temp),
                             paste0("up95:",up95_temp),
                             paste0("low95:",low95_temp)))

rownames(GSE29621.sample.23) <- GSE29621.sample.23$sample
test_cli1 = test_cli
rownames(test_cli1) = GSE29621.sample.23[rownames(test_cli),"gse17536"]

Age <- GSE17536_clinical[rownames(test_cli1),"Age"]
Sex <- GSE17536_clinical[rownames(test_cli1),"Sex"]
Stage <- factor(GSE17536_clinical[rownames(test_cli1),"Stage"], levels = c(2,3))
Group <- factor(group,levels=c("low","high"))

cox.info.gse29621 <- cbind(test_cli1,Age) %>% cbind(Sex) %>% cbind(Group) %>% cbind(Stage)
## group
uni.cox.group <- coxph(Surv(time, event) ~ Group, data = cox.info.gse29621)
summary(uni.cox.group)
## age
uni.cox.age <- coxph(Surv(time, event) ~ Age, data = cox.info.gse29621)
summary(uni.cox.age)
## sex
uni.cox.sex <- coxph(Surv(time, event) ~ Sex, data = cox.info.gse29621)
summary(uni.cox.sex)
####
uni.cox.stage <- coxph(Surv(time, event) ~ Stage, data = cox.info.gse29621)
summary(uni.cox.stage)
## age+sex+group+stage
mul.cox <- coxph(Surv(time, event) ~ Group+Age+Sex+Stage, data = cox.info.gse29621)
summary(mul.cox)
#基础森林图
ggforest(mul.cox,  #coxph得到的Cox回归结果
         data = cox.info.gse29621,  #数据集
         main = 'Hazard ratio of GSE29621',  #标题
         cpositions = c(0.05, 0.15, 0.35),  #前三列距离
         fontsize = 1, #字体大小
         refLabel = 'reference', #相对变量的数值标签，也可改为1
         noDigits = 3 #保留HR值以及95%CI的小数位数
)



##################GSE17537 testing cohort
#drug related information##########提取GSE17537中期使用过药物的样本
load("4.ACT/0.data/GSE17537_matrix.Rdata")
load("1.primary data/test/GSE17537.RData")
test_matrix <- test_matrix[probe2gene_loc,]
test_matrix <- test_matrix[sig_probe_loc,]
rownames(test_matrix) <- all.gene #assign gene symbols to matrix
#drug related information##########提取GSE17537中期使用过药物的样本
GSE17537_sample_ACT <- GSE17537_mid_drug$geoid
test_cli <- test_cli[GSE17537_sample_ACT,]
#predict testing set
test <- data.frame(t(test_matrix[hub_genes,GSE17537_sample_ACT]),stringsAsFactors = F)
pred <- predict(grow.RFS, newdata=test)
score <- pred$predicted
group <- rep(0, length(score))
for (i in 1:length(group)) {
  if( score[i]<=cut_off) {
    group[i]<-"low"
  } else{
    group[i]<-"high"
  }
}

{
  ##5-years
  test_cli[which(test_cli[,1]>60),2] = 0
  test_cli[which(test_cli[,1]>60),1] = 60
}

surv.df <- data.frame(test_cli[,1], test_cli[,2], score,stringsAsFactors = F)
colnames(surv.df)[1:2] <- c("surtime", "censor")
test.survdiff <- survdiff(Surv(surtime,censor)~group,data=surv.df)
HR_temp <- (test.survdiff$obs[2]/test.survdiff$exp[2])/(test.survdiff$obs[1]/test.survdiff$exp[1])
HR_temp <- round(HR_temp,4)
p.val_temp <- round(1 - pchisq(test.survdiff$chisq, 1),4)
up95_temp <- exp(log(as.numeric(HR_temp)) + qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
up95_temp <- round(up95_temp,4)
low95_temp <- exp(log(as.numeric(HR_temp)) - qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
low95_temp <- round(low95_temp,4)
c_index_temp <- rcorr.cens(-surv.df$score, Surv(surv.df$surtime, surv.df$censor))[1]
c_index_temp <- round(c_index_temp,4)
plot(survfit(Surv(surtime,censor)~group,data=surv.df),xlab="Survival in months", ylab="Survival probability",
     col=c("red","darkblue"),bty="l",lwd=4)
title(main=paste("GSE17537", "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
legend("topright", legend=c(paste0("high (N=",table(group)["high"],")"),paste0("low (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                             paste0("p valu:",p.val_temp),
                             paste0("HR:",HR_temp),
                             paste0("up95:",up95_temp),
                             paste0("low95:",low95_temp)))


GSE17537_clinical <- read.table(paste0("1.primary data/test/","GSE17537",".frma.sample.tsv"),
                                header = T,sep = "\t",stringsAsFactors = F)
rownames(GSE17537_clinical) <- GSE17537_clinical$Sample

Age <- GSE17537_clinical[rownames(test_cli),"Age"]
Sex <- GSE17537_clinical[rownames(test_cli),"Sex"]
Stage <- factor(GSE17537_clinical[rownames(test_cli),"Stage"], levels = c(2,3))
Group <- factor(group,levels=c("low","high"))

cox.info.gse17537 <- cbind(test_cli,Age) %>% cbind(Sex) %>% cbind(Group) %>% cbind(Stage)
## group
uni.cox.group <- coxph(Surv(time, event) ~ Group, data = cox.info.gse17537)
summary(uni.cox.group)
## age
uni.cox.age <- coxph(Surv(time, event) ~ Age, data = cox.info.gse17537)
summary(uni.cox.age)
## sex
uni.cox.sex <- coxph(Surv(time, event) ~ Sex, data = cox.info.gse17537)
summary(uni.cox.sex)
####
uni.cox.stage <- coxph(Surv(time, event) ~ Stage, data = cox.info.gse17537)
summary(uni.cox.stage)
## age+sex+group+stage
mul.cox <- coxph(Surv(time, event) ~ Group+Age+Sex+Stage, data = cox.info.gse17537)
summary(mul.cox)


##################GSE17538 testing cohort
load("1.primary data/test/GSE17538.RData")
test_matrix <- test_matrix[probe2gene_loc,]
test_matrix <- test_matrix[sig_probe_loc,]
rownames(test_matrix) <- all.gene #assign gene symbols to matrix
#drug related information##########GSE17538中期使用过药物的样本
GSE17538_sample_ACT <- c(GSE17536_sample_ACT,GSE17537_sample_ACT)
test_cli <- test_cli[GSE17538_sample_ACT,]
#predict testing set
test <- data.frame(t(test_matrix[hub_genes,GSE17538_sample_ACT]),stringsAsFactors = F)
pred <- predict(grow.RFS, newdata=test)
score <- pred$predicted
group <- rep(0, length(score))
for (i in 1:length(group)) {
  if( score[i]<=cut_off) {
    group[i]<-"low"
  } else{
    group[i]<-"high"
  }
}

{
  ##5-years
  test_cli[which(test_cli[,1]>60),2] = 0
  test_cli[which(test_cli[,1]>60),1] = 60
}

surv.df <- data.frame(test_cli[,1], test_cli[,2], score,stringsAsFactors = F)
colnames(surv.df)[1:2] <- c("surtime", "censor")
test.survdiff <- survdiff(Surv(surtime,censor)~group,data=surv.df)
HR_temp <- (test.survdiff$obs[2]/test.survdiff$exp[2])/(test.survdiff$obs[1]/test.survdiff$exp[1])
HR_temp <- round(HR_temp,4)
p.val_temp <- round(1 - pchisq(test.survdiff$chisq, 1),4)
up95_temp <- exp(log(as.numeric(HR_temp)) + qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
up95_temp <- round(up95_temp,4)
low95_temp <- exp(log(as.numeric(HR_temp)) - qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
low95_temp <- round(low95_temp,4)
c_index_temp <- rcorr.cens(-surv.df$score, Surv(surv.df$surtime, surv.df$censor))[1]
c_index_temp <- round(c_index_temp,4)
plot(survfit(Surv(surtime,censor)~group,data=surv.df),xlab="Survival in months", ylab="Survival probability",
     col=c("red","darkblue"),bty="l",lwd=4)
title(main=paste("GSE17538", "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
legend("topright", legend=c(paste0("high (N=",table(group)["high"],")"),paste0("low (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                             paste0("p valu:",p.val_temp),
                             paste0("HR:",HR_temp),
                             paste0("up95:",up95_temp),
                             paste0("low95:",low95_temp)))


GSE17538_clinical <- rbind(GSE17536_clinical,GSE17537_clinical)


Age <- GSE17538_clinical[rownames(test_cli),"Age"]
Sex <- GSE17538_clinical[rownames(test_cli),"Sex"]
Stage <- factor(GSE17538_clinical[rownames(test_cli),"Stage"], levels = c(2,3))
Group <- factor(group,levels=c("low","high"))

cox.info.gse17538 <- cbind(test_cli,Age) %>% cbind(Sex) %>% cbind(Group) %>% cbind(Stage) %>% cbind(score)
#
cox.info.gse17538.ACT = cox.info.gse17538
cox.info.gse17538.ACT$Group = plyr::revalue(cox.info.gse17538.ACT$Group, c("high"="ACT-nonbenefit", "low"="ACT-benefit"))
save(cox.info.gse17538.ACT, file = "Script/2022.11.01/GSE17538.10gene.RData")
## group
uni.cox.group <- coxph(Surv(time, event) ~ Group, data = cox.info.gse17538)
summary(uni.cox.group)
## score
uni.cox.score <- coxph(Surv(time, event) ~ score, data = cox.info.gse17538)
summary(uni.cox.score)
## age
uni.cox.age <- coxph(Surv(time, event) ~ Age, data = cox.info.gse17538)
summary(uni.cox.age)
## sex
uni.cox.sex <- coxph(Surv(time, event) ~ Sex, data = cox.info.gse17538)
summary(uni.cox.sex)
####
uni.cox.stage <- coxph(Surv(time, event) ~ Stage, data = cox.info.gse17538)
summary(uni.cox.stage)
## age+sex+group+stage
mul.cox <- coxph(Surv(time, event) ~ Group+Age+Sex+Stage, data = cox.info.gse17538)
summary(mul.cox)
#基础森林图
ggforest(mul.cox,  #coxph得到的Cox回归结果
         data = cox.info.gse17538,  #数据集
         main = 'Hazard ratio of GSE17538',  #标题
         cpositions = c(0.05, 0.15, 0.35),  #前三列距离
         fontsize = 1, #字体大小
         refLabel = 'reference', #相对变量的数值标签，也可改为1
         noDigits = 3 #保留HR值以及95%CI的小数位数
)


#######################################################################3
##Nanostring
{
  GSE39582_exp_matrix <- t(scale(t(GSE39582_exp_matrix[,sample_5FU])))
  GSE39582.RFS <- GSE39582_exp_matrix[hub_genes,]
  #prognosis model
  time <- clinical_infor_df[sample_5FU,]$time
  event <- clinical_infor_df[sample_5FU,]$event
  train.matrix <- na.omit(data.frame(time = time, event = event, t(GSE39582.RFS), stringsAsFactors = F))
  
  set.seed(1244746)
  grow.RFS=rfsrc(Surv(time,event) ~.,data =train.matrix,ntree = 1000,nsplit = 2)
  grow.RFS
  
  {
    GSE39582_rfs_score = grow.RFS$predicted
    summary(GSE39582_rfs_score)
    GSE39582_rfs_score.surv = cbind(train.matrix[,1:2],GSE39582_rfs_score)
    head(GSE39582_rfs_score.surv)
    #Determining the optimal cutpoint for GSE39582_rfs_score
    GSE39582_rfs_score.cut = surv_cutpoint(
      data = GSE39582_rfs_score.surv,
      time = "time",
      event = "event",
      variables = "GSE39582_rfs_score",
      minprop = 0.1
    )
    summary(GSE39582_rfs_score.cut)
    cut_off = GSE39582_rfs_score.cut$cutpoint$cutpoint
  }
}

library(YuGene)
##Nanostring
{
  #########################new-nanostring testing cohort
  library(openxlsx)
  #nanostring表达数据
  load("0.New Nanostring/去除批次效应/nanostring.mRNA.norm.RData")
  #nanostring临床数据
  nanostring_cli <- read.xlsx("0.New Nanostring/VUMCdataset_6.5.19_clinical_colon.xlsx",sheet = 6,colNames = T,rowNames = F,check.names=F)
  rownames(nanostring_cli) <- nanostring_cli$`NanoString.V3.#`
  nanostring_cli <- nanostring_cli[,-1]
  #nanostring 2期样本
  nanostring_2 <- nanostring_cli[grep("2",nanostring_cli$Stage.AJCC.2002),]
  #nanostring 23期样本
  nanostring_23 <- nanostring_cli[grep("2|3",nanostring_cli$Stage.AJCC.2002),]
  nanostring_23_cli <- nanostring_23[,c(4,3)]
  colnames(nanostring_23_cli) <- c("time","event")
  nanostring_23_cli$time <- nanostring_23_cli$time/(365/12)
  nanostring_23_cli$event[grep("N/A",nanostring_23_cli$event)] <- 0
  nanostring_23_cli$event[which(nanostring_23_cli$event != 0)] <- 1
  nanostring_23_cli$event <- as.numeric(nanostring_23_cli$event)
  
  #处理表达数据
  nanostring_23_exp <- combat_nanostring[,rownames(nanostring_23_cli)]
  nanostring_genes <- unique(sub("_.*","",rownames(nanostring_23_exp)))
  nanostring.exp.process <- function(nanostring.exp=NULL,method=c("average","median","min","max","MAD")){
    unique_exp <- t(apply(matrix(nanostring_genes,ncol = 1), 1, 
                          function(x){
                            gene_exp <- nanostring.exp[grep(x,rownames(nanostring.exp)),]
                            if(method == "average"){
                              return(colMeans(gene_exp))
                            }else if(method == "median"){
                              return(apply(gene_exp, 2, median))
                            }else if(method == "max"){
                              return(apply(gene_exp, 2, max))
                            }else if(method == "min"){
                              return(apply(gene_exp, 2, min))
                            }else if(method == "MAD"){
                              probe_var <- apply(gene_exp, 1, mad)
                              return(gene_exp[which.max(probe_var),])
                            }
                          }))
    unique_exp <- as.data.frame(unique_exp)
    rownames(unique_exp) <- nanostring_genes
    return(unique_exp)
  }
  
  ##1-2-3最大值表达
  nanostring_23_max <- nanostring.exp.process(nanostring.exp = nanostring_23_exp,
                                              method = "max")
}

#nanostring 用药样本
nanostring_23_drugs_all <- read.xlsx("0.New Nanostring/VUMCdataset_6.5.19_clinical_colon.xlsx",sheet = 5,colNames = T,rowNames = F,check.names=F)
rownames(nanostring_23_drugs_all) <- nanostring_23_drugs_all$`NanoString.V3.#`
nanostring_23_drugs <- read.xlsx("0.New Nanostring/VUMCdataset_6.5.19_clinical_colon.xlsx",sheet = 3,colNames = T,rowNames = F,check.names=F)
nanostring_23_drugs2 <- read.xlsx("0.New Nanostring/VUMCdataset_6.5.19_clinical_colon.xlsx",sheet = 4,colNames = T,rowNames = F,check.names=F)
test_cli <- nanostring_23_cli[nanostring_23_drugs$`NanoString.V3.#`,]
{
  ##5-years
  test_cli[which(test_cli[,1]>60),2] = 0
  test_cli[which(test_cli[,1]>60),1] = 60
}
test_matrix <- nanostring_23_max[,nanostring_23_drugs$`NanoString.V3.#`]
if(T){
  ##
  #test_matrix_z <- YuGene(test_matrix)
  test_matrix_z <- t(scale(t(test_matrix)))
  test <- data.frame(t(test_matrix_z[hub_genes,]),stringsAsFactors = F)
  
  pred <- predict(grow.RFS, newdata=test)
  score <- pred$predicted
  group <- rep(0, length(score))
  for (i in 1:length(group)) {
    if( score[i]<=cut_off) {
      group[i]<-"low"
    } else{
      group[i]<-"high"
    }
  }
  surv.df <- data.frame(test_cli[,1], test_cli[,2], score,stringsAsFactors = F)
  colnames(surv.df)[1:2] <- c("surtime", "censor")
  test.survdiff <- survdiff(Surv(surtime,censor)~group,data=surv.df)
  HR_temp <- (test.survdiff$obs[2]/test.survdiff$exp[2])/(test.survdiff$obs[1]/test.survdiff$exp[1])
  HR_temp <- round(HR_temp,4)
  p.val_temp <- round(1 - pchisq(test.survdiff$chisq, 1),4)
  up95_temp <- exp(log(as.numeric(HR_temp)) + qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
  up95_temp <- round(up95_temp,4)
  low95_temp <- exp(log(as.numeric(HR_temp)) - qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
  low95_temp <- round(low95_temp,4)
  c_index_temp <- rcorr.cens(-surv.df$score, Surv(surv.df$surtime, surv.df$censor))[1]
  c_index_temp <- round(c_index_temp,4)
  plot(survfit(Surv(surtime,censor)~group,data=surv.df),xlab="Survival in months", ylab="Survival probability",
       col=c("red","darkblue"),bty="l",lwd=4)
  title(main=paste("Nanostring", "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
  legend("topright", legend=c(paste0("high (N=",table(group)["high"],")"),paste0("low (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
}

Age <- nanostring_23[rownames(test_cli),]$Age.at.DX
Stage <- nanostring_23[rownames(test_cli),]$Stage.AJCC.2002
Stage[grep("2",Stage)] <- 2
Stage[grep("3",Stage)] <- 3
Stage <- factor(Stage,levels = c(2,3))
Group <- factor(group,levels=c("low","high"))
cox.info <- cbind(test_cli,Group) %>% cbind(Stage) %>% cbind(Age) %>% cbind(score) 

{
  MSI <- nanostring_23_drugs_all[rownames(cox.info),"MSI.Status"]
  cox.info$MSI <- factor(MSI, levels=c("MSS","MSI-H"))
  KRAS <- nanostring_23_drugs_all[rownames(cox.info),"KRAS"]
  KRAS[grep("p.",KRAS)] = "MT"
  cox.info$KRAS <- factor(KRAS, levels=c("wt","MT"))
  BRAF <- nanostring_23_drugs_all[rownames(cox.info),"BRAF"]
  BRAF[grep("p.",BRAF)] = "MT"
  cox.info$BRAF <- factor(BRAF, levels=c("wt","MT"))
  PIK3CA <- nanostring_23_drugs_all[rownames(cox.info),"PIK3CA"]
  PIK3CA[grep("p.",PIK3CA)] = "MT"
  cox.info$PIK3CA <- factor(PIK3CA, levels=c("wt","MT"))
  PTEN <- nanostring_23_drugs_all[rownames(cox.info),"PTEN"]
  PTEN[grep("p.",PTEN)] = "MT"
  cox.info$PTEN <- factor(PTEN, levels=c("wt","MT"))
  NRAS <- nanostring_23_drugs_all[rownames(cox.info),"NRAS"]
  NRAS[grep("p.",NRAS)] = "MT"
  cox.info$NRAS <- factor(NRAS, levels=c("wt","MT"))
  AKT1 <- nanostring_23_drugs_all[rownames(cox.info),"AKT1"]
  AKT1[grep("p.",AKT1)] = "MT"
  cox.info$AKT1 <- factor(AKT1, levels=c("wt","MT"))
  
}

## group
uni.cox.group <- coxph(Surv(time, event) ~ Group, data = cox.info)
summary(uni.cox.group)
## score
uni.cox.score <- coxph(Surv(time, event) ~ score, data = cox.info)
summary(uni.cox.score)
## age
uni.cox.age <- coxph(Surv(time, event) ~ Age, data = cox.info)
summary(uni.cox.age)
## stage
uni.cox.stage <- coxph(Surv(time, event) ~ Stage, data = cox.info)
summary(uni.cox.stage)

## MSI
uni.cox.MSI <- coxph(Surv(time, event) ~ MSI, data = cox.info)
summary(uni.cox.MSI)
## KRAS
uni.cox.KRAS <- coxph(Surv(time, event) ~ KRAS, data = cox.info)
summary(uni.cox.KRAS)
## BRAF
uni.cox.BRAF <- coxph(Surv(time, event) ~ BRAF, data = cox.info)
summary(uni.cox.BRAF)
## PIK3CA
uni.cox.PIK3CA <- coxph(Surv(time, event) ~ PIK3CA, data = cox.info)
summary(uni.cox.PIK3CA)

##
mul.cox <- coxph(Surv(time, event) ~ Group+Age+Stage+MSI, data = cox.info)
summary(mul.cox)
#基础森林图
ggforest(mul.cox,  #coxph得到的Cox回归结果
         data = cox.info,  #数据集
         main = 'Hazard ratio of Nanostring',  #标题
         cpositions = c(0.05, 0.15, 0.35),  #前三列距离
         fontsize = 1, #字体大小
         refLabel = 'reference', #相对变量的数值标签，也可改为1
         noDigits = 3 #保留HR值以及95%CI的小数位数
)


