##validated this prognosis model in the other seven independent cohorts with RFS time
c("GSE17536", "GSE17537", "GSE17538","GSE29621","GSE33113","GSE37892","GSE38832")
library(dplyr)
rm(list = ls())
setwd("H:/210E盘文档/colon cancer/0DataAndScript")

load("1.primary data/train/ESP.RData")
load("1.primary data/train/hub_probe.RData")
##construct RSF prognosis model according to the RFS time of GSE39582
library(randomForestSRC)
library(survival)
library(Hmisc)
hub_genes <- as.character(hub_probe$gene[1:18])
GSE39582.RFS <- GSE39582_exp_matrix[hub_genes,]
#prognosis model
time <- clinical_infor_df$time
event <- clinical_infor_df$event
train.matrix <- data.frame(time = time, event = event, t(GSE39582.RFS), stringsAsFactors = F)
set.seed(1244746)
grow.RFS=rfsrc(Surv(time,event) ~.,data =train.matrix,ntree = 1000,nsplit = 2)
grow.RFS

library(survminer)
#######use the model.grow$predicted, which is in-bag predicted value, 
#######to apply maxstat, then use the cutoff value on testing data;
GSE39582_rfs_score = grow.RFS$predicted
summary(GSE39582_rfs_score)
GSE39582_rfs_score.surv = cbind(na.omit(train.matrix[,1:2]),GSE39582_rfs_score)
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
#Categorize GSE39582_rfs_score variable
GSE39582_rfs_score.cat = surv_categorize(GSE39582_rfs_score.cut)

all(rownames(GSE39582_rfs_score.surv) == rownames(GSE39582_rfs_score.cat))
GSE39582.info = cbind(GSE39582_rfs_score.cat, GSE39582_rfs_score.surv$GSE39582_rfs_score)
colnames(GSE39582.info)[3:4] = c("group","score")
save(GSE39582_exp_matrix, GSE39582.info, file = "Script/2022.11.03/GSE39582.18gene.RData")


#############GSE17536####################################
gse.name <- "GSE17536"
test_data=paste0(gse.name,".RData")
c_index=c()
p_val=c()
HR=c()
up95=c()
low95=c()
for (k in 1:length(test_data)) {
  load(paste0("1.primary data/test/",test_data[k]))
  test_matrix <- test_matrix[probe2gene_loc,]
  test_matrix <- test_matrix[sig_probe_loc,]
  rownames(test_matrix) <- all.gene #assign gene symbols to matrix
  test <- data.frame(t(test_matrix[hub_genes,]),stringsAsFactors = F)
  #predict testing set
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
  title(main=paste(test_data[k], "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
  legend("topright", legend=c("high","low"), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  c_index=c(c_index,c_index_temp)
  p_val=c(p_val,p.val_temp)
  HR=c(HR,HR_temp)
  up95=c(up95,up95_temp)
  low95=c(low95,low95_temp)
}

GSE17536_clinical <- read.table(paste0("1.primary data/test/",gse.name,".frma.sample.tsv"),
                           header = T,sep = "\t",stringsAsFactors = F)
GSE17536_mid_drug<-read.table("4.ACT/0.data/GSE17536_mid_drug.txt",sep="\t",header=T,stringsAsFactors = F)
GSE17536_sample_ACT = GSE17536_mid_drug$geoid
rownames(GSE17536_clinical) <- GSE17536_clinical$Sample
Age <- GSE17536_clinical[rownames(test_cli),"Age"]
Sex <- GSE17536_clinical[rownames(test_cli),"Sex"]
Stage <- GSE17536_clinical[rownames(test_cli),"Stage"]
Group <- factor(group,levels=c("low","high"))

cox.info.gse17536 <- cbind(test_cli,Age) %>% cbind(Sex) %>% cbind(Group) %>% cbind(Stage)
cox.info.gse17536$ACT = "No"
cox.info.gse17536[intersect(rownames(cox.info.gse17536),GSE17536_sample_ACT),7] = "Yes"
## group
uni.cox.group <- coxph(Surv(time, event) ~ Group, data = cox.info.gse17536)
summary(uni.cox.group)
## score
uni.cox.score <- coxph(Surv(time, event) ~ score, data = cox.info.gse17536)
summary(uni.cox.score)
## age
uni.cox.age <- coxph(Surv(time, event) ~ Age, data = cox.info.gse17536)
summary(uni.cox.age)
## sex
uni.cox.sex <- coxph(Surv(time, event) ~ Sex, data = cox.info.gse17536)
summary(uni.cox.sex)
##stage
cox.info.gse17536$Stage <- factor(cox.info.gse17536$Stage,levels=c(2,3))
uni.cox.stage <- coxph(Surv(time, event) ~ Stage, data = cox.info.gse17536)
summary(uni.cox.stage)
##ACT
cox.info.gse17536$ACT <- factor(cox.info.gse17536$ACT,levels=c("No","Yes"))
uni.cox.act <- coxph(Surv(time, event) ~ ACT, data = cox.info.gse17536)
summary(uni.cox.act)
## age+sex+group+stage
mul.cox <- coxph(Surv(time, event) ~ Group+Age+Sex+Stage+ACT, data = cox.info.gse17536)
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



#############GSE17537#####
gse.name <- "GSE17537"
test_data=paste0(gse.name,".RData")
c_index=c()
p_val=c()
HR=c()
up95=c()
low95=c()
for (k in 1:length(test_data)) {
  load(paste0("1.primary data/test/",test_data[k]))
  test_matrix <- test_matrix[probe2gene_loc,]
  test_matrix <- test_matrix[sig_probe_loc,]
  rownames(test_matrix) <- all.gene #assign gene symbols to matrix
  test <- data.frame(t(test_matrix[hub_genes,]),stringsAsFactors = F)
  #predict testing set
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
  title(main=paste(test_data[k], "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
  legend("topright", legend=c("high","low"), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  c_index=c(c_index,c_index_temp)
  p_val=c(p_val,p.val_temp)
  HR=c(HR,HR_temp)
  up95=c(up95,up95_temp)
  low95=c(low95,low95_temp)
}

GSE17537_clinical <- read.table(paste0("1.primary data/test/",gse.name,".frma.sample.tsv"),
                                header = T,sep = "\t",stringsAsFactors = F)
#drug related information##########提取GSE17537中期使用过药物的样本
load("4.ACT/0.data/GSE17537_matrix.Rdata")
GSE17537_sample_ACT <- GSE17537_mid_drug$geoid
rownames(GSE17537_clinical) <- GSE17537_clinical$Sample
Age <- GSE17537_clinical[rownames(test_cli),"Age"]
Sex <- GSE17537_clinical[rownames(test_cli),"Sex"]
Group <- factor(group,levels=c("low","high"))
Stage <- factor(GSE17537_clinical[rownames(test_cli),"Stage"],levels = c(2,3))

cox.info.gse17537 <- cbind(test_cli,Age) %>% cbind(Sex) %>% cbind(Group) %>% cbind(Stage)
cox.info.gse17537$ACT = "No"
cox.info.gse17537[intersect(rownames(cox.info.gse17537),GSE17537_sample_ACT),7] = "Yes"
## group
uni.cox.group <- coxph(Surv(time, event) ~ Group, data = cox.info.gse17537)
summary(uni.cox.group)
## score
uni.cox.score <- coxph(Surv(time, event) ~ score, data = cox.info.gse17537)
summary(uni.cox.score)
## age
uni.cox.age <- coxph(Surv(time, event) ~ Age, data = cox.info.gse17537)
summary(uni.cox.age)
## sex
uni.cox.sex <- coxph(Surv(time, event) ~ Sex, data = cox.info.gse17537)
summary(uni.cox.sex)
##stage
uni.cox.stage <- coxph(Surv(time, event) ~ Stage, data = cox.info.gse17537)
summary(uni.cox.stage)
##ACT
cox.info.gse17537$ACT <- factor(cox.info.gse17537$ACT,levels=c("No","Yes"))
uni.cox.act <- coxph(Surv(time, event) ~ ACT, data = cox.info.gse17537)
summary(uni.cox.act)

## age+sex+group
mul.cox <- coxph(Surv(time, event) ~ Group+Age+Sex+Stage+ACT, data = cox.info.gse17537)
summary(mul.cox)
#基础森林图
ggforest(mul.cox,  #coxph得到的Cox回归结果
         data = cox.info.gse17537,  #数据集
         main = 'Hazard ratio of GSE17537',  #标题
         cpositions = c(0.05, 0.15, 0.35),  #前三列距离
         fontsize = 1, #字体大小
         refLabel = 'reference', #相对变量的数值标签，也可改为1
         noDigits = 3 #保留HR值以及95%CI的小数位数
)



#############GSE17538
gse.name <- "GSE17538"
test_data=paste0(gse.name,".RData")
c_index=c()
p_val=c()
HR=c()
up95=c()
low95=c()
for (k in 1:length(test_data)) {
  load(paste0("1.primary data/test/",test_data[k]))
  test_matrix <- test_matrix[probe2gene_loc,]
  test_matrix <- test_matrix[sig_probe_loc,]
  rownames(test_matrix) <- all.gene #assign gene symbols to matrix
  test <- data.frame(t(test_matrix[hub_genes,]),stringsAsFactors = F)
  #predict testing set
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
  title(main=paste(test_data[k], "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
  legend("topright", legend=c("high","low"), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  c_index=c(c_index,c_index_temp)
  p_val=c(p_val,p.val_temp)
  HR=c(HR,HR_temp)
  up95=c(up95,up95_temp)
  low95=c(low95,low95_temp)
}

Group <- factor(group,levels=c("low","high"))
Score <- score
cox.info.gse17538 <- rbind(cox.info.gse17536,cox.info.gse17537)
cox.info.gse17538$Group <- Group
cox.info.gse17538$Score <- Score
# save(test_matrix,cox.info.gse17538, file = "Script/2022.11.01/GSE17538.18gene.RData")

## group
uni.cox.group <- coxph(Surv(time, event) ~ Group, data = cox.info.gse17538)
summary(uni.cox.group)
## score
uni.cox.score <- coxph(Surv(time, event) ~ Score, data = cox.info.gse17538)
summary(uni.cox.score)
## age
uni.cox.age <- coxph(Surv(time, event) ~ Age, data = cox.info.gse17538)
summary(uni.cox.age)
## sex
uni.cox.sex <- coxph(Surv(time, event) ~ Sex, data = cox.info.gse17538)
summary(uni.cox.sex)
##ACT
cox.info.gse17538$ACT <- factor(cox.info.gse17538$ACT,levels=c("No","Yes"))
uni.cox.act <- coxph(Surv(time, event) ~ ACT, data = cox.info.gse17538)
summary(uni.cox.act)
##
uni.cox.stage <- coxph(Surv(time, event) ~ Stage, data = cox.info.gse17538)
summary(uni.cox.stage)

## age+sex+group
mul.cox <- coxph(Surv(time, event) ~ Group+Age+Sex+Stage+ACT, data = cox.info.gse17538)
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


#####GSE29621
gse.name <- "GSE29621"
test_data=paste0(gse.name,".RData")
c_index=c()
p_val=c()
HR=c()
up95=c()
low95=c()
for (k in 1:length(test_data)) {
  load(paste0("1.primary data/test/",test_data[k]))
  test_matrix <- test_matrix[probe2gene_loc,]
  test_matrix <- test_matrix[sig_probe_loc,]
  rownames(test_matrix) <- all.gene #assign gene symbols to matrix
  test <- data.frame(t(test_matrix[hub_genes,]),stringsAsFactors = F)
  #predict testing set
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
  title(main=paste(test_data[k], "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
  legend("topright", legend=c("high","low"), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  c_index=c(c_index,c_index_temp)
  p_val=c(p_val,p.val_temp)
  HR=c(HR,HR_temp)
  up95=c(up95,up95_temp)
  low95=c(low95,low95_temp)
}

clinal_info <- read.csv("1.primary data/test/GSE29621.sample.txt",header = T,sep = "\t",stringsAsFactors = F)
stageII_III <- which(1<clinal_info$"stage" & clinal_info$"stage"<4)
clinal_info <- clinal_info[stageII_III,c(1:4,6:7)]
rownames(clinal_info) <- clinal_info$sample
clinal_info <- clinal_info[,-1]

clinal_info <- clinal_info[rownames(test_cli),]
clinal_info <- cbind(clinal_info,group)

clinal_info$group <- factor(clinal_info$group,levels=c("low","high"))
clinal_info$stage <- factor(clinal_info$stage,levels=c(2,3))
cox.info <- clinal_info

colnames(cox.info) = c("Stage","time","event","Age","Sex","Group")
## group
uni.cox.group <- coxph(Surv(time, event) ~ Group, data = cox.info)
summary(uni.cox.group)
## score
uni.cox.score <- coxph(Surv(time, event) ~ score, data = cox.info)
summary(uni.cox.score)
## age
uni.cox.age <- coxph(Surv(time, event) ~ Age, data = cox.info)
summary(uni.cox.age)
## sex
uni.cox.sex <- coxph(Surv(time, event) ~ Sex, data = cox.info)
summary(uni.cox.sex)
## stage
uni.cox.stage <- coxph(Surv(time, event) ~ Stage, data = cox.info)
summary(uni.cox.stage)

## age+sex+group
mul.cox <- coxph(Surv(time, event) ~ Group+Age+Sex+Stage, data = cox.info)
summary(mul.cox)
#基础森林图
ggforest(mul.cox,  #coxph得到的Cox回归结果
         data = cox.info,  #数据集
         main = 'Hazard ratio of GSE29621',  #标题
         cpositions = c(0.05, 0.15, 0.35),  #前三列距离
         fontsize = 1, #字体大小
         refLabel = 'reference', #相对变量的数值标签，也可改为1
         noDigits = 3 #保留HR值以及95%CI的小数位数
)


#####GSE33113
gse.name <- "GSE33113"
test_data=paste0(gse.name,".RData")
c_index=c()
p_val=c()
HR=c()
up95=c()
low95=c()
for (k in 1:length(test_data)) {
  load(paste0("1.primary data/test/",test_data[k]))
  test_matrix <- test_matrix[probe2gene_loc,]
  test_matrix <- test_matrix[sig_probe_loc,]
  rownames(test_matrix) <- all.gene #assign gene symbols to matrix
  test <- data.frame(t(test_matrix[hub_genes,]),stringsAsFactors = F)
  #predict testing set
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
  title(main=paste(test_data[k], "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
  legend("topright", legend=c("high","low"), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  c_index=c(c_index,c_index_temp)
  p_val=c(p_val,p.val_temp)
  HR=c(HR,HR_temp)
  up95=c(up95,up95_temp)
  low95=c(low95,low95_temp)
}

GSE_clinical <- read.table(paste0("1.primary data/test/",gse.name,".frma.sample.tsv"),
                           header = T,sep = "\t",stringsAsFactors = F)
rownames(GSE_clinical) <- GSE_clinical$Sample
Age <- as.numeric(GSE_clinical[rownames(test_cli),"Age"])
Sex <- factor(GSE_clinical[rownames(test_cli),"Sex"], levels=c("f","m"))
Group <- factor(group,levels=c("low","high"))
cox.info <- cbind(test_cli,Age) %>% cbind(Sex) %>% cbind(Group) %>% cbind(score)
## group
uni.cox.group <- coxph(Surv(time, event) ~ Group, data = cox.info)
summary(uni.cox.group)
# score
uni.cox.score <- coxph(Surv(time, event) ~ score, data = cox.info)
summary(uni.cox.score)
## age
uni.cox.age <- coxph(Surv(time, event) ~ Age, data = cox.info)
summary(uni.cox.age)
## sex
uni.cox.sex <- coxph(Surv(time, event) ~ Sex, data = cox.info)
summary(uni.cox.sex)
## age+sex+group
mul.cox <- coxph(Surv(time, event) ~ Group+Age+Sex, data = cox.info)
summary(mul.cox)
#基础森林图
ggforest(mul.cox,  #coxph得到的Cox回归结果
         data = cox.info,  #数据集
         main = 'Hazard ratio of GSE33113',  #标题
         cpositions = c(0.05, 0.15, 0.35),  #前三列距离
         fontsize = 1, #字体大小
         refLabel = 'reference', #相对变量的数值标签，也可改为1
         noDigits = 3 #保留HR值以及95%CI的小数位数
)


#####GSE37892
gse.name <- "GSE37892"
test_data=paste0(gse.name,".RData")
c_index=c()
p_val=c()
HR=c()
up95=c()
low95=c()
for (k in 1:length(test_data)) {
  load(paste0("1.primary data/test/",test_data[k]))
  test_matrix <- test_matrix[probe2gene_loc,]
  test_matrix <- test_matrix[sig_probe_loc,]
  rownames(test_matrix) <- all.gene #assign gene symbols to matrix
  test <- data.frame(t(test_matrix[hub_genes,]),stringsAsFactors = F)
  #predict testing set
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
  title(main=paste(test_data[k], "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
  legend("topright", legend=c("high","low"), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  c_index=c(c_index,c_index_temp)
  p_val=c(p_val,p.val_temp)
  HR=c(HR,HR_temp)
  up95=c(up95,up95_temp)
  low95=c(low95,low95_temp)
}

GSE_clinical <- read.table(paste0("1.primary data/test/",gse.name,".frma.sample.tsv"),
                           header = T,sep = ",",stringsAsFactors = F)
rownames(GSE_clinical) <- GSE_clinical$sample
Age <- as.numeric(GSE_clinical[rownames(test_cli),"age"])
Sex <- GSE_clinical[rownames(test_cli),"gender"]
Stage <- factor(GSE_clinical[rownames(test_cli),"stage"],levels = c(2,3))
Group <- factor(group,levels=c("low","high"))
cox.info <- cbind(test_cli,Age) %>% cbind(Sex) %>% cbind(Group) %>% cbind(Stage) %>% cbind(score)
## group
uni.cox.group <- coxph(Surv(time, event) ~ Group, data = cox.info)
summary(uni.cox.group)
## score
uni.cox.score <- coxph(Surv(time, event) ~ score, data = cox.info)
summary(uni.cox.score)
## age
uni.cox.age <- coxph(Surv(time, event) ~ Age, data = cox.info)
summary(uni.cox.age)
## sex
uni.cox.sex <- coxph(Surv(time, event) ~ Sex, data = cox.info)
summary(uni.cox.sex)
## stage
uni.cox.stage <- coxph(Surv(time, event) ~ Stage, data = cox.info)
summary(uni.cox.stage)

## age+sex+group
mul.cox <- coxph(Surv(time, event) ~ Group+Age+Sex+Stage, data = cox.info)
summary(mul.cox)
#基础森林图
ggforest(mul.cox,  #coxph得到的Cox回归结果
         data = cox.info,  #数据集
         main = 'Hazard ratio of GSE37892',  #标题
         cpositions = c(0.05, 0.15, 0.35),  #前三列距离
         fontsize = 1, #字体大小
         refLabel = 'reference', #相对变量的数值标签，也可改为1
         noDigits = 3 #保留HR值以及95%CI的小数位数
)


#####GSE38832
gse.name <- "GSE38832"
test_data=paste0(gse.name,".RData")
c_index=c()
p_val=c()
HR=c()
up95=c()
low95=c()
for (k in 1:length(test_data)) {
  load(paste0("1.primary data/test/",test_data[k]))
  test_matrix <- test_matrix[probe2gene_loc,]
  test_matrix <- test_matrix[sig_probe_loc,]
  rownames(test_matrix) <- all.gene #assign gene symbols to matrix
  test <- data.frame(t(test_matrix[hub_genes,]),stringsAsFactors = F)
  #predict testing set
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
  title(main=paste(test_data[k], "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
  legend("topright", legend=c("high","low"), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  c_index=c(c_index,c_index_temp)
  p_val=c(p_val,p.val_temp)
  HR=c(HR,HR_temp)
  up95=c(up95,up95_temp)
  low95=c(low95,low95_temp)
}

GSE_clinical <- read.table(paste0("1.primary data/test/",gse.name,".sample.txt"),
                           header = T,sep = "\t",stringsAsFactors = F)
rownames(GSE_clinical) <- GSE_clinical$Sample
Stage <- factor(GSE_clinical[rownames(test_cli),"Stage"],levels = c(2,3))
Group <- factor(group,levels=c("low","high"))
cox.info <- cbind(test_cli,Group) %>% cbind(Stage) %>% cbind(score)
## group
uni.cox.group <- coxph(Surv(cox.info$DFSurvivalTime,cox.info$DFSurvivalEvent) ~ Group, data = cox.info)
summary(uni.cox.group)
## score
uni.cox.score <- coxph(Surv(cox.info$DFSurvivalTime,cox.info$DFSurvivalEvent) ~ score, data = cox.info)
summary(uni.cox.score)
# stgae
uni.cox.stgae <- coxph(Surv(cox.info$DFSurvivalTime,cox.info$DFSurvivalEvent) ~ Stage, data = cox.info)
summary(uni.cox.stgae)

mul.cox <- coxph(Surv(cox.info$DFSurvivalTime,cox.info$DFSurvivalEvent) ~ Group+Stage, data = cox.info)
summary(mul.cox)
#基础森林图
ggforest(mul.cox,  #coxph得到的Cox回归结果
         data = cox.info,  #数据集
         main = 'Hazard ratio of GSE38832',  #标题
         cpositions = c(0.05, 0.15, 0.35),  #前三列距离
         fontsize = 1, #字体大小
         refLabel = 'reference', #相对变量的数值标签，也可改为1
         noDigits = 3 #保留HR值以及95%CI的小数位数
)


####
rm(list = ls())
library(dplyr)
setwd("H:/210E盘文档/colon cancer/0DataAndScript")
load("1.primary data/train/ESP.RData")
load("1.primary data/train/hub_probe.RData")

##construct RSF prognosis model according to the RFS time of GSE39582
library(randomForestSRC)
library(survival)
library(Hmisc)
library(survminer)
library(maxstat)

{
  hub_genes <- as.character(hub_probe$gene[1:18])
  
  GSE39582_exp_matrix <- t(scale(t(GSE39582_exp_matrix)))
  GSE39582.RFS <- GSE39582_exp_matrix[hub_genes,]
  #prognosis model
  time <- clinical_infor_df$time
  event <- clinical_infor_df$event
  train.matrix <- data.frame(time = time, event = event, t(GSE39582.RFS), stringsAsFactors = F)
  
  set.seed(1244746)
  grow.RFS=rfsrc(Surv(time,event) ~.,data =train.matrix,ntree = 1000,nsplit = 2)
  grow.RFS
  
  {
    GSE39582_rfs_score = grow.RFS$predicted
    summary(GSE39582_rfs_score)
    GSE39582_rfs_score.surv = cbind(na.omit(train.matrix[,1:2]),GSE39582_rfs_score)
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

##Nanostring
{
  #########################new-nanostring testing cohort
  library(openxlsx)
  #nanostring表达数据
  load("0.New Nanostring/去除批次效应/nanostring.mRNA.norm.RData")
  #nanostring临床数据
  nanostring_cli <- read.xlsx("0.New Nanostring/VUMCdataset_6.5.19_clinical_colon.xlsx",sheet = 2,colNames = T,rowNames = F,check.names=F)
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

test_cli <- nanostring_23_cli
{
  ##5-years
  test_cli[which(test_cli[,1]>60),2] = 0
  test_cli[which(test_cli[,1]>60),1] = 60
}
test_matrix <- nanostring_23_max
{
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
  legend("topright", legend=c("high","low"), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
}

Age <- nanostring_23$Age.at.DX
Stage <- nanostring_23$Stage.AJCC.2002
Stage[grep("2",Stage)] <- 2
Stage[grep("3",Stage)] <- 3
Stage <- factor(Stage,levels = c(2,3))
Group <- factor(group,levels=c("low","high"))
cox.info <- cbind(test_cli,Group) %>% cbind(Stage) %>% cbind(Age) %>% cbind(score)

#nanostring 用药样本
nanostring_23_drugs_all <- read.xlsx("0.New Nanostring/VUMCdataset_6.5.19_clinical_colon.xlsx",sheet = 5,colNames = T,rowNames = F,check.names=F)
rownames(nanostring_23_drugs_all) <- nanostring_23_drugs_all$`NanoString.V3.#`
nanostring_23_drugs <- read.xlsx("0.New Nanostring/VUMCdataset_6.5.19_clinical_colon.xlsx",sheet = 3,colNames = T,rowNames = F,check.names=F)
nanostring_23_drugs2 <- read.xlsx("0.New Nanostring/VUMCdataset_6.5.19_clinical_colon.xlsx",sheet = 4,colNames = T,rowNames = F,check.names=F)
cox.info$ACT = "No"
cox.info[intersect(rownames(cox.info),nanostring_23_drugs$`NanoString.V3.#`),7] = "Yes"


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
##ACT
cox.info$ACT <- factor(cox.info$ACT,levels=c("No","Yes"))
uni.cox.act <- coxph(Surv(time, event) ~ ACT, data = cox.info)
summary(uni.cox.act)

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
mul.cox <- coxph(Surv(time, event) ~ Group+Age+Stage+ACT+MSI, data = cox.info)
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


################################################################## not run