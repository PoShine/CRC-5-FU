rm(list = ls()); gc()
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
#hub_genes <- as.character(hub_probe$gene[1:24])[c(1,2,4,7:11,13,15,16,18,19,23)]
#hub_genes <- as.character(hub_probe$gene[1:6])
#select the sample of training set
hub_genes <- as.character(hub_probe$gene[1:18])

#hub_genes <- hub_genes[c(1:2,4:6,8:10,13,15:16)]
hub_genes
###
hub_genes_combn <- combn(hub_genes,10)
hub_genes <- hub_genes_combn[,10651]
hub_genes

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

##################GSE14333  not run
{
  load("1.primary data/test/GSE14333.RData")
  test_matrix <- test_matrix[probe2gene_loc,]
  test_matrix <- test_matrix[sig_probe_loc,]
  rownames(test_matrix) <- all.gene #assign gene symbols to matrix
  
  gse14333_drug = read.table("1.primary data/test/GSE14333_mid_drug.txt",
                             header = T, sep = "\t",stringsAsFactors = F)
  gse1433_withACT = gse14333_drug$GEOID
  gse1433_withoutACT = setdiff(colnames(test_matrix),gse1433_withACT)
  
  test_cli = test_cli[gse1433_withACT,]
  
  test <- data.frame(t(test_matrix[hub_genes,gse1433_withACT]),stringsAsFactors = F)
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
  HR_temp <- (test.survdiff$obs[1]/test.survdiff$exp[1])/(test.survdiff$obs[2]/test.survdiff$exp[2])
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
  title(main="GSE14333",cex.main=1, adj=0.4)
  #length(grep("high",group))
  legend("topright", legend=c(paste0("ACT-nonbenefit (N=",table(group)["high"],")"),paste0("ACT-benefit (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
}

##################GSE17536 testing cohort
{
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
  HR_temp <- (test.survdiff$obs[1]/test.survdiff$exp[1])/(test.survdiff$obs[2]/test.survdiff$exp[2])
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
  legend("topright", legend=c(paste0("ACT-nonbenefit (N=",table(group)["high"],")"),paste0("ACT-benefit (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
}

{#high
  group.536 = readRDS("3.18genes/GSE17536.RDS")
  group.536.high = subset(group.536, group=="high")
  group.536.high.act = intersect(rownames(group.536.high), GSE17536_sample_ACT)
  test_cli_high <- test_cli[group.536.high.act,]
  #predict testing set
  test <- data.frame(t(test_matrix[hub_genes,group.536.high.act]),stringsAsFactors = F)
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
  
  surv.df <- data.frame(test_cli_high[,1], test_cli_high[,2], score,stringsAsFactors = F)
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
  legend("topright", legend=c(paste0("no response (N=",table(group)["high"],")"),paste0("response (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
}
{#low
  group.536 = readRDS("3.18genes/GSE17536.RDS")
  group.536.low = subset(group.536, group=="low")
  group.536.low.act = intersect(rownames(group.536.low), GSE17536_sample_ACT)
  test_cli_low <- test_cli[group.536.low.act,]
  #predict testing set
  test <- data.frame(t(test_matrix[hub_genes,group.536.low.act]),stringsAsFactors = F)
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
  
  surv.df <- data.frame(test_cli_low[,1], test_cli_low[,2], score,stringsAsFactors = F)
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
  legend("topright", legend=c(paste0("no response (N=",table(group)["high"],")"),paste0("response (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
}

library(dplyr)
###################GSE29621 testing cohort
{
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
  HR_temp <- (test.survdiff$obs[1]/test.survdiff$exp[1])/(test.survdiff$obs[2]/test.survdiff$exp[2])
  HR_temp <- round(HR_temp,3)
  p.val_temp <- round(1 - pchisq(test.survdiff$chisq, 1),4)
  up95_temp <- exp(log(as.numeric(HR_temp)) + qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
  up95_temp <- round(up95_temp,3)
  low95_temp <- exp(log(as.numeric(HR_temp)) - qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
  low95_temp <- round(low95_temp,3)
  c_index_temp <- rcorr.cens(-surv.df$score, Surv(surv.df$surtime, surv.df$censor))[1]
  c_index_temp <- round(c_index_temp,4)
  plot(survfit(Surv(surtime,censor)~group,data=surv.df),xlab="Survival in months", ylab="Survival probability",
       col=c("red","darkblue"),bty="l",lwd=4)
  title(main=paste("GSE29621", "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
  legend("topright", legend=c(paste0("ACT-nonbenefit (N=",table(group)["high"],")"),paste0("ACT-benefit (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
  #生存图
  fit = survfit(Surv(surtime,censor)~group,data=surv.df)
  surv.p <- ggsurvplot(fit = fit,
                       data = surv.df,
                       pval = TRUE, 
                       pval.coord = c(40,0.1),
                       legend = c(0.8,0.3),
                       legend.title = "Group",
                       legend.labs = c("ACT-nonbenefit", "ACT-benefit"),
                       palette = c("red","darkblue"),
                       title = paste("GSE29621", "(pvalue=", p.val_temp, ")", sep="" ),
                       xlab = "Survival in months",
                       risk.table = TRUE,       # show risk table.
                       tables.height = 0.2,
                       risk.table.y.text.col = T, # colour risk table text annotations.
                       risk.table.y.text = FALSE, # show bars instead of names in text annotations
                       tables.theme = theme_cleantable(),
                       linetype = 1,
                       surv.median.line = "hv",
                       ggtheme = theme_survminer()
  )
  surv.p$plot = surv.p$plot+annotate("text", x=0, y=0.1,
                                     label=paste0(paste0("P value:",p.val_temp),"\n",
                                                  paste0(paste0("HR:",HR_temp)," (95%CI,",low95_temp,"-",up95_temp,")")),
                                     size = 5,hjust = 0)
  pdf(paste0("Script/revise4.26/figure6_","GSE29621",".pdf"), width = 6, height = 7)
  print(surv.p)
  dev.off()
  
}
{#high
  group.536 = readRDS("3.18genes/GSE29621.RDS")
  group.536.high = subset(group.536, group=="high")
  group.536.high.act = intersect(rownames(group.536.high), GSE29621.sample.drug21)
  test_cli_high <- test_cli[group.536.high.act,]
  #predict testing set
  test <- data.frame(t(test_matrix[hub_genes,group.536.high.act]),stringsAsFactors = F)
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
  
  surv.df <- data.frame(test_cli_high[,1], test_cli_high[,2], score,stringsAsFactors = F)
  colnames(surv.df)[1:2] <- c("surtime", "censor")
  test.survdiff <- survdiff(Surv(surtime,censor)~group,data=surv.df)
  HR_temp <- (test.survdiff$obs[1]/test.survdiff$exp[1])/(test.survdiff$obs[2]/test.survdiff$exp[2])
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
  title(main="GSE29621(high-risk with ACT)",cex.main=1, adj=0.4)
  legend("topright", legend=c(paste0("no response (N=",table(group)["high"],")"),paste0("response (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
}
{#low
  group.536 = readRDS("3.18genes/GSE29621.RDS")
  group.536.low = subset(group.536, group=="low")
  group.536.low.act = intersect(rownames(group.536.low), GSE29621.sample.drug21)
  test_cli_low <- test_cli[group.536.low.act,]
  #predict testing set
  test <- data.frame(t(test_matrix[hub_genes,group.536.low.act]),stringsAsFactors = F)
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
  
  surv.df <- data.frame(test_cli_low[,1], test_cli_low[,2], score,stringsAsFactors = F)
  colnames(surv.df)[1:2] <- c("surtime", "censor")
  test.survdiff <- survdiff(Surv(surtime,censor)~group,data=surv.df)
  HR_temp <- (test.survdiff$obs[1]/test.survdiff$exp[1])/(test.survdiff$obs[2]/test.survdiff$exp[2])
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
  title(main="GSE29621(low-risk with ACT)",cex.main=1, adj=0.4)
  legend("topright", legend=c(paste0("no response (N=",table(group)["high"],")"),paste0("response (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
}

##################GSE17537 testing cohort
{
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
  legend("topright", legend=c(paste0("ACT-nonbenefit (N=",table(group)["high"],")"),paste0("ACT-benefit (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
}


##################GSE17538 testing cohort
{
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
  HR_temp <- (test.survdiff$obs[1]/test.survdiff$exp[1])/(test.survdiff$obs[2]/test.survdiff$exp[2])
  HR_temp <- round(HR_temp,3)
  p.val_temp <- round(1 - pchisq(test.survdiff$chisq, 1),4)
  up95_temp <- exp(log(as.numeric(HR_temp)) + qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
  up95_temp <- round(up95_temp,3)
  low95_temp <- exp(log(as.numeric(HR_temp)) - qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
  low95_temp <- round(low95_temp,3)
  c_index_temp <- rcorr.cens(-surv.df$score, Surv(surv.df$surtime, surv.df$censor))[1]
  c_index_temp <- round(c_index_temp,4)
  plot(survfit(Surv(surtime,censor)~group,data=surv.df),xlab="Survival in months", ylab="Survival probability",
       col=c("red","darkblue"),bty="l",lwd=4)
  title(main=paste("GSE17538", "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
  legend("topright", legend=c(paste0("ACT-nonbenefit (N=",table(group)["high"],")"),paste0("ACT-benefit (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
  #生存图
  fit = survfit(Surv(surtime,censor)~group,data=surv.df)
  surv.p <- ggsurvplot(fit = fit,
                       data = surv.df,
                       pval = TRUE, 
                       pval.coord = c(40,0.1),
                       legend = c(0.8,0.3),
                       legend.title = "Group",
                       legend.labs = c("ACT-nonbenefit", "ACT-benefit"),
                       palette = c("red","darkblue"),
                       title = paste("GSE17538", "(pvalue=", p.val_temp, ")", sep="" ),
                       xlab = "Survival in months",
                       risk.table = TRUE,       # show risk table.
                       tables.height = 0.2,
                       risk.table.y.text.col = T, # colour risk table text annotations.
                       risk.table.y.text = FALSE, # show bars instead of names in text annotations
                       tables.theme = theme_cleantable(),
                       linetype = 1,
                       surv.median.line = "hv",
                       ggtheme = theme_survminer()
  )
  surv.p$plot = surv.p$plot+annotate("text", x=0, y=0.1,
                                     label=paste0(paste0("P value:",p.val_temp),"\n",
                                                  paste0(paste0("HR:",HR_temp)," (95%CI,",low95_temp,"-",up95_temp,")")),
                                     size = 5,hjust = 0)
  pdf(paste0("Script/revise4.26/figure6_","GSE17538",".pdf"), width = 6, height = 7)
  print(surv.p)
  dev.off()
  
}

{#high
  group.538 = readRDS("3.18genes/GSE17538.RDS")
  group.538.high = subset(group.538, group=="high")
  group.538.high.act = intersect(rownames(group.538.high), GSE17538_sample_ACT)
  test_cli_high <- test_cli[group.538.high.act,]
  #predict testing set
  test <- data.frame(t(test_matrix[hub_genes,group.538.high.act]),stringsAsFactors = F)
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
  
  surv.df <- data.frame(test_cli_high[,1], test_cli_high[,2], score,stringsAsFactors = F)
  colnames(surv.df)[1:2] <- c("surtime", "censor")
  test.survdiff <- survdiff(Surv(surtime,censor)~group,data=surv.df)
  HR_temp <- (test.survdiff$obs[1]/test.survdiff$exp[1])/(test.survdiff$obs[2]/test.survdiff$exp[2])
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
  title(main="GSE17538(high-risk with ACT)",cex.main=1, adj=0.4)
  legend("topright", legend=c(paste0("no response (N=",table(group)["high"],")"),paste0("response (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
}
{#low
  group.538 = readRDS("3.18genes/GSE17538.RDS")
  group.538.low = subset(group.538, group=="low")
  group.538.low.act = intersect(rownames(group.538.low), GSE17538_sample_ACT)
  test_cli_low <- test_cli[group.538.low.act,]
  #predict testing set
  test <- data.frame(t(test_matrix[hub_genes,group.538.low.act]),stringsAsFactors = F)
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
  
  surv.df <- data.frame(test_cli_low[,1], test_cli_low[,2], score,stringsAsFactors = F)
  colnames(surv.df)[1:2] <- c("surtime", "censor")
  test.survdiff <- survdiff(Surv(surtime,censor)~group,data=surv.df)
  HR_temp <- (test.survdiff$obs[1]/test.survdiff$exp[1])/(test.survdiff$obs[2]/test.survdiff$exp[2])
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
  title(main="GSE17538(low-risk with ACT)",cex.main=1, adj=0.4)
  legend("topright", legend=c(paste0("no response (N=",table(group)["high"],")"),paste0("response (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
}

#######################################################################3
##TCGA and Nanostring
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
  
  ##1-2-3均值表达
  nanostring_23_expAverage <- nanostring.exp.process(nanostring.exp = nanostring_23_exp,
                                                     method = "average")
  ##1-2-3中值表达
  nanostring_23_expMedian <- nanostring.exp.process(nanostring.exp = nanostring_23_exp,
                                                    method = "median")
  ##1-2-3最大值表达
  nanostring_23_max <- nanostring.exp.process(nanostring.exp = nanostring_23_exp,
                                              method = "max")
  ##1-2-3最小值表达
  nanostring_23_min <- nanostring.exp.process(nanostring.exp = nanostring_23_exp,
                                              method = "min")
  ##1-2-3 MAD最大
  nanostring_23_mad <- nanostring.exp.process(nanostring.exp = nanostring_23_exp,
                                              method = "MAD")
}

#nanostring 用药样本
nanostring_23_drugs <- read.xlsx("0.New Nanostring/VUMCdataset_6.5.19_clinical_colon.xlsx",sheet = 3,colNames = T,rowNames = F,check.names=F)
nanostring_23_drugs2 <- read.xlsx("0.New Nanostring/VUMCdataset_6.5.19_clinical_colon.xlsx",sheet = 4,colNames = T,rowNames = F,check.names=F)

test_cli <- nanostring_23_cli[nanostring_23_drugs$`NanoString.V3.#`,]
{
  ##5-years
  test_cli[which(test_cli[,1]>60),2] = 0
  test_cli[which(test_cli[,1]>60),1] = 60
}

if(F){##############not run
  test_matrix <- nanostring_23_expMedian[,nanostring_23_drugs$`NanoString.V3.#`]
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
  
  {
    #Nanostring微卫星不稳定性和基因突变情况
    response_group = nanostring_23_drugs$`NanoString.V3.#`[which(group == "low")]
    no_response_group = nanostring_23_drugs$`NanoString.V3.#`[which(group == "high")]
    nanostring.genomics = t(nanostring_cli[c(response_group,no_response_group),9:14])
    nanostring.genomics[grep("p.",nanostring.genomics)] = "MT"
    nanostring.genomics[grep("Unknown",nanostring.genomics)] = -1
    nanostring.genomics[grep("wt",nanostring.genomics)] = 0
    nanostring.genomics[grep("MT",nanostring.genomics)] = 1
    nanostring.genomics = apply(nanostring.genomics, 2, as.numeric)
    rownames(nanostring.genomics) = colnames(nanostring_cli)[9:14]
    #
    nanostring.genomics.response = nanostring.genomics[,response_group]
    nanostring.genomics.response.ph = pheatmap(nanostring.genomics.response,cluster_rows = F,cellwidth = 8, cellheight = 8,
                                               color = c("#ecf2f9","#1f3c88","#e00543"))
    response_group1 = response_group[nanostring.genomics.response.ph$tree_col$order]
    #
    nanostring.genomics.noresponse = nanostring.genomics[,no_response_group]
    nanostring.genomics.noresponse.ph = pheatmap(nanostring.genomics.noresponse,cluster_rows = F,cellwidth = 8, cellheight = 8,
                                                 color = c("#ecf2f9","#1f3c88","#e00543"))
    no_response_group1 = no_response_group[nanostring.genomics.noresponse.ph$tree_col$order]
    
    nanostring.genomics1 = nanostring.genomics[,c(response_group1,no_response_group1)]
    
    library(pheatmap)
    annotation_col = data.frame(
      GroupType = factor(rep(c("Response Group", "No Response Group"), c(length(response_group1),length(no_response_group1))))
    )
    rownames(annotation_col) = c(response_group1,no_response_group1)
    
    pheatmap(nanostring.genomics1, cluster_rows = F, cluster_cols = F, cellwidth = 12, cellheight = 12,
             fontsize=8, fontsize_row=8,color = c("#ecf2f9","#1f3c88","#e00543"),
             annotation_col = annotation_col)
    
    response.unknown = length(which(nanostring.genomics.response == -1))
    response.wt = length(which(nanostring.genomics.response == 0))
    response.mt = length(which(nanostring.genomics.response == 1))
    nonresponse.unknown = length(which(nanostring.genomics.noresponse == -1))
    nonresponse.wt = length(which(nanostring.genomics.noresponse == 0))
    nonresponse.mt = length(which(nanostring.genomics.noresponse == 1))
    
    high.low.mut <- matrix(c(response.unknown,response.wt,response.mt,nonresponse.unknown,nonresponse.wt,nonresponse.mt),nrow = 3,
                           dimnames = list(c("unknown","wt","mt"),
                                           c("response","non response")))
    high.low.mut_fishertest <- fisher.test(high.low.mut)
    print(high.low.mut_fishertest)
    #p value:0.001113
  }
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
  HR_temp <- (test.survdiff$obs[1]/test.survdiff$exp[1])/(test.survdiff$obs[2]/test.survdiff$exp[2])
  HR_temp <- round(HR_temp,3)
  p.val_temp <- round(1 - pchisq(test.survdiff$chisq, 1),4)
  up95_temp <- exp(log(as.numeric(HR_temp)) + qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
  up95_temp <- round(up95_temp,3)
  low95_temp <- exp(log(as.numeric(HR_temp)) - qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
  low95_temp <- round(low95_temp,3)
  c_index_temp <- rcorr.cens(-surv.df$score, Surv(surv.df$surtime, surv.df$censor))[1]
  c_index_temp <- round(c_index_temp,4)
  plot(survfit(Surv(surtime,censor)~group,data=surv.df),xlab="Survival in months", ylab="Survival probability",
       col=c("red","darkblue"),bty="l",lwd=4)
  title(main=paste("Nanostring", "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
  legend("topright", legend=c(paste0("ACT-nonbenefit (N=",table(group)["high"],")"),paste0("ACT-nonbenefit (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
  #生存图
  fit = survfit(Surv(surtime,censor)~group,data=surv.df)
  surv.p <- ggsurvplot(fit = fit,
                       data = surv.df,
                       pval = TRUE, 
                       pval.coord = c(40,0.1),
                       legend = c(0.8,0.3),
                       legend.title = "Group",
                       legend.labs = c("ACT-nonbenefit", "ACT-benefit"),
                       palette = c("red","darkblue"),
                       title = paste("VUMC", "(pvalue=", p.val_temp, ")", sep="" ),
                       xlab = "Survival in months",
                       risk.table = TRUE,       # show risk table.
                       tables.height = 0.2,
                       risk.table.y.text.col = T, # colour risk table text annotations.
                       risk.table.y.text = FALSE, # show bars instead of names in text annotations
                       tables.theme = theme_cleantable(),
                       linetype = 1,
                       surv.median.line = "hv",
                       ggtheme = theme_survminer()
  )
  surv.p$plot = surv.p$plot+annotate("text", x=0, y=0.1,
                                     label=paste0(paste0("P value:",p.val_temp),"\n",
                                                  paste0(paste0("HR:",HR_temp)," (95%CI,",low95_temp,"-",up95_temp,")")),
                                     size = 5,hjust = 0)
  pdf(paste0("Script/revise4.26/figure7_2","VUMC",".pdf"), width = 6, height = 7)
  print(surv.p)
  dev.off()
}
library(pheatmap)
{
  #Nanostring微卫星不稳定性和基因突变情况
  response_group = nanostring_23_drugs$`NanoString.V3.#`[which(group == "low")]
  no_response_group = nanostring_23_drugs$`NanoString.V3.#`[which(group == "high")]
  nanostring.genomics = t(nanostring_cli[c(response_group,no_response_group),9:14])
  nanostring.genomics[grep("p.",nanostring.genomics)] = "MT"
  nanostring.genomics[grep("Unknown",nanostring.genomics)] = -1
  nanostring.genomics[grep("wt",nanostring.genomics)] = 0
  nanostring.genomics[grep("MT",nanostring.genomics)] = 1
  nanostring.genomics = apply(nanostring.genomics, 2, as.numeric)
  rownames(nanostring.genomics) = colnames(nanostring_cli)[9:14]
  #
  nanostring.genomics.response = nanostring.genomics[,response_group]
  nanostring.genomics.response.ph = pheatmap(nanostring.genomics.response,cluster_rows = F,cellwidth = 8, cellheight = 8,
                                             color = c("#ecf2f9","#1f3c88","#e00543"))
  response_group1 = response_group[nanostring.genomics.response.ph$tree_col$order]
  #
  nanostring.genomics.noresponse = nanostring.genomics[,no_response_group]
  nanostring.genomics.noresponse.ph = pheatmap(nanostring.genomics.noresponse,cluster_rows = F,cellwidth = 8, cellheight = 8,
                                               color = c("#ecf2f9","#1f3c88","#e00543"))
  no_response_group1 = no_response_group[nanostring.genomics.noresponse.ph$tree_col$order[c(1,9:15,2:8)]]
  
  nanostring.genomics1 = nanostring.genomics[,c(response_group1,no_response_group1)]
  
  library(pheatmap)
  annotation_col = data.frame(
    GroupType = factor(rep(c("Response Group", "No Response Group"), c(length(response_group1),length(no_response_group1))))
  )
  rownames(annotation_col) = c(response_group1,no_response_group1)
  
  pheatmap(nanostring.genomics1, cluster_rows = F, cluster_cols = F, cellwidth = 12, cellheight = 12,
           fontsize=8, fontsize_row=8,color = c("#ecf2f9","#1f3c88","#e00543"),
           annotation_col = annotation_col)
  
  response.unknown = length(which(nanostring.genomics.response == -1))
  response.wt = length(which(nanostring.genomics.response == 0))
  response.mt = length(which(nanostring.genomics.response == 1))
  nonresponse.unknown = length(which(nanostring.genomics.noresponse == -1))
  nonresponse.wt = length(which(nanostring.genomics.noresponse == 0))
  nonresponse.mt = length(which(nanostring.genomics.noresponse == 1))
  
  high.low.mut <- matrix(c(response.unknown,response.wt,response.mt,nonresponse.unknown,nonresponse.wt,nonresponse.mt),nrow = 3,
                         dimnames = list(c("unknown","wt","mt"),
                                         c("response","non response")))
  high.low.mut_fishertest <- fisher.test(high.low.mut)
  print(high.low.mut_fishertest)
  #p value:p-value = 1.148e-06
  
  #KRAS
  high.low.mut <- matrix(c(30,7,4,7,2,6),nrow = 3,
                         dimnames = list(c("unknown","wt","mt"),
                                         c("response","non response")))
  high.low.mut_fishertest <- fisher.test(high.low.mut)
  print(high.low.mut_fishertest)
  #p value:p-value = 0.04367
  
  #PIK3CA
  high.low.mut <- matrix(c(31,8,2,8,4,3),nrow = 3,
                         dimnames = list(c("unknown","wt","mt"),
                                         c("response","non response")))
  high.low.mut_fishertest <- fisher.test(high.low.mut)
  print(high.low.mut_fishertest)
  #p value:p-value = 0.04367
}

{#high
  group.nona = readRDS("3.18genes/Nanostring.RDS")
  group.nona.high = subset(group.nona, group=="high")
  group.nona.high.act = intersect(rownames(group.nona.high), rownames(test_cli))
  test_cli_high <- test_cli[group.nona.high.act,]
  #predict testing set
  test <- data.frame(t(test_matrix_z[hub_genes,group.nona.high.act]),stringsAsFactors = F)
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
  
  surv.df <- data.frame(test_cli_high[,1], test_cli_high[,2], score,stringsAsFactors = F)
  colnames(surv.df)[1:2] <- c("surtime", "censor")
  test.survdiff <- survdiff(Surv(surtime,censor)~group,data=surv.df)
  HR_temp <- (test.survdiff$obs[1]/test.survdiff$exp[1])/(test.survdiff$obs[2]/test.survdiff$exp[2])
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
  title(main="Nanostring(high-risk with ACT)",cex.main=1, adj=0.4)
  legend("topright", legend=c(paste0("ACT-nonbenefit (N=",table(group)["high"],")"),paste0("ACT-benefit (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
}
{#low
  group.nona = readRDS("3.18genes/Nanostring.RDS")
  group.nona.low = subset(group.nona, group=="low")
  group.nona.low.act = intersect(rownames(group.nona.low), rownames(test_cli))
  test_cli_low <- test_cli[group.nona.low.act,]
  #predict testing set
  test <- data.frame(t(test_matrix_z[hub_genes,group.nona.low.act]),stringsAsFactors = F)
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
  
  surv.df <- data.frame(test_cli_low[,1], test_cli_low[,2], score,stringsAsFactors = F)
  colnames(surv.df)[1:2] <- c("surtime", "censor")
  test.survdiff <- survdiff(Surv(surtime,censor)~group,data=surv.df)
  HR_temp <- (test.survdiff$obs[1]/test.survdiff$exp[1])/(test.survdiff$obs[2]/test.survdiff$exp[2])
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
  title(main="Nanostring(low-risk with ACT)",cex.main=1, adj=0.4)
  legend("topright", legend=c(paste0("ACT-nonbenefit (N=",table(group)["high"],")"),paste0("ACT-benefit (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
}

if(F){### not run
  ##################TCGA testing cohort
  ##TCGA
  #TCGA 用药样本
  {
    load("0.TCGA_XenaGDC/TCGA-COAD.Rdata")
    tcga_drug = readRDS("0.TCGA_XenaGDC/chemo_5fu_COAD.RDS")
    tcga_drug = intersect(paste0(tcga_drug,"-01A"),colnames(test_matrix))
    test_cli = test_cli[tcga_drug,2:3]
    test_matrix_z <- t(scale(t(test_matrix[,tcga_drug])))
    
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
    
    {
      ##5-years
      test_cli[which(test_cli[,1]>60),2] = 0
      test_cli[which(test_cli[,1]>60),1] = 60
    }
    
    surv.df <- data.frame(test_cli[,1], test_cli[,2], score,stringsAsFactors = F)
    colnames(surv.df)[1:2] <- c("surtime", "censor")
    test.survdiff <- survdiff(Surv(surtime,censor)~group,data=surv.df)
    HR_temp <- (test.survdiff$obs[1]/test.survdiff$exp[1])/(test.survdiff$obs[2]/test.survdiff$exp[2])
    HR_temp <- round(HR_temp,3)
    p.val_temp <- round(1 - pchisq(test.survdiff$chisq, 1),4)
    up95_temp <- exp(log(as.numeric(HR_temp)) + qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
    up95_temp <- round(up95_temp,3)
    low95_temp <- exp(log(as.numeric(HR_temp)) - qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
    low95_temp <- round(low95_temp,3)
    c_index_temp <- rcorr.cens(-surv.df$score, Surv(surv.df$surtime, surv.df$censor))[1]
    c_index_temp <- round(c_index_temp,4)
    plot(survfit(Surv(surtime,censor)~group,data=surv.df),xlab="Survival in months", ylab="Survival probability",
         col=c("red","darkblue"),bty="l",lwd=4)
    title(main=paste("TCGA", "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
    legend("topright", legend=c(paste0("high (N=",table(group)["high"],")"),paste0("low (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
    legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                                 paste0("p valu:",p.val_temp),
                                 paste0("HR:",HR_temp),
                                 paste0("up95:",up95_temp),
                                 paste0("low95:",low95_temp)))
    
    #生存图
    fit = survfit(Surv(surtime,censor)~group,data=surv.df)
    surv.p <- ggsurvplot(fit = fit,
                         data = surv.df,
                         pval = TRUE, 
                         pval.coord = c(40,0.1),
                         legend = c(0.8,0.3),
                         legend.title = "Group",
                         legend.labs = c("high risk", "low risk"),
                         palette = c("red","darkblue"),
                         title = paste("TCGA−COAD", "(pvalue=", p.val_temp, ")", sep="" ),
                         xlab = "Survival in months",
                         risk.table = TRUE,       # show risk table.
                         tables.height = 0.2,
                         risk.table.y.text.col = T, # colour risk table text annotations.
                         risk.table.y.text = FALSE, # show bars instead of names in text annotations
                         tables.theme = theme_cleantable(),
                         linetype = 1,
                         surv.median.line = "hv",
                         ggtheme = theme_survminer()
    )
    surv.p$plot = surv.p$plot+annotate("text", x=0, y=0.1,
                                       label=paste0(paste0("P value:",p.val_temp),"\n",
                                                    paste0(paste0("HR:",HR_temp)," (95%CI,",low95_temp,"-",up95_temp,")")),
                                       size = 5,hjust = 0)
    pdf(paste0("Script/revise4.26/figure6_","TCGA_COAD",".pdf"), width = 6, height = 7)
    print(surv.p)
    dev.off()
  }
  
  ##################TCGA CNV
  tcga_CNV <- read.table("0.TCGA_XenaGDC/colon_cancer_CNV.by_genes",header = T,sep = "\t",row.names = 1,
                         stringsAsFactors = F,check.names = F)
  colnames(tcga_CNV) = paste0(colnames(tcga_CNV),"A")
  mid_drug_CNV <- tcga_CNV[hub_genes,intersect(tcga_drug,colnames(tcga_CNV))]
  high_risk_sample <- tcga_drug[which(group == "high")]
  low_risk_sample <- tcga_drug[which(group == "low")]
  high_risk_CNV <- as.matrix(mid_drug_CNV[,high_risk_sample])
  high_risk_CNV[which(high_risk_CNV > 1)] <- 1
  low_risk_CNV <- as.matrix(mid_drug_CNV[,low_risk_sample])
  low_risk_CNV[which(low_risk_CNV > 1)] <- 1
  
  hub_gene_order <- c("DEPDC1","FRMD6","FBN1","AURKB","THY1","NAP1L3","BUB1","STON1","HSD17B2",
                      "ATAD2","TPX2")
  library(pheatmap)
  annotation_col1 = data.frame(
    GroupType = factor(rep("Bad Response", 27))
  )
  rownames(annotation_col1) = colnames(high_risk_CNV)
  high_risk_ph <- pheatmap(high_risk_CNV,cluster_cols = FALSE,cluster_rows = F,cellwidth =8, cellheight =8, 
                           fontsize=8, fontsize_row=8,color = c("#1f3c88","#4f98ca","#ecf2f9","#e00543"),
                           annotation_col = annotation_col1)
  
  annotation_col2 = data.frame(
    GroupType = factor(rep("Bad Response", 54))
  )
  rownames(annotation_col2) = colnames(low_risk_CNV)
  low_risk_ph <- pheatmap(low_risk_CNV,cluster_cols = FALSE,cluster_rows = F,cellwidth =8, cellheight =8, 
                          fontsize=8, fontsize_row=8,color = c("#1f3c88","#4f98ca","#ecf2f9","#e00543"),
                          annotation_col = annotation_col2)
  
  high.low.cnv <- matrix(c(94,176,247,293),nrow = 2,
                         dimnames = list(c("alternation","non alternation"),
                                         c("high risk","low risk")))
  high.low.cnv_fishertest <- fisher.test(high.low.cnv)
  #p value:7.778663e-06
  
  library(maftools)
  high_risk_sample.Barcodes = substr(high_risk_sample,1,12)
  low_risk_sample.Barcodes = substr(low_risk_sample,1,12)
  #mutation
  maf.coad = read.maf("0.TCGA_XenaGDC/TCGA.COAD.varscan.8177ce4f-02d8-4d75-a0d6-1c5450ee08b0.DR-10.0.somatic.maf",isTCGA = T)
  ##maf.coadread@clinical.data添加group信息
  maf.coad@clinical.data$Tumor_Sample_Barcode1 = paste0(maf.coad@clinical.data$Tumor_Sample_Barcode,"-01A")
  maf.coad@clinical.data$group = "NO"
  maf.coad@clinical.data[which(maf.coad@clinical.data$Tumor_Sample_Barcode1 %in% high_risk_sample), "group"] = "bad response"
  maf.coad@clinical.data[which(maf.coad@clinical.data$Tumor_Sample_Barcode1 %in% low_risk_sample), "group"] = "good response"
  
  ##针对不同group绘制突变图谱
  maf.coad@data$Tumor_Sample_Barcode1 = paste0(maf.coad@data$Tumor_Sample_Barcode,"-01A")
  #设置更改变异类型的颜色
  vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
  names(vc_cols) = c(
    'Frame_Shift_Del',
    'Missense_Mutation',
    'Nonsense_Mutation',
    'Multi_Hit',
    'Frame_Shift_Ins',
    'In_Frame_Ins',
    'Splice_Site',
    'In_Frame_Del'
  )
  
  gene = c("KRAS","BRAF","PIK3CA","PTEN","NRAS","AKT1")
  #group1 bad response
  maf.coad.group1 = subsetMaf(maf.coad, tsb = high_risk_sample.Barcodes)
  oncoplot(maf = maf.coad.group1, colors = vc_cols, clinicalFeatures = 'group', 
           top = 20, 
           genes = gene,
           logColBar = T)
  
  #group2 good response
  maf.coad.group2 = subsetMaf(maf.coad, tsb = low_risk_sample.Barcodes)
  oncoplot(maf = maf.coad.group2, colors = vc_cols, clinicalFeatures = 'group', 
           top = 20, 
           genes = gene,
           logColBar = T)
  
}

##################TCGA testing cohort
##TCGA
#TCGA 用药样本
{
  load("0.TCGA_XenaGDC/TCGA-COAD.Rdata")
  load("0.TCGA_XenaGDC/TCGA-READ.Rdata")
  test_cli = rbind(test_cli,test_cli_read)
  test_matrix = cbind(test_matrix,test_matrix_read)
  
  tcga_drug = c(readRDS("0.TCGA_XenaGDC/chemo_5fu_COAD.RDS"),
                readRDS("0.TCGA_XenaGDC/chemo_5fu_READ.RDS"))
  
  tcga_drug = intersect(paste0(tcga_drug,"-01A"),colnames(test_matrix))
  test_cli = test_cli[tcga_drug,2:3]
  test_matrix_z <- t(scale(t(test_matrix[,tcga_drug])))
  
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
  
  {
    ##5-years
    test_cli[which(test_cli[,1]>60),2] = 0
    test_cli[which(test_cli[,1]>60),1] = 60
  }
  
  surv.df <- data.frame(test_cli[,1], test_cli[,2], score,stringsAsFactors = F)
  colnames(surv.df)[1:2] <- c("surtime", "censor")
  test.survdiff <- survdiff(Surv(surtime,censor)~group,data=surv.df)
  HR_temp <- (test.survdiff$obs[1]/test.survdiff$exp[1])/(test.survdiff$obs[2]/test.survdiff$exp[2])
  HR_temp <- round(HR_temp,3)
  p.val_temp <- round(1 - pchisq(test.survdiff$chisq, 1),4)
  up95_temp <- exp(log(as.numeric(HR_temp)) + qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
  up95_temp <- round(up95_temp,3)
  low95_temp <- exp(log(as.numeric(HR_temp)) - qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
  low95_temp <- round(low95_temp,3)
  c_index_temp <- rcorr.cens(-surv.df$score, Surv(surv.df$surtime, surv.df$censor))[1]
  c_index_temp <- round(c_index_temp,4)
  plot(survfit(Surv(surtime,censor)~group,data=surv.df),xlab="Survival in months", ylab="Survival probability",
       col=c("red","darkblue"),bty="l",lwd=4)
  title(main="TCGA-COADREAD",cex.main=1, adj=0.4)
  legend("topright", legend=c(paste0("ACT-nonbenefit (N=",table(group)["high"],")"),paste0("ACT-benefit (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
  #生存图
  fit = survfit(Surv(surtime,censor)~group,data=surv.df)
  surv.p <- ggsurvplot(fit = fit,
                       data = surv.df,
                       pval = TRUE, 
                       pval.coord = c(40,0.1),
                       legend = c(0.8,0.3),
                       legend.title = "Group",
                       legend.labs = c("ACT-nonbenefit", "ACT-benefit"),
                       palette = c("red","darkblue"),
                       title = paste("TCGA−COAD", "(pvalue=", p.val_temp, ")", sep="" ),
                       xlab = "Survival in months",
                       risk.table = TRUE,       # show risk table.
                       tables.height = 0.2,
                       risk.table.y.text.col = T, # colour risk table text annotations.
                       risk.table.y.text = FALSE, # show bars instead of names in text annotations
                       tables.theme = theme_cleantable(),
                       linetype = 1,
                       surv.median.line = "hv",
                       ggtheme = theme_survminer()
  )
  surv.p$plot = surv.p$plot+annotate("text", x=0, y=0.1,
                                     label=paste0(paste0("P value:",p.val_temp),"\n",
                                                  paste0(paste0("HR:",HR_temp)," (95%CI,",low95_temp,"-",up95_temp,")")),
                                     size = 5,hjust = 0)
  pdf(paste0("Script/revise4.26/figure6_","TCGA_COAD",".pdf"), width = 6, height = 7)
  print(surv.p)
  dev.off()
}


high_risk_sample <- tcga_drug[which(group == "high")]
low_risk_sample <- tcga_drug[which(group == "low")]

library(maftools)
high_risk_sample.Barcodes = substr(high_risk_sample,1,12)
low_risk_sample.Barcodes = substr(low_risk_sample,1,12)
#mutation
maf.coad = read.maf("0.TCGA_XenaGDC/TCGA.COAD.varscan.8177ce4f-02d8-4d75-a0d6-1c5450ee08b0.DR-10.0.somatic.maf",isTCGA = T)
maf.read = read.maf("0.TCGA_XenaGDC/TCGA.READ.varscan.b2689e8f-3b64-4214-8a87-dc7e7cf6fe5e.DR-10.0.somatic.maf",isTCGA = T)

maf.coadread = merge_mafs(list(maf.coad,maf.read))
##maf.coadread@clinical.data添加group信息
maf.coadread@clinical.data$Tumor_Sample_Barcode1 = paste0(maf.coadread@clinical.data$Tumor_Sample_Barcode,"-01A")
maf.coadread@clinical.data$group = "NO"
maf.coadread@clinical.data[which(maf.coadread@clinical.data$Tumor_Sample_Barcode1 %in% high_risk_sample), "group"] = "bad response"
maf.coadread@clinical.data[which(maf.coadread@clinical.data$Tumor_Sample_Barcode1 %in% low_risk_sample), "group"] = "good response"

##针对不同group绘制突变图谱
maf.coadread@data$Tumor_Sample_Barcode1 = paste0(maf.coadread@data$Tumor_Sample_Barcode,"-01A")
#设置更改变异类型的颜色
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

gene = c("KRAS","BRAF","PIK3CA","PTEN","NRAS","AKT1")
#group1 bad response
maf.coadread.group1 = subsetMaf(maf.coadread, tsb = high_risk_sample.Barcodes)
oncoplot(maf = maf.coadread.group1, colors = vc_cols, clinicalFeatures = 'group', 
         top = 20, 
         genes = gene,
         logColBar = T)

#group2 good response
maf.coadread.group2 = subsetMaf(maf.coadread, tsb = low_risk_sample.Barcodes)
oncoplot(maf = maf.coadread.group2, colors = vc_cols, clinicalFeatures = 'group', 
         top = 20, 
         genes = gene,
         logColBar = T)









