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
nanostring_23_drugs <- read.xlsx("0.New Nanostring/VUMCdataset_6.5.19_clinical_colon.xlsx",sheet = 3,colNames = T,rowNames = F,check.names=F)

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







