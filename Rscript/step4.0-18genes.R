rm(list = ls()); gc()
library(dplyr)
setwd("H:/210E盘文档/colon cancer/0DataAndScript")
load("1.primary data/train/ESP.RData")
load("1.primary data/train/hub_probe.RData")
##construct RSF prognosis model according to the RFS time of GSE39582
library(randomForestSRC)
library(survival)
library(Hmisc)
library(survminer)
hub_genes <- as.character(hub_probe$gene[1:18])
hub_genes
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

if(F){
  GSE39582_rfs_score.cat.high = subset(GSE39582_rfs_score.cat,GSE39582_rfs_score=="high")
  gse39582.highRisk = rownames(GSE39582_rfs_score.cat.high)
  save(gse39582.highRisk,file = "3.18genes/gse39582.highRisk.RData")
  
  fit <- survfit(Surv(time, event) ~ GSE39582_rfs_score,
                 data = GSE39582_rfs_score.cat)
  ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    surv.median.line = "hv", # Add medians survival
    risk.table = TRUE,       # show risk table.
    tables.height = 0.2,
    tables.theme = theme_cleantable(),
    pval = TRUE,             # show p-value of log-rank test.
    ggtheme = theme_bw(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = FALSE, # show bars instead of names in text annotations
    palette = c("red", "darkblue")
  )
}

dev.off()
##validated this prognosis model in the other seven independent cohorts with RFS time
test_data=paste0(c("GSE14333","GSE17536", "GSE17537", "GSE17538","GSE29621","GSE33113","GSE37892","GSE38832"),".RData")
c_index=c()
p_val=c()
HR=c()
up95=c()
low95=c()
for (k in 1:length(test_data)) {
  #load(paste0("/media/xp/文档/colon cancer/primary data/test/",test_data[k]))
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
  
  #test_cli$group = group
  #print(paste0(substr(test_data[k],1,8),".RDS"))
  #saveRDS(test_cli,file = paste0("3.18genes/",substr(test_data[k],1,8),".RDS"))
  
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
  c_index_temp <- round(c_index_temp,3)
  #plot(survfit(Surv(surtime,censor)~group,data=surv.df),xlab="Survival in months", ylab="Survival probability",
  #     col=c("red","darkblue"),bty="l",lwd=4)
  #title(main=paste(test_data[k], "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
  #length(grep("high",group))
  #legend("topright", legend=c(paste0("high (N=",table(group)["high"],")"),paste0("low (N=",table(group)["low"],")")), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  #legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
  #                             paste0("p valu:",p.val_temp),
  #                             paste0("HR:",HR_temp),
  #                             paste0("up95:",up95_temp),
  #                             paste0("low95:",low95_temp)))
  #生存图
  fit = survfit(Surv(surtime,censor)~group,data=surv.df)
  surv.p <- ggsurvplot(fit = fit,
                       data = surv.df,
                       pval = TRUE, 
                       pval.coord = c(40,0.25),
                       legend = c(0.8,0.5),
                       legend.title = "Group",
                       legend.labs = c("high risk", "low risk"),
                       palette = c("red","darkblue"),
                       title = paste(test_data[k], "(pvalue=", p.val_temp, ")", sep="" ),
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
                                     label=paste0(paste0("P value:",p.val_temp),"/n",
                                                  paste0(paste0("HR:",HR_temp)," (95%CI,",low95_temp,"-",up95_temp,")")),
                                     size = 5,hjust = 0)
  pdf(paste0("Script/revise4.26/figure3_",test_data[k],".pdf"), width = 6, height = 7)
  print(surv.p)
  dev.off()
  c_index=c(c_index,c_index_temp)
  p_val=c(p_val,p.val_temp)
  HR=c(HR,HR_temp)
  up95=c(up95,up95_temp)
  low95=c(low95,low95_temp)
}
pre_result=rbind(c_index=c_index,p_val=p_val,HR=HR,low95=low95,up95=up95)
colnames(pre_result)=c("GSE14333","GSE17536", "GSE17537", "GSE17538","GSE29621","GSE33113","GSE37892","GSE38832")
print(t(pre_result))


#######################################################################3
##TCGA and Nanostring
{
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
  #setwd("F:/colon cancer/0DataAndScript/0.New Nanostring")
  #nanostring表达数据
  load("0.New Nanostring/去除批次效应/nanostring.mRNA.norm.RData")
  #nanostring临床数据
  nanostring_cli <- read.xlsx("0.New Nanostring/VUMCdataset_6.5.19_clinical_colon.xlsx",sheet = 2,colNames = T,rowNames = F,check.names=F)
  rownames(nanostring_cli) <- nanostring_cli$`NanoString.V3.#`
  nanostring_cli <- nanostring_cli[,-1]
  #nanostring 2、3期样本
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

##nanostring_23_expAverage 0.0863;
##nanostring_23_expMedian 0.0393;
##nanostring_23_max 0.0124;
test_cli <- nanostring_23_cli
{
  ##5-years
  test_cli[which(test_cli[,1]>60),2] = 0
  test_cli[which(test_cli[,1]>60),1] = 60
}
test_matrix <- nanostring_23_max

if(T){
  #提取colon,right,left
  vumc_colon = unlist(read.table("0.New Nanostring/colon1.txt",header = F,stringsAsFactors = F)[,1])
  vumc_colon = intersect(vumc_colon, rownames(test_cli))
  test_cli = test_cli[vumc_colon,]
  test_matrix = test_matrix[,vumc_colon]
}

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
  
  test_cli$group = group
  #print(paste0(substr(test_data[k],1,8),".RDS"))
  #saveRDS(test_cli,file = "Nanostring.RDS")
  
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
                       pval.coord = c(40,0.25),
                       legend = c(0.8,0.5),
                       legend.title = "Group",
                       legend.labs = c("high risk", "low risk"),
                       palette = c("red","darkblue"),
                       title = paste(test_data[k], "(pvalue=", p.val_temp, ")", sep="" ),
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
                                     label=paste0(paste0("P value:",p.val_temp),"/n",
                                                  paste0(paste0("HR:",HR_temp)," (95%CI,",low95_temp,"-",up95_temp,")")),
                                     size = 5,hjust = 0)
  pdf(paste0("Script/revise4.26/figure7_1","VUMC",".pdf"), width = 6, height = 7)
  print(surv.p)
  dev.off()
}


#####TCGA
{
  load("F:/210E盘文档/colon cancer/TCGA/XenaGDC/TCGA-COAD.Rdata")
  #tcga_2 = rownames(test_cli)[which(test_cli$pathologic_stage == 'II')]
  #test_matrix = test_matrix[,tcga_2]
  test_cli = test_cli[,2:3]
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
  title(main=paste("TCGA", "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
  legend("topright", legend=c("high","low"), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
}


#####TCGA
{
  load("H:/210E盘文档/colon cancer/0DataAndScript/0.TCGA_XenaGDC/TCGA-COAD.Rdata")
  load("H:/210E盘文档/colon cancer/0DataAndScript/0.TCGA_XenaGDC/TCGA-READ.Rdata")
  test_cli = rbind(test_cli,test_cli_read)
  test_matrix = cbind(test_matrix,test_matrix_read)
  #tcga_2 = rownames(test_cli)[which(test_cli$pathologic_stage == 'II')]
  #test_matrix = test_matrix[,tcga_2]
  test_cli = test_cli[,2:3]
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
  title(main=paste("TCGA", "(pvalue=", p.val_temp, ")", sep="" ),cex.main=1, adj=0.4)
  legend("topright", legend=c("high","low"), lty=1,col=c("red", "darkblue"),bty = "n",lwd = 4)
  legend("bottomleft",legend=c(paste0("C index:",c_index_temp),
                               paste0("p valu:",p.val_temp),
                               paste0("HR:",HR_temp),
                               paste0("up95:",up95_temp),
                               paste0("low95:",low95_temp)))
  
}







