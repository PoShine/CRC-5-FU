rm(list = ls())
setwd("F:/colon cancer/0DataAndScript")
load("1.primary data/train/ESP.RData")
load("1.primary data/train/train_combat.RData")
load("1.primary data/train/hub_probe.RData")
cat("The dimension of training set=",n=dim(train_combat))
train_matrix=train_combat[probe2gene_loc,]
train_matrix=train_matrix[sig_probe_loc,]
rownames(train_matrix)=all.gene
head(train_matrix[1:10,1:20])
#no fdr, p<=0.05
p_value_list=matrix(p2g.probe.p[sig_probe_loc])
rownames(p_value_list)=all.gene
loc_index=sort(p_value_list,index.return =T)
gene_rank=rownames(p_value_list)[loc_index$ix]
num=which(loc_index$x<=0.05)
net_probe=gene_rank[num]
#survival relatives genes(COX)
SR_gene <- p_value_list[net_probe,]
cat("survival relatives genes: ",length(SR_gene)) #3091
#####################
library(bioDist)
hub = as.character(hub_probe$gene[1:18])
sr = as.character(names(SR_gene)[1:18])
expr.hub = GSE39582_exp_matrix[hub,]
head(t(expr.hub))
expr.sr <- GSE39582_exp_matrix[sr,]
head(t(expr.sr))

## compare principle component ##
fit.hub <- prcomp(t(expr.hub))
fit.sr <- prcomp(t(expr.sr))
color <- c("#C82121","#1A3263")
gene_sets <- c("Hub gene set","Top-ranked gene set")
par(mfrow = c(1,2))
plot(fit.sr, ylim=c(0,7), type = "l", cex.lab=1.5, main="",col=color[2], lwd = 4)
plot(fit.hub, ylim=c(0,7), type = "l", cex.lab=1.5, main="",col=color[1], lwd = 4)
legend("topright", gene_sets,lty=1, cex = 0.7, col = color, bty = "n",lwd = 4)

## compare the mutualInfo distance ##
xdist.hub <- mutualInfo(as.matrix(expr.hub))
xdist.sr <- mutualInfo(as.matrix(expr.sr))
ks_test <- ks.test(xdist.sr, xdist.hub)
ecdf.hub <- ecdf(xdist.hub)
ecdf.sr <- ecdf(xdist.sr)
plot(ecdf.hub, col=color[1], lty=1, lwd=4, ylab="Cumulative distribution",
     xlab="Mutual Informatin Distance", cex.lab=1.5, pch=".",main="")
lines(ecdf.sr, col=color[2], lty=1, lwd=4, pch=".")
legend(0.18, 0.30, gene_sets, col=color, lwd=4, pch=".",
       bty="n", cex=1)
text(0.3,0.4,"KS test, p<0.00001")

## compare the median values of C-Indexes ##
library(randomForestSRC)
library(survival)
library(Hmisc)
library(YuGene)
sr_genes = as.character(names(SR_gene)[1:18])
GSE39582.RFS <- GSE39582_exp_matrix[sr_genes,]
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


##validated this prognosis model in the other seven independent cohorts with RFS time
test_data=paste0(c("GSE17536", "GSE17537", "GSE17538","GSE29621","GSE33113","GSE37892","GSE38832"),".RData")
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
  test <- data.frame(t(test_matrix[sr_genes,]),stringsAsFactors = F)
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
pre_result=rbind(c_index=c_index,p_val=p_val,HR=HR,low95=low95,up95=up95)
colnames(pre_result)=c("GSE17536", "GSE17537", "GSE17538","GSE29621","GSE33113","GSE37892","GSE38832")
print(t(pre_result))
top_genes.18.survival.result <- t(pre_result)[3:7,]
hub_genes.18.survival.result <- c(0.005,0.02,0.003,0,0.0314)
#
#load("hub_genes_18_survival_result.RData")
boxplot(-log10(hub_genes.18.survival.result), -log10(top_genes.18.survival.result[,2]),
        ylab="-log10(p value)",
        names=gene_sets, col=color,range=1)
median(-log10(hub_genes.18.survival.result))      ###2.30103
median(-log10(top_genes.18.survival.result[,2]))  ###0.5975667



##compare the gene expression patten##
expr.hub = GSE39582_exp_matrix[hub,]
dim(expr.hub)
expr.sr <- GSE39582_exp_matrix[sr,]
dim(expr.sr)
library(pheatmap)
hub_heatmap <- pheatmap(t(expr.hub),cluster_rows = F,cluster_cols = F,show_rownames=F)
top_heatmap <- pheatmap(t(expr.sr),cluster_rows = F,cluster_cols = F,show_rownames=F)
library(patchwork)
hub_heatmap + top_heatmap
expr.hub.sr.t.test <- t.test(expr.hub,expr.sr)


########################not run####
##compare 18 hub genes with randomly selected the same number of 18 genes##
###times: the time of randomly selected the same number of 18 genes from all 20514 genes
library(survminer)
prediction_performance <- function(numbers){
  random_result <- NULL
  test_data=paste0(c("GSE17536", "GSE17537", "GSE17538","GSE29621","GSE33113","GSE37892","GSE38832"),".RData")
  for (number in 1:numbers) {
    set.seed(number)
    random_genes <- rownames(GSE39582_exp_matrix)[sample(1:20514,18)]
    GSE39582.RFS <- GSE39582_exp_matrix[random_genes,]
    #prognosis model
    time <- clinical_infor_df$time
    event <- clinical_infor_df$event
    train.matrix <- data.frame(time = time, event = event, t(GSE39582.RFS), stringsAsFactors = F)
    set.seed(1244746)
    grow.RFS=rfsrc(Surv(time,event) ~.,data =train.matrix,ntree = 1000,nsplit = 2)
    grow.RFS
    
    #library(survminer)
    #######use the model.grow$predicted, which is in-bag predicted value, 
    #######to apply maxstat, then use the cutoff value on testing data;
    GSE39582_rfs_score = grow.RFS$predicted
    #summary(GSE39582_rfs_score)
    GSE39582_rfs_score.surv = cbind(na.omit(train.matrix[,1:2]),GSE39582_rfs_score)
    #head(GSE39582_rfs_score.surv)
    #Determining the optimal cutpoint for GSE39582_rfs_score
    GSE39582_rfs_score.cut = surv_cutpoint(
      data = GSE39582_rfs_score.surv,
      time = "time",
      event = "event",
      variables = "GSE39582_rfs_score",
      progressbar = F
    )
    #summary(GSE39582_rfs_score.cut)
    cut_off = GSE39582_rfs_score.cut$cutpoint$cutpoint
    #Plot the cutpoint for GSE39582_rfs_score
    #plot(GSE39582_rfs_score.cut, "GSE39582_rfs_score", palette = "npg")
    #Categorize GSE39582_rfs_score variable
    GSE39582_rfs_score.cat = surv_categorize(GSE39582_rfs_score.cut)
    ##validated this prognosis model in the other seven independent cohorts with RFS time
    c_index=c()
    for (k in 1:length(test_data)) {
      load(paste0("1.primary data/test/",test_data[k]))
      test_matrix <- test_matrix[probe2gene_loc,]
      test_matrix <- test_matrix[sig_probe_loc,]
      rownames(test_matrix) <- all.gene #assign gene symbols to matrix
      test <- data.frame(t(test_matrix[random_genes,]),stringsAsFactors = F)
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
      
      if(length(table(group)) < 2){
        c_index_temp = NA
      }else{
        surv.df <- data.frame(test_cli[,1], test_cli[,2], score,stringsAsFactors = F)
        colnames(surv.df)[1:2] <- c("surtime", "censor")
        test.survdiff <- survdiff(Surv(surtime,censor)~group,data=surv.df)
        c_index_temp <- rcorr.cens(-surv.df$score, Surv(surv.df$surtime, surv.df$censor))[1]
        c_index_temp <- round(c_index_temp,4)
        }
      c_index=c(c_index,c_index_temp)
    }
    random_result <- rbind(random_result,c_index)
    cat(paste0("has done numbers: ",number,"\n"))
  }
  colnames(random_result) <- c("GSE17536", "GSE17537", "GSE17538","GSE29621","GSE33113","GSE37892","GSE38832")
  rownames(random_result) = paste0("c_index",1:1000)
  return(random_result)
}
#random_18genes_result <- prediction_performance(1000)
random_result <- NULL
test_data=paste0(c("GSE17536", "GSE17537", "GSE17538","GSE29621","GSE33113","GSE37892","GSE38832"),".RData")
for (number in 1:1000) {
  set.seed(number)
  random_genes <- rownames(GSE39582_exp_matrix)[sample(1:20514,18)]
  GSE39582.RFS <- GSE39582_exp_matrix[random_genes,]
  #prognosis model
  time <- clinical_infor_df$time
  event <- clinical_infor_df$event
  train.matrix <- data.frame(time = time, event = event, t(GSE39582.RFS), stringsAsFactors = F)
  set.seed(1244746)
  grow.RFS=rfsrc(Surv(time,event) ~.,data =train.matrix,ntree = 1000,nsplit = 2)
  grow.RFS
  
  #library(survminer)
  #######use the model.grow$predicted, which is in-bag predicted value, 
  #######to apply maxstat, then use the cutoff value on testing data;
  GSE39582_rfs_score = grow.RFS$predicted
  #summary(GSE39582_rfs_score)
  GSE39582_rfs_score.surv = cbind(na.omit(train.matrix[,1:2]),GSE39582_rfs_score)
  #head(GSE39582_rfs_score.surv)
  #Determining the optimal cutpoint for GSE39582_rfs_score
  GSE39582_rfs_score.cut = surv_cutpoint(
    data = GSE39582_rfs_score.surv,
    time = "time",
    event = "event",
    variables = "GSE39582_rfs_score",
    progressbar = F
  )
  #summary(GSE39582_rfs_score.cut)
  cut_off = GSE39582_rfs_score.cut$cutpoint$cutpoint
  #Plot the cutpoint for GSE39582_rfs_score
  #plot(GSE39582_rfs_score.cut, "GSE39582_rfs_score", palette = "npg")
  #Categorize GSE39582_rfs_score variable
  GSE39582_rfs_score.cat = surv_categorize(GSE39582_rfs_score.cut)
  ##validated this prognosis model in the other seven independent cohorts with RFS time
  c_index=c()
  for (k in 1:length(test_data)) {
    load(paste0("1.primary data/test/",test_data[k]))
    test_matrix <- test_matrix[probe2gene_loc,]
    test_matrix <- test_matrix[sig_probe_loc,]
    rownames(test_matrix) <- all.gene #assign gene symbols to matrix
    test <- data.frame(t(test_matrix[random_genes,]),stringsAsFactors = F)
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
    if(length(table(group))<2 | min(table(group))<5){
      c_index_temp = NA
    }else{
      surv.df <- data.frame(test_cli[,1], test_cli[,2], score,stringsAsFactors = F)
      colnames(surv.df)[1:2] <- c("surtime", "censor")
      test.survdiff <- survdiff(Surv(surtime,censor)~group,data=surv.df)
      c_index_temp <- rcorr.cens(-surv.df$score, Surv(surv.df$surtime, surv.df$censor))[1]
      c_index_temp <- round(c_index_temp,4)
    }
    c_index=c(c_index,c_index_temp)
  }
  random_result <- rbind(random_result,c_index)
  cat(paste0("has done numbers: ",number,"\n"))
}
colnames(random_result) <- c("GSE17536", "GSE17537", "GSE17538","GSE29621","GSE33113","GSE37892","GSE38832")
rownames(random_result) = paste0("c_index",1:1000)
save(random_result,file = "random_18genes_result_5years.RData")

random_result_df = na.omit(as.data.frame(random_result))

load("random_18genes_result_5years.RData")
hub_genes.18.survival.result <- c(0.5999,0.6800,0.6199,0.7255,0.7228,0.8700,0.7292)
boxplot(hub_genes.18.survival.result, unlist(random_result_df), ylim=c(0.4,1),
        names=c("Hub gene set","Random 18 genes"), col=color)
boxplot(hub_genes.18.survival.result, colMeans(random_result_df), ylim=c(0.4,1),
        names=c("Hub gene set","Random 18 genes"), col=color,range=1.5)

par(mfrow = c(1,2))
