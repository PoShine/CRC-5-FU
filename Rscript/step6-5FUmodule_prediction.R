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

## NCI-60 rna-seq
gene_zscore<-read.csv(file="5.NCI-60/rna-seq.csv",header = T, row.names = 1, stringsAsFactors = F)
gene_zscore <- scale(gene_zscore[,hub_genes])

## gi50
GI50=as.matrix(read.csv(file="5.NCI-60/5FU_gi50.csv", header = T, check.names=FALSE,na.strings ="-"))
GI50.mean <- apply(GI50, 2,
                   function(x){
                     mean(na.omit(x))
                   })

######################################################################################################
########################### 预测60个细胞系对5FU的治疗获益 ############################################
###########################
test <- data.frame(gene_zscore,stringsAsFactors = F)
pred <- predict(grow.RFS, newdata=test)
score <- pred$predicted
group <- rep(0, length(score))
for (i in 1:length(group)) {
  if( score[i]<=cut_off) {
    group[i]<-"ACT-benefit"
  } else{
    group[i]<-"ACT-nonbenefit"
  }
}

pred.df <- data.frame(cell.lines=rownames(test), pred.group=group, GI50.mean, stringsAsFactors = F)
pred.df$group <- ifelse(pred.df$pred.group=="ACT-benefit", 1, -1)

library(ggplot2)
ggplot(pred.df, aes(y=factor(cell.lines, levels = rev(cell.lines)), x=group, fill=GI50.mean)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks=c(-0.5,0.5), labels = c("ACT-nonbenefit","ACT-benefit"))+
  xlab("")+
  ylab("cell lines") +
  scale_fill_gradient(low = "#56B1F7",
                      high = "#132B43") +
  guides(fill=guide_colorbar(title = "−log10(GI50)",
                            title.position = "left",
                            title.theme = element_text(angle = 90))) +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank())
ggsave(filename = "5.NCI-60/2021.12.8revise/5FU_cellLines_modulePredict.pdf", width = 5, height = 8)


##############   5FU −log10(GI50) 比较 60细胞系
pred.df <- data.frame(cell.lines=rownames(test), pred.group=group, stringsAsFactors = F)

ACT.benefit <- rownames(test)[which(group=="ACT-benefit")]
ACT.nonbenefit <- rownames(test)[which(group=="ACT-nonbenefit")]

class1_matrix=GI50[,which(colnames(GI50) %in% ACT.benefit)]  ############组别一
class2_matrix=GI50[,which(colnames(GI50) %in% ACT.nonbenefit)]  ###########组别二

temp=c()
class1_list=c()
for (i in 1:dim(class1_matrix)[1]) {
  temp=mean(na.omit(class1_matrix[i,]))
  class1_list=c(class1_list,temp)
}

temp=c()
class2_list=c()
for (i in 1:dim(class2_matrix)[1]) {
  temp=mean(na.omit(class2_matrix[i,]))
  class2_list=c(class2_list,temp)
}

c<-wilcox.test(class1_list,class2_list, paired=TRUE,correct = TRUE) 
c$p.value

final_data<-data.frame(class1_list,class2_list)
final_data<-na.omit(final_data)

#############vioplot小提琴图
library(sm)
library (vioplot)
x1 <- final_data$class1_list
x2 <- final_data$class2_list
pdf("5.NCI-60/2021.12.8revise/5FU_60cc.pdf", width = 5, height = 8)
vioplot(x2,x1,names=c("ACT-nonbenefit","ACT-benefit"),col = c("#DB2D43","#1D97C1"), ylab = "−log10(GI50)")
title(main="5-Fluorouracil",cex.main=1, adj=0.4)
text(1,6.3,paste("wilcox.test, p-value < 2.2E-16"), cex = 0.8)
dev.off()


##############   5FU −log10(GI50) 只用比较 colon细胞系
test <- data.frame(gene_zscore[grep("CO:", rownames(gene_zscore)),hub_genes],stringsAsFactors = F)
pred <- predict(grow.RFS, newdata=test)
score <- pred$predicted
group <- rep(0, length(score))
for (i in 1:length(group)) {
  if( score[i]<=cut_off) {
    group[i]<-"ACT-benefit"
  } else{
    group[i]<-"ACT-nonbenefit"
  }
}

pred.df <- data.frame(cell.lines=rownames(test), pred.group=group, stringsAsFactors = F)

ACT.benefit <- rownames(test)[which(group=="ACT-benefit")]
ACT.nonbenefit <- rownames(test)[which(group=="ACT-nonbenefit")]

class1_matrix=GI50[,which(colnames(GI50) %in% ACT.benefit)]  ############组别一
class2_matrix=GI50[,which(colnames(GI50) %in% ACT.nonbenefit)]  ###########组别二

temp=c()
class1_list=c()
for (i in 1:dim(class1_matrix)[1]) {
  temp=mean(na.omit(class1_matrix[i,]))
  class1_list=c(class1_list,temp)
}

temp=c()
class2_list=c()
for (i in 1:dim(class2_matrix)[1]) {
  temp=mean(na.omit(class2_matrix[i,]))
  class2_list=c(class2_list,temp)
}

c<-wilcox.test(class1_list,class2_list, paired=TRUE,correct = TRUE) 
c$p.value

final_data<-data.frame(class1_list,class2_list)
final_data<-na.omit(final_data)

#############vioplot小提琴图
library(sm)
library (vioplot)
x1 <- final_data$class1_list
x2 <- final_data$class2_list
pdf("5.NCI-60/2021.12.8revise/5FU_7coloncc.pdf", width = 5, height = 8)
vioplot(x2,x1,names=c("ACT-nonbenefit","ACT-benefit"),col = c("#DB2D43","#1D97C1"), ylab = "−log10(GI50)")
title(main="5-Fluorouracil",cex.main=1, adj=0.4)
text(1,7,paste("wilcox.test, p-value < 2.2E-16"), cex = 0.8)
dev.off()








