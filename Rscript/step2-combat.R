rm(list = ls())
setwd("F:/colon cancer/0DataAndScript/1.primary data/train")
#setwd("E:/colon cancer/primary data/train")
library(sva)
library(dplyr)
#merge the 6 expressiong matrixs
#exact the stage II-III patients
#GSE2109
GSE2109.frma <- read.table("GSE2109.frma.tsv",header = T,sep = "\t",stringsAsFactors = F)
rownames(GSE2109.frma) <- GSE2109.frma[,1]
GSE2109.frma <- GSE2109.frma[,-1]
GSE2109.clinal <- read.table("GSE2109.frma.sample.tsv",header = T,sep = "\t",stringsAsFactors = F)
GSE2109.sample <- GSE2109.clinal$Sample[c(grep("2",GSE2109.clinal$Stage),grep("3",GSE2109.clinal$Stage))]
GSE2109.frma <- GSE2109.frma[,GSE2109.sample]
#GSE18088  all patients are stage II
GSE18088.frma <- read.table("GSE18088.frma.tsv",header = T,sep = "\t",stringsAsFactors = F)
rownames(GSE18088.frma) <- GSE18088.frma[,1]
GSE18088.frma <- GSE18088.frma[,-1]
GSE18088.clinal <- read.table("GSE18088.frma.sample.tsv",header = T,sep = "\t",stringsAsFactors = F)
#GSE21510
GSE21510.frma <- read.table("GSE21510.frma.tsv",header = T,sep = "\t",stringsAsFactors = F)
rownames(GSE21510.frma) <- GSE21510.frma[,1]
GSE21510.frma <- GSE21510.frma[,-1]
GSE21510.clinal <- read.table("GSE21510.frma.sample.tsv",header = T,sep = "\t",stringsAsFactors = F)
GSE21510.sample <- GSE21510.clinal$Sample[c(grep("2",GSE21510.clinal$Stage),grep("3",GSE21510.clinal$Stage))]
GSE21510.frma <- GSE21510.frma[,GSE21510.sample]
#GSE26906  all patients are stage II
GSE26906.frma <- read.table("GSE26906.frma.tsv",header = T,sep = "\t",stringsAsFactors = F)
rownames(GSE26906.frma) <- GSE26906.frma[,1]
GSE26906.frma <- GSE26906.frma[,-1]
GSE26906.clinal <- read.table("GSE26906.frma.sample.tsv",header = T,sep = "\t",stringsAsFactors = F)
#GSE27854
GSE27854.frma <- read.table("GSE27854.frma.tsv",header = T,sep = "\t",stringsAsFactors = F)
rownames(GSE27854.frma) <- GSE27854.frma[,1]
GSE27854.frma <- GSE27854.frma[,-1]
GSE27854.clinal <- read.table("GSE27854.frma.sample.tsv",header = T,sep = "\t",stringsAsFactors = F)
GSE27854.sample <- GSE27854.clinal$Sample[c(grep("2",GSE27854.clinal$Stage),grep("3",GSE27854.clinal$Stage))]
GSE27854.frma <- GSE27854.frma[,GSE27854.sample]
#GSE39582
GSE39582.frma <- read.table("GSE39582.frma.tsv",header = T,sep = "\t",stringsAsFactors = F)
rownames(GSE39582.frma) <- GSE39582.frma[,1]
GSE39582.frma <- GSE39582.frma[,-1]
GSE39582.clinal <- read.table("GSE39582.frma.sample.tsv",header = T,sep = "\t",stringsAsFactors = F)
GSE39582.sample <- GSE39582.clinal$Sample[c(grep("2",GSE39582.clinal$Stage),grep("3",GSE39582.clinal$Stage))]
GSE39582.frma <- GSE39582.frma[,GSE39582.sample]

merge_expression <- cbind(GSE2109.frma,GSE18088.frma) %>% 
  cbind(GSE21510.frma) %>% 
  cbind(GSE26906.frma) %>%
  cbind(GSE27854.frma) %>%
  cbind(GSE39582.frma)
#save(merge_expression, file="merge_expression.RData")
#read the batch of samples
batch <- read.csv("batch.csv", row.names=1)
rownames(batch) <- batch[,1]
batch <- batch[colnames(merge_expression),]

#pca analysis
pca_result <- prcomp(t(merge_expression))
#draw image
colors = rainbow(14)
pointnames = c(97:110)
BN <- as.character(batch$batchID)
pdf(file="PCA_training_batch.pdf", width=9, height=9)
plot(pca_result$x[,1],pca_result$x[,2],xlab="PC1",ylab="PC2",main="PCA", col=colors[BN],type="n")
text(pca_result$x[,1],pca_result$x[,2],col=colors[as.numeric(BN)],BN, cex=0.8)
legend("bottomleft",14,legend=c(1:14),lty=1, col=colors)
dev.off()

#combat
modcombat = model.matrix(~1, data=batch)
combat_edata = ComBat(dat=as.matrix(merge_expression), batch=batch$batchID, mod=modcombat, 
                      par.prior=TRUE, prior.plots=TRUE)
train_combat <- data.frame(combat_edata)
save(train_combat, file="train_combat.RData")

#pca analysis after combat
pca_result_combat <- prcomp(t(combat_edata))
#draw image
colors = rainbow(14)
pointnames = c(97:110)
BN <- as.character(batch$batchID)
pdf(file="PCA_training_combat.pdf", width=9, height=9)
plot(pca_result_combat$x[,1],pca_result_combat$x[,2],xlab="PC1",ylab="PC2",main="PCA", col=colors[BN],type="n")
text(pca_result_combat$x[,1],pca_result_combat$x[,2],col=colors[as.numeric(BN)],BN, cex=0.8)
legend("bottomleft",14,legend=c(1:14),lty=1, col=colors)
dev.off()




