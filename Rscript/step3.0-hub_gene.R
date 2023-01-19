rm(list = ls())
setwd("F:/colon cancer/0DataAndScript/1.primary data/train")
##1.GSE39582 was used to identify expression signature of probesets (ESP) by using Cox hazard analysis
#clinical information
GSE39582.clinal <- read.table("GSE39582.frma.sample.tsv",header = T,sep = "\t",stringsAsFactors = F)
GSE39582.clinal[1:5,1:8]
stageII_III <- which(1<GSE39582.clinal$"Stage" & GSE39582.clinal$"Stage"<4)
clinical_infor=cbind(GSE39582.clinal$DFSurvivalTime,GSE39582.clinal$DFSurvivalEvent)
clinical_infor=clinical_infor[stageII_III,]
colnames(clinical_infor)=c("time","event")
rownames(clinical_infor)=GSE39582.clinal$Sample[stageII_III]

#obtain stage II and III express data
GSE39582.frma <- read.table("GSE39582.frma.tsv",header = T,sep = "\t",stringsAsFactors = F)
rownames(GSE39582.frma) <- GSE39582.frma[,1]
GSE39582.frma <- GSE39582.frma[,-1]
GSE39582.frma23 <- GSE39582.frma[,GSE39582.clinal$Sample[stageII_III]]
dim(GSE39582.frma23)
GSE39582.frma23[1:5,1:8]

#remove the internal batch effect of GSE39582
library(sva)
batch <- read.csv("batch.csv", row.names=1)
rownames(batch) <- batch[,1]
batch <- batch[colnames(GSE39582.frma23),]
modcombat = model.matrix(~1, data=batch)
combat_edata = ComBat(dat=as.matrix(GSE39582.frma23), batch=batch$batchID, mod=modcombat, 
                      par.prior=TRUE, prior.plots=FALSE)
GSE39582.frma23.reba <- data.frame(combat_edata)
GSE39582.frma23.reba[1:5,1:8]

#identify expression signature of probesets (ESP)
library(survival)
#calculate cox p value for all probes 
clinical_infor_df <- as.data.frame(clinical_infor)
all.probe.p <- apply(as.matrix(GSE39582.frma23.reba),1,
                     function(prob_exp){
                       surv_curve=coxph(Surv(clinical_infor_df$time,clinical_infor_df$event)~as.numeric(prob_exp),data=clinical_infor_df)
                       coeffs=coef(summary(surv_curve))
                       p=coeffs[,5]
                       return(p)
                     })
head(all.probe.p)
cat("The number of probe sets<=0.05:",n=length(which(as.numeric(all.probe.p)<=0.05)))
cat("The number of probe sets<=0.01:",n=length(which(as.numeric(all.probe.p)<=0.01)))

#Convert to a gene symbol list corresponded to probe id
library("hgu133plus2.db")
all.sym <- as.list(hgu133plus2SYMBOL[mappedkeys(hgu133plus2SYMBOL)])
loc=c()
for(i in 1:length(all.sym)) {
  temp_loc=match(names(all.sym)[i],rownames(GSE39582.frma23.reba))
  loc=rbind(loc,temp_loc)
}
probe2gene_loc=loc  #54675 probes to 42323 probes with gene symbols, for train and test sets
GSE39582.frma23.reba.p2g <- GSE39582.frma23.reba[probe2gene_loc,]
dim(GSE39582.frma23.reba.p2g)

p2g.probe.p <- all.probe.p[probe2gene_loc]
cat("The number of probe sets<=0.05:",n=length(which(as.numeric(p2g.probe.p)<=0.05)))
cat("The number of probe sets<=0.01:",n=length(which(as.numeric(p2g.probe.p)<=0.01)))

all.gene=unique(as.character(all.sym))
cat("The number of genes=",n=length(unique(as.character(all.sym))))

#the significant probesets
sig_probe_loc=c() #42323 probes to 20514 genes, genes with significant P value, for train and test sets
for (i in 1:length(all.gene)) {
  list_loc=c() #locations of each gene to  probes, if multiples, record in "loc"
  loc=c()
  list_loc=which(as.character(all.sym) %in% all.gene[i]) 
  loc=as.numeric(which(p2g.probe.p[list_loc]==min(p2g.probe.p[list_loc])))  
  sig_probe_loc=c(sig_probe_loc,list_loc[loc])
}
GSE39582_exp_matrix=GSE39582.frma23.reba.p2g[sig_probe_loc,]
rownames(GSE39582_exp_matrix)=all.gene
dim(GSE39582_exp_matrix)
save.image(file = "ESP.RData")

##2.load the merge expression matrix, using space to find hub-gene
load("ESP.RData")
load("train_combat.RData")
cat("The dimension of training set=",n=dim(train_combat))
train_matrix=train_combat[probe2gene_loc,]
train_matrix=train_matrix[sig_probe_loc,]
rownames(train_matrix)=all.gene
head(train_matrix[1:10,1:20])

#########################
#p_value_list: p score list (must be matrix type) calculated by cox survival anlaysis; 
#method=1 or another none 1;
#p=0.05 or 0.01;
#degree=7 or other, p value computed by using space package
#fdr_p, p value for selecting genes in method 2
#########################
hub.extr=function(p_value_list,method,p, fdr_p) {
  if (method=="1") {
    ##select gene with significant cox-p vaulue
    loc_index=sort(p_value_list,index.return =T)
    gene_rank=rownames(p_value_list)[loc_index$ix]
    num=which(loc_index$x<=p)
    net_probe=gene_rank[num]
  } else { 
    ###select genes by FDR
    p_value_list_adjust <- p.adjust(p_value_list,method = "fdr")
    index <- which(p_value_list_adjust <= 0.25)
    net_probe=rownames(p_value_list)[index]
  }
  ##extract hub-gene matrix
  all_gene=net_probe
  temp=c()
  for(i in 1:length(all_gene)) {
    temp_loc=match(all_gene[i],rownames(train_matrix))
    temp=rbind(temp,temp_loc)
  }
  netgene_matrix=train_matrix[temp,]
  ##use space package to screen hub genes
  dat=t(netgene_matrix)  ## matrix: row-genes, col-samples
  library(space)
  names=rownames(netgene_matrix)
  n=nrow(dat)
  p=ncol(dat)
  alpha=1
  coef_lam=1.61
  l1=1/sqrt(n)*qnorm(1-alpha/(2*p^2))
  iter=3
  result=space.joint(dat, lam1=l1*n*coef_lam, lam2=0, weight=2, iter=iter)   #!!!check the matirix for "space.joint" matching: row-samples, col-genes!!!
  diag(result$ParCor)=0
  ###0.05
  fit.adj=abs(result$ParCor)>0.05
  sum(fit.adj==1)/2
  index <- which(fit.adj == T, arr.ind=T)
  index <- index[which(index[,1] > index[,2]),]
  tb <- sort(table(names[as.vector(index)]),decreasing=T)
  hub <- data.frame(gene = names(tb), connection=as.numeric(tb))
  return (hub)
}
#########################

p_value_list=matrix(p2g.probe.p[sig_probe_loc])
rownames(p_value_list)=all.gene
hub_probe=hub.extr(p_value_list,method="1",p=0.05,fdr_p=0.25)
hub_probe[1:25,]
save(hub_probe,file = "hub_probe.RData")

