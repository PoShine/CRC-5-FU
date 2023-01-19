rm(list = ls())
setwd("F:/colon cancer/0DataAndScript/")
load("1.primary data/train/ESP.RData")
load("1.primary data/train/train_combat.RData")
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
net_probe=gene_rank[num] ##3091
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
rownames(fit.adj) <- names
colnames(fit.adj) <- names
library(igraph)
fit.adj.igraph <- graph_from_adjacency_matrix(fit.adj,mode = "undirected")
gene_interaction <- as_edgelist(fit.adj.igraph)
save.image("2.SPACE/space_result.RData")
write.table(gene_interaction,file = "2.SPACE/space_gene_interaction.txt",quote = F,sep = "\t",row.names=F,col.names=F)
sum(fit.adj==1)/2
index <- which(fit.adj == T, arr.ind=T)
index <- index[which(index[,1] > index[,2]),]
tb <- sort(table(names[as.vector(index)]),decreasing=T)
hub <- data.frame(gene = names(tb), connection=as.numeric(tb))
write.table(hub,file = "2.SPACE/space_gene_degree.txt",quote = F,sep = "\t",row.names=F,col.names=F)







