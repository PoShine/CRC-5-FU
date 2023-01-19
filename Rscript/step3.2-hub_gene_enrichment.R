rm(list = ls())
setwd("F:/colon cancer/0DataAndScript/")
############################THY1 neighbour genes
library(clusterProfiler)
library(ggplot2)
genes <- unlist(read.table("2.SPACE/THY1_neighbour_genes.txt",stringsAsFactors = F))
##GO
go_result <- enrichGO(gene = genes,"org.Hs.eg.db",
                      keyType = "SYMBOL",
                      ont = "BP")
View(summary(go_result))
go_bp <- summary(go_result)[,c(2,5:7)]
sig_go_bp <- go_bp[1:15,]
sig_go_bp$Description <- factor(sig_go_bp$Description,levels = sig_go_bp$Description[sort(1:15,decreasing = T)])
gg <- ggplot(data = sig_go_bp,aes(x=Description, y=-log10(qvalue)))+
  geom_bar(stat = "identity",fill="#729d39",width = 0.5)+
  coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1,colour = "black"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))+
  ylab("Enrichment Score[-log10(Pvalue)]")+
  xlab(NULL)+
  labs(title="GO")
gg

##KEGG
genes <- bitr(genes,"SYMBOL","ENTREZID","org.Hs.eg.db")
kegg_result <- enrichKEGG(gene = genes$ENTREZID)
View(summary(kegg_result))
kegg <- summary(kegg_result)[,c(2,5:7)]
kegg$Description <- factor(kegg$Description,levels = kegg$Description[sort(1:4,decreasing = T)])
gk <- ggplot(data = kegg,aes(x=Description, y=-log10(qvalue)))+
  geom_bar(stat = "identity",fill="#729d39",width = 0.5)+
  coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1,colour = "black"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))+
  ylab("Enrichment Score[-log10(Pvalue)]")+
  xlab(NULL)+
  labs(title="GO")
gk

library(gridExtra)
result <- list(gg,gk)
grid.arrange(arrangeGrob(grobs = result,nrow = 2))



#######################18 hub genes
rm(list = ls())
setwd("F:/colon cancer/0DataAndScript/")
load("1.primary data/train/hub_probe.RData")
genes <- as.character(hub_probe$gene[1:18])
##GO
go_result <- enrichGO(gene = genes,"org.Hs.eg.db",
                      keyType = "SYMBOL",
                      ont = "BP")
View(summary(go_result))
go_bp <- summary(go_result)[,c(2,5:7)]
sig_go_bp <- go_bp[1:15,]
sig_go_bp$Description <- factor(sig_go_bp$Description,levels = sig_go_bp$Description[sort(1:15,decreasing = T)])
gg <- ggplot(data = sig_go_bp,aes(x=Description, y=-log10(pvalue)))+
  geom_bar(stat = "identity",fill="#729d39",width = 0.5)+
  coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1,colour = "black"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))+
  ylab("Enrichment Score[-log10(Pvalue)]")+
  xlab(NULL)+
  labs(title="GO")
gg

##KEGG
load("aaa.RData")
View(summary(aaa))

genes <- bitr(genes,"SYMBOL","ENTREZID","org.Hs.eg.db")
kegg_result <- enrichKEGG(gene = genes$ENTREZID)
View(summary(kegg_result))

kegg <- summary(aaa)[1:7,c(2,5:7)]
kegg$Description <- factor(kegg$Description,levels = kegg$Description[sort(1:7,decreasing = T)])
gk <- ggplot(data = kegg,aes(x=Description, y=-log10(pvalue)))+
  geom_bar(stat = "identity",fill="#729d39",width = 0.5)+
  coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1,colour = "black"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))+
  ylab("Enrichment Score[-log10(Pvalue)]")+
  xlab(NULL)+
  labs(title="KEGG")
gk

library(gridExtra)
result <- list(gg,gk)
grid.arrange(arrangeGrob(grobs = result,nrow = 2))


