#########################
## using frma package for normalizing expression datasets
########################
rm(list = ls())
#GSE29621
array_name <- "GSE29621"
main_path = "F:/colon cancer/0DataAndScript/0.GEOdatasets/"
library(affy)
library(limma)
setwd(paste0(main_path,array_name))
cel_data <- ReadAffy()

library(frma)
library(hgu133plus2frmavecs)
data("hgu133plus2frmavecs")
exp_data <- frma(cel_data,summarize = "robust_weighted_average",input.vecs = hgu133plus2frmavecs)
exp_matrix <- exprs(exp_data)
colnames(exp_matrix)=unlist(strsplit(colnames(exp_matrix),".CEL.gz"))
exp_matrix <- round(exp_matrix,2)
write.table(exp_matrix,file = paste0(array_name,".frma.expression.txt"),sep = "\t",quote = F)








































