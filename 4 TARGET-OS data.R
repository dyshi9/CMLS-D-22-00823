# For the TARGET-OS RNA-Seq data and clinical information, We did not use TCGAbiolinks.
# We obtained them form the GDC TARGET website directly because of update (clinical information version 20200928)

# Download all RNA-Seqexpression txt file, and put them in working dir
#
file_names <- dir("\\CRAN\\datasets\\TARGET-OS\\RNA-seq\\gene", 
                  pattern = "*.txt", recursive = F, full.names = T)
#txt file name
file_id = as.data.frame(file_names)
#
library(stringr)
file_id[,1] = str_match(file_id[,1], "(gene/TARGET-[^-]*-[^-]*-[^-]*)")[,2]
file_id[,1] = str_match(file_id[,1], "(TARGET-[^-]*-[^-]*-[^-]*)")[,2]
#
merge.data <- (read.table(file_names[1],header = T,row.names = 1))
#column 3-TPM,column 4-counts
merge.data <- (read.table(file_names[1],header = T,row.names = 1))[3]
#从第2个txt文件开始合并
for (i in 2:length(file_names)) {
  new.data = read.table(file_names[i],header = T,row.names = 1)[3]
  merge.data <- cbind(merge.data, new.data)
  }
#TPM
expr_TPM <- merge.data
colnames(expr_TPM) <- as.vector(file_id[,1])
expr_TPM = as.data.frame(expr_TPM)
save(expr_TPM,file = "Target_OS.98_TPMRAW.Rdata")
#raw counts
expr_counts <- merge.data
colnames(expr_counts) <- as.vector(file_id[,1])
expr_counts = as.data.frame(expr_counts)
save(expr_counts,file = "Target_OS.98_countsRAW.Rdata")

#
require(dplyr)
require(tidyr)
load("Homo_sapiens.GRCh38.104.chr.gtf.Rda")
expr_TPM$gene_id = row.names(expr_TPM)
expr_TPM <- expr_TPM %>%
            select(gene_id, everything())
expr_TPM$gene_id = sub("\\..*","",expr_TPM$gene_id)
#转换全体id
expr.all_TPM <- gtf_df %>% 
  dplyr::filter(type=="gene") %>% 
  dplyr::select(c(gene_name,gene_id,gene_biotype)) %>% 
  dplyr::inner_join(expr_TPM,by ="gene_id") %>%   #exprfile读入表达矩阵
  tidyr::unite(gene_id,gene_name)
#check duplication
dup<-data.frame(table(expr.mRNA_TPM$gene_id)) 
#calculate mean expression if there are duplicated gene names
library("limma")
rt=as.matrix(expr.all_TPM)
rownames(rt)=rt[,1]
exp=rt[,3:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
save(data,file="expr.all_TPM.Rdata")
write.table(data,file="expr.all_TPM.txt",sep="\t",quote=F)
###############################################################################
expr_counts$gene_id = row.names(expr_counts)
expr_counts <- expr_counts %>%
  select(gene_id, everything())
expr_counts$gene_id = sub("\\..*","",expr_counts$gene_id)
expr.all_counts <- gtf_df %>% 
  dplyr::filter(type=="gene") %>% 
  dplyr::select(c(gene_name,gene_id,gene_biotype)) %>% 
  dplyr::inner_join(expr_counts,by ="gene_id") %>%   #exprfile读入表达矩阵
  tidyr::unite(gene_id,gene_name)
library("limma")
rt=as.matrix(expr.all_counts)
rownames(rt)=rt[,1]
exp=rt[,3:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=round(data)
save(data,file="expr.all_counts.Rdata")
write.table(data,file="expr.all_counts.txt",sep="\t",quote=F)

#----------------------------------------------------------------------------
library(DESeq2)
rawdata=TARGET.OS_n98_count.symbol[rowMeans(TARGET.OS_n98_count.symbol)>0,]
rawdata = round(rawdata)
dim(rawdata)

rawdata <- rawdata[rowSums(rawdata)>5,] 
head(rawdata)

group <- c(rep('Tumor',98))
condition = factor(group)

sample <- data.frame(row.names = colnames(rawdata), condition)

ddsCountTable <- DESeqDataSetFromMatrix(countData = rawdata,
                                        colData = sample,
                                        design= ~ 1) #~ condition; ~ 1
ddsCountTable$condition<- relevel(ddsCountTable$condition, 
                                  ref = "Tumor") 
ddsCountTable <- estimateSizeFactors(ddsCountTable)
dim(ddsCountTable)

dds <- DESeq(ddsCountTable)

normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)
save(normalized_counts, file="DESeq2.normalized_counts_n622.rdata")
vsd <- vst(dds,blind=T)
vsd <- assay(vsd)
colnames(vsd) = gsub("-01A","",colnames(vsd))
vsd = as.data.frame(vsd)
save(vsd, file="vst_TARGET.OS_n96.rdata")
