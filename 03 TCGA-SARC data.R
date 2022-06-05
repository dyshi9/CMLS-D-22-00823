# For the RNA-Seq data
#-----------------------------------------------------------------------------------
library(TCGAbiolinks)
library(stringr)
library(readr)
library(SummarizedExperiment)
projectid = "TCGA-SARC"
query <- GDCquery(project = projectid, 
                  legacy = F, 
                  experimental.strategy = "RNA-Seq", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification",
                  sample.type = c("Primary Tumor"),     
                  workflow.type = "HTSeq - Counts")
GDCdownload(query)
dataAssy = GDCprepare(query, summarizedExperiment = F)
dataAssy_1 = dataAssy[,2:ncol(dataAssy)]
dataAssy_1 = data.matrix(dataAssy_1) # as.matrix
dataAssy_1 = as.data.frame(dataAssy_1)
rownames(dataAssy_1) = as.vector(dataAssy$X1)
dataAssy = dataAssy_1
class(dataAssy)
View(dataAssy[1:5,1:5])
rm(dataAssy_1)
###################################
colnames(dataAssy)[1]
colnames(dataAssy) = str_match(colnames(dataAssy), "(TCGA-[^-]*-[^-]*-[^-]*)")[,2]
dataAssyout = cbind(rownames(dataAssy), dataAssy)
colnames(dataAssyout)[1] = "Symbol"
dataAssyout$Symbol=as.character(dataAssyout$Symbol)
my_function=function(x) {x=str_split(x,"\\.")[[1]][1]}
dataAssyout$Symbol=apply(data.frame(dataAssyout$Symbol),1,my_function)
dataAssyout2=dataAssyout[-c(60484:60488),]
save(dataAssyout2,file = "TCGA-SARC_n259_rawcount.rdata")
################################################################################
query <- GDCquery(project = projectid, 
                  legacy = FALSE, 
                  experimental.strategy = "RNA-Seq", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification",
                  sample.type = c("Primary Tumor"),
                  workflow.type = "HTSeq - FPKM")
GDCdownload(query)
dataAssy = GDCprepare(query, summarizedExperiment = F)
dataAssy_1 = dataAssy[,2:ncol(dataAssy)]
dataAssy_1 = data.matrix(dataAssy_1) 
dataAssy_1 = as.data.frame(dataAssy_1)
rownames(dataAssy_1) = as.vector(dataAssy$X1)
dataAssy = dataAssy_1
class(dataAssy)
View(dataAssy[1:5,1:5])
rm(dataAssy_1)
colnames(dataAssy)[1]
colnames(dataAssy) = str_match(colnames(dataAssy), "(TCGA-[^-]*-[^-]*-[^-]*)")[,2]
dataAssyout = cbind(rownames(dataAssy), dataAssy)
colnames(dataAssyout)[1] = "Symbol"
dataAssyout$Symbol=as.character(dataAssyout$Symbol)
my_function=function(x) {x=str_split(x,"\\.")[[1]][1]}
dataAssyout$Symbol=apply(data.frame(dataAssyout$Symbol),1,my_function)
dataAssyout2=dataAssyout[-c(60484:60488),] #check row names
save(dataAssyout2,file = "TCGA-SARC_n259_fpkm.rdata")
#-----------------------------------------------------------------------------------
# Ensmeble ID to gene symbol
load("D:/CRAN/datasets/Homo_sapiens.GRCh38.104.chr.gtf.Rdata")
gtf_df[is.na(gtf_df$gene_name),"gene_name"] = gtf_df[is.na(gtf_df$gene_name),"gene_id"]
require(dplyr)
require(tidyr)
load("TCGA-SARC_n259_fpkm.rdata") # example / same for TCGA-SARC_n259_rawcount.rdata
GENEexpr = dataAssyout2
colnames(GENEexpr)[1] = "gene_id"
RNA_exprSet <- gtf_df %>% 
  dplyr::filter(type=="gene") %>% #筛选gene
  dplyr::select(c(gene_name,gene_id,gene_biotype)) %>% 
  dplyr::inner_join(GENEexpr,by ="gene_id") %>%   #exprfile读入表达矩阵
  tidyr::unite(gene_id,gene_name)
dup<-data.frame(table(RNA_exprSet$gene_id)) 
library("limma") # calculate mean if there are duplicated gene names
rt = RNA_exprSet
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,3:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
TCGA-SARC_n259_fpkm.symbol = data
save(SARC_n259_fpkm.symbol,file = "SARC_n259_fpkm.symbol.rdata")

# FPKM to TPM
#row-gene，col-sample
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
SARC_n259_tpm.symbol <- apply(SARC_n259_fpkm.symbol,2,fpkmToTpm)
#save


#-----------------------------------------------------------------------------------
# normalized raw count and output VST data via DESEQ2
library(DESeq2)
#load pan.sarcoma count
rawdata=TCGA.SARC_n259_count.symbol[rowMeans(TCGA.SARC_n259_count.symbol)>0,]
rawdata = round(rawdata)
dim(rawdata)
# exclude genes of extremely low abundance
rawdata <- rawdata[rowSums(rawdata)>5,]  #
dim(rawdata)
#
group <- c(rep('Tumor',259))
condition = factor(group)
sample <- data.frame(row.names = colnames(rawdata), condition)
ddsCountTable <- DESeqDataSetFromMatrix(countData = rawdata,
                                        colData = sample,
                                        design= ~ 1) 
ddsCountTable$condition<- relevel(ddsCountTable$condition, 
                                  ref = "Tumor") 
ddsCountTable <- estimateSizeFactors(ddsCountTable)
dim(ddsCountTable)
dds <- DESeq(ddsCountTable)
normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)
save(normalized_counts, file="normalized.counts_TCGA.SARC_n259.rdata")
# vst
vsd <- vst(dds,blind=T)
vsd <- assay(vsd)
colnames(vsd) = gsub("-01A","",colnames(vsd))
colnames(vsd) = gsub("-01B","",colnames(vsd))
vsd = as.data.frame(vsd)
save(vsd, file="vst_TCGA.SARC_n259.rdata")


#-----------------------------------------------------------------------------------
#TCGA-SARC subtypes information 
#type The Cancer Genome Atlas (TCGA) Research Network
library(TCGAbiolinks)
SARC.subtype <- TCGAquery_subtype(tumor = "SARC")
write.table(SARC.subtype,"SARC.subtype.txt",quote = F,sep = "\t",row.names = F)
#-----------------------------------------------------------------------------------
#Latest clinical information was downloaded directly form the GDC TCGA portal website
