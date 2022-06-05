#prepare all available clinical information
#combine and put it in a csv or txt file.
sample.batch = read.csv("clin.csv",header = T)
colnames(sample.batch)[1] = "ID"
write.csv(sample.batch,"clin.csv",row.names = F)

# put all prepared microarray expression datasets in working dir

filelist = c("E.MEXP.1922_affy.RMA.rdata",
             "E.MEXP.3628_affy.RMA.rdata",
             "E.MEXP.964_affy.RMA.rdata",
             "E.TABM.1202_affy.RMA.rdata",
             "GSE13433_affy.RMA.rdata",
             "GSE142162_affy.RMA.rdata",
             "GSE14827_affy.RMA.rdata",
             "GSE17618_affy.RMA.rdata",
             "GSE20196_affy.RMA.rdata",
             "GSE20559_affy.RMA.rdata",
             "GSE23980_affy.RMA.rdata",
             "GSE34620_affy.RMA.rdata",
             "GSE34800_affy.RMA.rdata",
             "GSE37371_affy.RMA.rdata",
             "GSE66533_affy.RMA.rdata",
             "GSE71118_affy.RMA.rdata",
             "GSE87437_affy.RMA.rdata")
for (i in 1:17)
{load(filelist[i])}

rt1=E.MEXP.1922_affy.RMA
rt2=E.MEXP.3628_affy.RMA
rt3=E.MEXP.964_affy.RMA
rt4=E.TABM.1202_affy.RMA
rt5=GSE13433_affy.RMA
rt6=GSE142162_affy.RMA
rt7=GSE14827_affy.RMA
rt8=GSE17618_affy.RMA
rt9=GSE20196_affy.RMA
rt10=GSE20559_affy.RMA
rt11=GSE23980_affy.RMA
rt12=GSE34620_affy.RMA
rt13=GSE34800_affy.RMA
rt14=GSE37371_affy.RMA
rt15=GSE66533_affy.RMA
rt16=GSE71118_affy.RMA
rt17=GSE87437_affy.RMA

#check gene ids
rownames.common1 = intersect(row.names(rt1),row.names(rt2))
rownames.common2 = intersect(row.names(rt1),row.names(rt18))
#......
rownames.common = intersect(rownames.common1,rownames.common2)

rt1=rt1[rownames.common,]
rt2=rt2[rownames.common,]
rt3=rt3[rownames.common,]
rt4=rt4[rownames.common,]
rt5=rt5[rownames.common,]
rt6=rt6[rownames.common,]
rt7=rt7[rownames.common,]
rt8=rt8[rownames.common,]
rt9=rt9[rownames.common,]
rt10=rt10[rownames.common,]
rt11=rt11[rownames.common,]
rt12=rt12[rownames.common,]
rt13=rt13[rownames.common,]
rt14=rt14[rownames.common,]
rt15=rt15[rownames.common,]
rt16=rt16[rownames.common,]
rt17=rt17[rownames.common,]
rt18=rt18[rownames.common,]

pan.sarcoma.array.n1085 = cbind(rt1,rt2,rt3,rt4,rt5,rt6,rt7,rt8,rt9,
                                rt10,rt11,rt12,rt13,rt14,rt15,rt16,rt17,rt18)
dim(pan.sarcoma.array.n1085)

sample.common = intersect(colnames(pan.sarcoma.array.n1085),sample.batch$ID)

sample.batch = sample.batch[sample.batch$ID %in% sample.common,] # keep ID with expr
pan.sarcoma.array.n1085 = pan.sarcoma.array.n1085[,sample.batch$ID] #reorder cols for PCA
sample.batch.n1085 = sample.batch
save(pan.sarcoma.array.n1085,file = "pan.sarcoma.array.n1085.rdata")
save(sample.batch.n1085,file = "sample.batch.n1085.rdata")

#-----------------------------------------------------------------------------
load("pan.sarcoma.array.n1085.rdata")
load("sample.batch.n1085.rdata")
#-----------------------------------------------------------------------------

load("Homo_sapiens.GRCh38.104.chr.gtf.Rdata")
gtf_df[is.na(gtf_df$gene_name),"gene_name"] = gtf_df[is.na(gtf_df$gene_name),"gene_id"]
require(dplyr)
require(tidyr)

GENEexpr = pan.sarcoma.array.n1085
GENEexpr$gene_id = rownames(GENEexpr)

#id
exprSet <- gtf_df %>% 
  dplyr::filter(type=="gene") %>% #筛选gene
  dplyr::select(c(gene_name,gene_id,gene_biotype)) %>% 
  dplyr::inner_join(GENEexpr,by ="gene_id") %>%   #exprfile读入表达矩阵
  tidyr::unite(gene_id,gene_name)

#
dup<-data.frame(table(exprSet$gene_id)) 
library("limma")
rt =exprSet
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,3:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
pan.sarcoma.array.symbol.n1085 = data
pan.sarcoma.array.symbol.n1085 = as.data.frame(pan.sarcoma.array.symbol.n1085)
save(pan.sarcoma.array.symbol.n1085,file = "pan.sarcoma.array.symbol.n1085.rdata")

#-----------------------------------------------------------------------------

#check PCA
#
library(pals)
data = pan.sarcoma.array.symbol.n1085
#
data <- data[apply(data, 1, var)!=0,]
mads <- apply(data, 1, mad)
data <- data[rev(order(mads)),]
data <- data[1:500,]
dim(data)
data_t <- t(data)
variableL <- ncol(data_t)
data_t_m <- merge(data_t, sample.batch.n1085, by.x="row.names",by.y ="ID")
rownames(data_t_m) <- data_t_m$Row.names
data_t <- data_t_m[,-1]
library(factoextra)
pca <- prcomp(data_t[,1:variableL], scale=T)
col = as.vector(kelly(19))
col = col[2:19]
#
fviz_pca_ind(pca, col.ind=data_t$DataSet_ID, mean.point=F, addEllipses = F, label= "none", repel =T,
             geom = "point",
             palette = col, 
             legend.title="Project") +
  scale_shape_manual(values=rep(20,18)) + guides(color=guide_legend(override.aes=list(size = 3)))+
  xlab("Dim 1") + ylab("Dim 2") + theme_light()
#
fviz_pca_ind(pca, col.ind=data_t$DataSet_ID, mean.point=T, addEllipses = T, label= "none", repel =T,
             geom = "point",
             palette = col,
             ellipse.type="confidence", 
             ellipse.level=0.95,
             legend.title="Project") +  scale_shape_manual(values=rep(20,18))+ 
  guides(color=guide_legend(override.aes=list(size = 1)))+
  theme_light()


#------------------------------------------------------------------------------------
expr = pan.sarcoma.array.symbol.n1085
# row-gene col-sample
library(sva)
batchsample = sample.batch.n1085$DataSet_ID
expr_sva = ComBat(expr,batch = batchsample,mean.only=T)
# mean.only=T 
# the variances are expected to be different across batches due to the biology
dim(expr_sva)
#--------------------------------------------------------
rowneg = which(rowSums(expr_sva<=0)>0) #
colneg = which(colSums(expr_sva<=0)>0)
negmatrix = expr_sva[rowneg,colneg]
for(i in 1:length(rowneg))
{
  for (j in 1:length(colneg))
  {
    if(expr_sva[rowneg[i],colneg[j]]<0)
    {expr_sva[rowneg[i],colneg[j]]=0}
  }
}
table(expr_sva<=0)
table(expr_sva<0)
#--------------------------------------------------------
pan.sarcoma.array.symbol.n1085.sva = expr_sva
pan.sarcoma.array.symbol.n1085.sva = as.data.frame(pan.sarcoma.array.symbol.n1085.sva)
save(pan.sarcoma.array.symbol.n1085.sva,file = "pan.sarcoma.array.symbol.n1085.sva.rdata")

#check PCA
#
library(pals)
data = pan.sarcoma.array.symbol.n1085.sva

#
data <- data[apply(data, 1, var)!=0,]
mads <- apply(data, 1, mad)
data <- data[rev(order(mads)),]
data <- data[1:500,]
dim(data)

data_t <- t(data)
variableL <- ncol(data_t)
data_t_m <- merge(data_t, sample.batch.n1085, by.x="row.names",by.y ="ID")
rownames(data_t_m) <- data_t_m$Row.names
data_t <- data_t_m[,-1]
library(factoextra)
pca <- prcomp(data_t[,1:variableL], scale=T)

col = as.vector(kelly(19))
col = col[2:19]

fviz_pca_ind(pca, col.ind=data_t$DataSet_ID, mean.point=F, addEllipses = F, label= "none", repel =T,
             geom = "point",
             palette = col, 
             legend.title="Project") +
  scale_shape_manual(values=rep(20,18)) + guides(color=guide_legend(override.aes=list(size = 3)))+
  xlab("Dim 1") + ylab("Dim 2") + theme_light()

fviz_pca_ind(pca, col.ind=data_t$DataSet_ID, mean.point=T, addEllipses = T, label= "none", repel =T,
             geom = "point",
             palette = col,
             ellipse.type="confidence", 
             ellipse.level=0.95,
             legend.title="Project") +  scale_shape_manual(values=rep(20,18))+ 
  guides(color=guide_legend(override.aes=list(size = 1)))+ theme_light()
