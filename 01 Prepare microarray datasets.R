# set working direction
library(limma)
library(affyPLM)
library(arrayQualityMetrics)

# Take GSE142162 as an example
# Download .CEL files from GEO
Data<-ReadAffy()
eset.rma<-rma(Data) #RMA
exprs<-exprs(eset.rma) #
colnames(exprs)[1] #
colnames(exprs)<-substr(colnames(exprs),1,10) #
exprs<-as.data.frame(exprs)
#prepare clinical information of GSE142162 in txt file
inf0<-read.table("GSE142162_info.txt",sep = "\t",header = T,row.names = 1)
#
sampleinter = intersect(colnames(exprs),row.names(inf0))
exprs = exprs[,sampleinter]
inf0 = inf0[sampleinter,]
exprs<-exprs[,match(rownames(inf0),colnames(exprs))]  

save(exprs,file = "GSE142162.n79.probeid.rdata")

#############################
#arrayQualityMetrics
Exp<-ExpressionSet(as.matrix(exprs),AnnotatedDataFrame(inf0))
preparedData = prepdata(expressionset =Exp,
                        intgroup = c(),
                        do.logtransform = F)  #
arrayQualityMetrics(Exp)

#delete outiliar samples if find
exprs=exprs[,-c(79)]
inf0= inf0[colnames(exprs),]
write.table(inf0,"GSE142162_info_n78.txt",sep = "\t",quote = F)

#
library(hgu133plus2.db)
ids=toTable(hgu133plus2ENSEMBL)
#
probeid<-intersect(rownames(exprs),ids$probe_id)  #
length(probeid)   #
exprs<-exprs[rownames(exprs)%in%probeid,]  #
ids<-ids[ids$probe_id %in%probeid,]
length(unique(ids$ensembl_id)) #

table(rownames(exprs)%in%ids$probe_id)  #
exprs<-exprs[rownames(exprs)%in%ids$probe_id,] #
dim(exprs)
ids<-ids[match(rownames(exprs),ids$probe_id),]  #
dim(ids)
head(ids)
head(exprs)
tmp<-by(exprs,ids$ensembl_id,function(x) rownames(x)[which.max(rowMeans(x))]) 
#
probes<-as.character(tmp)
exprs<-exprs[rownames(exprs) %in% probes,] #
dim(exprs)
ids<-ids[ids$probe_id %in%rownames(exprs),]

exprs$probe_id = row.names(exprs)
geneMatrix<-merge(ids,exprs,by="probe_id")
rownames(geneMatrix)<-geneMatrix$ensembl_id

geneMatrix<-geneMatrix[,-c(1,2)]

GSE142162_affy.RMA.log2_n78 = geneMatrix
save(GSE142162_affy.RMA.log2_n78,file = "GSE142162_affy.RMA.log2_n78.rdata")
