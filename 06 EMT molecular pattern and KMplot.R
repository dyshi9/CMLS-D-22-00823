load("pan.sarcoma.array.symbol.n1085.sva.rdata") # expression data
load("sample.batch.n1085.rdata") # clinical information
load("EMTgenes.rdata") # EMT gene list
#-----------------------------------------------------------------
EMTsig1 = EMTsig
# 
EMTsig1 = EMTsig
#-----------------------------------------------------------------
library(ConsensusClusterPlus)
workDir="./cluster"
#
data = pan.sarcoma.array.symbol.n1085.sva
data = data[row.names(data) %in% EMTsig1$true.symbol,]
data = as.matrix(data)
dim(data)
# data: row-gene, col-sample
#
maxK=7
results = ConsensusClusterPlus(data,
                               maxK=maxK,
                               reps=1000,
                               pItem=0.8,
                               pFeature=1,
                               title=workDir,
                               clusterAlg="km",
                               distance="euclidean",
                               seed=520,
                               verbose=T,
                               plot= "png")
#
clusterNum=7 # save results form 2 to 7, and test eac.
cluster=results[[clusterNum]][["consensusClass"]]
outTab=cbind(colnames(data),cluster)
write.table(outTab,file="EMT_cluster.txt",sep="\t",quote=F)

#---------------------------------------------------------------------------------
#silhouette
library(cluster)
library(factoextra)
library(CancerSubtypes)

#CancerSubtypes
cluster.s = silhouette_SimilarityMatrix(results[[7]]$consensusClass,
                                        results[[7]]$consensusMatrix)
fviz_silhouette(cluster.s, 
                #c("lightblue", "dodgerblue4"), #color
                label = F, 
                print.summary = TRUE) 

#---------------------------------------------------------------------------------
# We kept EMT cluster number of 2.
#KM
library(survival)
rt=read.table("EMT_cluster.txt",header=T,sep="\t",check.names=F)
colnames(rt)[1] = "ID"
rownames(rt) = 1:nrow(rt)
#load 
sample.batch_mfs = sample.batch.n1085[!is.na(sample.batch.n1085$MFS.Status),]
sample.batch_mfs = sample.batch_mfs[sample.batch_mfs$DataSet_ID=="GSE71118",]
rt=merge(rt,sample.batch_mfs,by.x = "ID",by.y="ID")
#
rt.mfs = rt[,c("ID","cluster","DataSet_ID","MFS.Time","MFS.Status")]
rt.mfs$MFS.Time = rt.mfs$MFS.Time/365   
clusterNum=length(levels(factor(rt.mfs$cluster)))

diff=survdiff(Surv(MFS.Time, MFS.Status) ~cluster,data = rt.mfs)
pValue=1-pchisq(diff$chisq,df=clusterNum-1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = F)

fit <- survfit(Surv(MFS.Time, MFS.Status) ~ cluster, data = rt.mfs)


library(survminer)
pdf(file="survival MFS.pdf",width = 5.5,height =5.3)
ggsurvplot(fit, 
           data=rt.mfs,
           size = 1.28,
           fontsize = 2.5,
           title="Metastasis-free Survival",
           conf.int= F,
           pval=paste0("p= ",pValue),
           pval.coord = c(11,1),
           pval.size=5.0,
           surv.median.line = "hv",
           risk.table= T,
           palette = c("steelblue","salmon"),
           legend.labs=c("EMT_C1", "EMT_C2"),
           legend.title="EMT Signature",
           xlab="Follow up time (years)",
           #break.time.by = 365,
           ggtheme = theme_light(),
           risk.table.y.text.col = T,
           risk.table.height = 0.18,
           risk.table.y.text = F,
           #ncensor.plot = T,
           #ncensor.plot.height = 0.18,
           conf.int.style = "ribbon")
dev.off()


#--------------------------------------------------------------------------------
#KM
library(survival)
rt=read.table("EMT_cluster.txt",header=T,sep="\t",check.names=F)
colnames(rt)[1] = "ID"
rownames(rt) = 1:nrow(rt)
#load clin
sample.batch_os = sample.batch.n1085[!is.na(sample.batch.n1085$OS.Status),]
table(sample.batch_os$DataSet_ID)

sample.batch_os = sample.batch_os[sample.batch_os$DataSet_ID=="E-TABM-1202",]

rt=merge(rt,sample.batch_os,by.x = "ID",by.y="ID")

#
rt.os = rt[,c("ID","cluster","DataSet_ID","OS.Time","OS.Status")]
rt.os$OS.Time = rt.os$OS.Time/365   
clusterNum=length(levels(factor(rt.os$cluster)))

diff=survdiff(Surv(OS.Time, OS.Status) ~cluster,data = rt.os)
pValue=1-pchisq(diff$chisq,df=clusterNum-1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = F)

fit <- survfit(Surv(OS.Time, OS.Status) ~ cluster, data = rt.os)


library(survminer)
pdf(file="survival OS.pdf",width = 5.5,height =5.3)
ggsurvplot(fit, 
           data=rt.os,
           size = 1.28,
           fontsize = 2.5,
           title="Overall Survival",
           conf.int= F,
           pval=paste0("p= ",pValue),
           pval.coord = c(8.75,1),
           pval.size=5.0,
           surv.median.line = "hv",
           risk.table= T,
           #palette = c("forestgreen","steelblue","mediumorchid","salmon"),
           #palette = c("steelblue","salmon"),
           palette = c("steelblue","salmon"),
           legend.labs=c("EMT_C1", "EMT_C2"),
           legend.title="EMT Signature",
           xlab="Follow up time (years)",
           #break.time.by = 365,
           ggtheme = theme_light(),
           risk.table.y.text.col = T,
           risk.table.height = 0.18,
           risk.table.y.text = F,
           #ncensor.plot = T,
           #ncensor.plot.height = 0.18,
           conf.int.style = "ribbon")
dev.off()


#------------------------------------------------------------
#

rt=read.table("EMT_cluster.txt",header=T,sep="\t",check.names=F)
colnames(rt)[1] = "ID"
rownames(rt) = 1:nrow(rt)
#load clin
sample.batch_mfs.stat = sample.batch.n1085[!is.na(sample.batch.n1085$Metastasis),]
rt=merge(rt,sample.batch_mfs.stat,by.x = "ID",by.y="ID")
colnames(rt)

library(ggstatsplot)
library(ggpubr)
library(pals)

rt[rt$cluster == "1","cluster"] = "EMT_C1"
rt[rt$cluster == "2","cluster"] = "EMT_C2"

rt1 = rt[rt$DataSet_ID == "GSE71118",]


pdf("EMT_Subtype metastasis.pdf",width = 8.85,height = 3.25)
ggbarstats(rt1,
           Metastasis, 
           cluster,
           bar.proptest = T, 
           #package = "pals",
           palette = 'Set1',#
           results.subtitle = T,
           type = "nonparametric",
           proportion.test = F
) + labs(x = '',y = '',title="Metastatic Disease") + theme_light()+
  theme(axis.text.x = element_text(size = 11.5)) +
  coord_flip() #
dev.off()
#8.75*3.5


