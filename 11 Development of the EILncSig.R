#load WGCNA
wgcna = read.csv("geneinfo.lnc.csv",header = T,row.names = 1)
core.me = c("blue","purple","grey60","royalblue","darkorange")
wgcna.core = wgcna[wgcna$moduleColor %in% core.me,1:2]
emt.lnc=wgcna.core$geneSymbol
#load immune gene cor
immune.gene = read.table("cor_lnc&immune.gene_0.3.txt",header = T,sep = "\t")
#load immune pathway cor
immune.pathway = read.table("cor_lnc&immport_0.3.txt",header = T,sep = "\t")
#load tme cell cor
tme.cell = read.table("cor_lnc&tme.cell_origninal_p.txt",header = T,sep = "\t")
tme.cell = tme.cell[(abs(tme.cell$Cor) >=0.3)&(tme.cell$FDR<0.05),]
table(duplicated(tme.cell$Gene))
#intersect
immune.lnc = intersect(immune.gene$Gene,immune.pathway$Gene)
immune.lnc = intersect(immune.lnc,tme.cell$Gene)
EMTILncRNA = intersect(immune.lnc,emt.lnc)
table(duplicated(EMTILncRNA))
#------------------------------------------------------------------------------
#load expr matrix
load("pan.sarcoma.array.symbol.n1085.sva.rdata")
#load clinical info
load("sample.batch.n1085.rdata")
#
expr = pan.sarcoma.array.symbol.n1085.sva[EMTILncRNA,]
clin = sample.batch.n1085
#GSE71118 has the largest sample size with available clinical follow up data (mfs)
clin.mfs = clin[!is.na(clin$MFS.Status),]
expr.mfs = as.data.frame(t(expr))
expr.mfs$ID = rownames(expr.mfs)
data.mfs = merge(clin.mfs[,c(1,11,12)],expr.mfs,by.x = "ID",by.y="ID")
#unicox
library(survival)
allOutTab=data.frame()
UnicoxTab=data.frame()
VarNames = as.character(colnames(data.mfs[,4:ncol(data.mfs)]))
for(i in VarNames){
  rt=cbind(data.mfs[,2:3],data.mfs[,i])
  colnames(rt)[ncol(rt)]=i
  cox <- coxph(Surv(MFS.Time, MFS.Status) ~ rt[,ncol(rt)], data = rt)
  coxSummary = summary(cox)
  allOutTab=rbind(allOutTab,
                  cbind(id=i,
                        HR=coxSummary$conf.int[,"exp(coef)"],
                        HR.95L=coxSummary$conf.int[,"lower .95"],
                        HR.95H=coxSummary$conf.int[,"upper .95"],
                        pvalue=coxSummary$coefficients[,"Pr(>|z|)"]) 
  )}
#num
UnicoxTab = allOutTab[,-1]
UnicoxTab = as.data.frame(lapply(UnicoxTab, as.numeric))
rownames(UnicoxTab) = allOutTab[,1]
UnicoxTab[,1:3] = round(UnicoxTab[,1:3],3)
UnicoxTab[,4] = round(UnicoxTab[,4],5)
write.table(UnicoxTab,"unicox_mfs.txt",sep = "\t")
################################################################################
library(survival)
unicox.gene = read.table("unicox_mfs.txt",header = T,sep = "\t")
unicox.gene = row.names(unicox.gene[unicox.gene$pvalue <0.05,])
uniSigExp = cbind(data.mfs[,1:3],data.mfs[,unicox.gene])
rownames(uniSigExp) = uniSigExp[,1]
uniSigExp = uniSigExp[,-1]
save(uniSigExp,file = "uniSigExp.rdata")
################################################################################
load("uniSigExp.rdata")
multiCox=coxph(Surv(MFS.Time, MFS.Status) ~ ., data = uniSigExp)
multiCoxSum=summary(multiCox)
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
multicox.beta = outTab[,-1]
write.table(multicox.beta,file="multicox.coef.txt",sep="\t",row.names=T,quote=F)
################################################################################
load("uniSigExp.rdata")
library(survival)
multiCox=coxph(Surv(MFS.Time, MFS.Status) ~ ., data = uniSigExp)  #MFS.Time, MFS.Status
multiCoxSum=summary(multiCox)
multi_var_coefs <- multiCoxSum$coefficients
multi_sign_gene <- rownames(multi_var_coefs)
multi_var_coefs <- multi_var_coefs[,'coef']
length(multi_sign_gene)
multi_sign_gene 
multi_sign_gene=gsub("`","",multi_sign_gene)
#
Risk_score = predict(multiCox,type="risk",newdata=uniSigExp)
coxGene = rownames(multiCoxSum$coefficients)
coxGene = gsub("`","",coxGene)
outCol = c("MFS.Time","MFS.Status",coxGene)
write.table(cbind(id=rownames(cbind(uniSigExp[,outCol],Risk_score)),cbind(uniSigExp[,outCol],Risk_score)),
            file="Risktable_training.txt",
            sep="\t",
            quote=F,
            row.names=F)
rm(list=ls())
#
library(survivalROC)
risktable = read.table("Risktable_training.txt",header = T,check.names = F,sep = "\t")
#1-3-5-median AUC
n.span = 0.05*311^(-.2)
roc = survivalROC(Stime=risktable$MFS.Time, status=risktable$MFS.Status, marker = risktable$Risk_score, 
                  predict.time = 365*3, method="NNE",span = n.span)
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
#
#ROC
pdf(file="ROC.pdf",width=5.75,height=5.75)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc1=survivalROC(Stime=risktable$MFS.Time, status=risktable$MFS.Status, marker = risktable$Risk_score, 
                 predict.time = 365*1, method="NNE",span = n.span)
roc2=survivalROC(Stime=risktable$MFS.Time, status=risktable$MFS.Status, marker = risktable$Risk_score, 
                 predict.time = 365*3, method="NNE",span = n.span)
roc3=survivalROC(Stime=risktable$MFS.Time, status=risktable$MFS.Status, marker = risktable$Risk_score, 
                 predict.time = 365*5, method="NNE",span = n.span)
plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='forestgreen', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("Time-dependent ROC curve for MFS prediction"), 
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1)
abline(0,1)
lines(roc2$FP,roc2$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="steelblue",lwd=2)
lines(roc3$FP,roc3$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="salmon",lwd=2)
legend("bottomright", 
       c("1-year AUC: 0.714","3-year AUC: 0.684","5-year AUC: 0.680"),
       lwd=2,
       bty = "n",
       col=c("forestgreen","steelblue","salmon"))
dev.off()
#surv_cutpoint
library(survminer)
risktable = read.table("Risktable_training.txt",header = T,check.names = F,sep = "\t")
surv_cutpoint(
  risktable,
  time = "MFS.Time",
  event = "MFS.Status",
  "Risk_score",
  minprop = 0.3, #min30%
  progressbar = TRUE
)
#cutoff 
risktable = as.data.frame(risktable)
risktable$Risk_level = as.vector(ifelse(risktable$Risk_score > 1.245959,"High","Low"))
colnames(risktable)[11] = "EILncSig_Score"
colnames(risktable)[12] = "EILncSig_Level"
write.table(risktable,"Risktable_training_1.245959.txt",sep="\t",quote=F)
#KM
library(survival)
library(survminer)
risktable=read.table("Risktable_training_1.245959.txt",header=T,sep="\t",check.names=F)
risktable$MFS.Time = risktable$MFS.Time/365   
diff=survdiff(Surv(MFS.Time, MFS.Status) ~EILncSig_Level,data = risktable)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = T)
fit <- survfit(Surv(MFS.Time, MFS.Status) ~EILncSig_Level, data = risktable)

pdf(file="survival mfs training.pdf",width = 5.5,height =5.75)
ggsurvplot(fit, 
           data=risktable,
           size = 1.28,
           fontsize = 2.5,
           title = "Metastasis-free Survival \nChibon et al.'s Sarcoma cohort (GSE71118)",
           conf.int= F,
           pval=paste0("p= ",pValue),
           pval.coord = c(10.5,1),
           pval.size=5.0,
           surv.median.line = "hv",
           risk.table= T,
           #palette = c("forestgreen","steelblue","mediumorchid","salmon"),
           #palette = c("steelblue","salmon"),
           palette = c("salmon","steelblue"),
           legend.labs=c("High", "Low"),
           legend.title="EILncSig Risk level",
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
#Riskplot 
library(ggthemes)
library(pheatmap)

risktable=read.table("Risktable_training_1.245959.txt",sep="\t",header=T,row.names = 1,check.names=F)
risktable=risktable[order(risktable$EILncSig_Score),]
riskClass=risktable[,"EILncSig_Level"]
lowLength=length(riskClass[riskClass=="Low"])
highLength=length(riskClass[riskClass=="High"])
line=risktable[,"EILncSig_Score"]
line[line>10]=10

cutoff = 1.245959

#RiskScore
#pdf 6.5*4.5
plot(line,
     type="p",
     pch=20,
     cex=0.7,
     xlab="Patients ordered by EILncSig Score (Red: High Risk level, Blue: Low Risk level)",
     ylab="EILncSig Score",
     col=c(rep("steelblue",lowLength),  #lighblue
           rep("salmon",highLength)))  #red
abline(h=cutoff,v=lowLength,lty=2,lwd=2,col="brown")

#SurvStat
color=as.vector(risktable$MFS.Status)
color[color==1]="salmon" #red
color[color==0]="steelblue" #lighblue
#pdf(file="SurvStat.pdf",width = 12,height = 5)
plot(risktable$MFS.Time/365,
     pch=19,
     cex=0.7,
     xlab="Patients ordered by EILncSig Score (Red: Metastatic, Blue: Non-Metastatic)",
     ylab="MFS time (years)",
     col=color)
abline(v=lowLength,lty=2,lwd=2,col="brown")
#dev.off()
#Heatmap
# colnames(risktable)[8] = "WWP1-AS1"
# colnames(risktable)[10] = "AFTPH-DT"
risktable1=risktable[,4:10]
rownames(risktable1) = risktable$id
risktable1=t(risktable1)
risktable2=as.matrix(risktable1)

annotation=data.frame(EILncSig_Level=risktable[,ncol(risktable)])
rownames(annotation)=risktable$id

my_colour = list(EILncSig_Level = c(High = "salmon", Low = "steelblue"))

library(pheatmap)
pdf(file="heatmap 7 gene.pdf",height = 3.5,width = 7.5) 
bk = unique(c(seq(-1.65,1.65, length= 50)))
pheatmap(risktable2, 
         breaks = bk,
         annotation_col = annotation,
         #color = colorRampPalette(c("blue", "white", "red"))(100),
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_cols = F,
         scale = "row",
         main = "Chibon et al.'s Sarcoma cohort (GSE71118)",
         #gaps_col = 46,
         fontsize_row=11,
         fontsize_col=3,
         annotation_colors = my_colour,
         show_colnames= F) 
dev.off()

################################################################################
#C-indec

library(survcomp)
rt=read.table("Risktable_training_1.245959.txt",header=T,sep="\t",check.names=F,row.names=1)
cindex <- concordance.index(x=rt$EILncSig_Score,
                            surv.time = rt$MFS.Time, 
                            surv.event = rt$MFS.Status,
                            method = "noether")
cindex$c.index
cindex$se
cindex$lower
cindex$upper
cindex$p.value

################################################################################

#COX-Regression with clin.feature
library(survival)
#load clin and risktable
risktable = read.table("Risktable_training_1.245959.txt",header = T,sep = "\t")
clin = read.table("GSE71118_info_n311.txt",header = T,sep = "\t")
clin$id = rownames(clin)
clin = clin[,c(7,5)]
risktable =merge(risktable,clin,by = "id")

risktable$EILncSig_Level = gsub("Low","01_Low",risktable$EILncSig_Level)
risktable$EILncSig_Level = gsub("High","02_High",risktable$EILncSig_Level)
#
colnames(risktable)
#unicox
f <- as.formula(Surv(MFS.Time,MFS.Status) ~ EILncSig_Score)
cox <- coxph(f,data=risktable)
coxsum = summary(cox)
outtable = cbind(HR = coxsum$conf.int[,"exp(coef)"],
                LCI = coxsum$conf.int[,"lower .95"],
                UCI = coxsum$conf.int[,"upper .95"],
                P = coxsum$coefficients[,"Pr(>|z|)"])
outtable

################################################################################
library(survival)
#load clin and risktable
risktable = read.table("Risktable_training_1.245959.txt",header = T,sep = "\t")
clin = read.table("GSE71118_info_n311.txt",header = T,sep = "\t")
clin$id = rownames(clin)
clin = clin[,c(7,5)]
risktable =merge(risktable,clin,by = "id")
colnames(risktable)[13] = "CINSARC"
library(ggplot2)
library(ggpubr)
ggplot(risktable, 
       aes(x = CINSARC,
           y = EILncSig_Score)) +
  geom_boxplot(aes(fill=CINSARC),
               width=0.65, 
               position=position_dodge(width = 1),
               outlier.shape=NA) +
  scale_fill_manual(values = c("forestgreen","brown3"))+
  stat_compare_means(
    label = "p.format",
    #method = "wilcox.test",
    hide.ns = T,
    size = 4.5, 
    label.y =3.12)+
  theme_light() + xlab(NULL) + 
  ylab("EILncSig Score") + 
  ggtitle("Association of CINSARC subgroup and EILncSig Score \nChibon et al.'s Sarcoma cohort (GSE71118)") +
  theme(axis.text.x = element_text(size = 0)) + coord_cartesian(ylim=c(0, 3.2))
#pdf 5.0*3.75
