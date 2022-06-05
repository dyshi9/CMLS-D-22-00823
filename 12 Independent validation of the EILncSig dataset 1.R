#Recheck model - gene and coefs value
load("uniSigExp.rdata") # export from the model testing process (same file)
library(survival)
multiCox=coxph(Surv(MFS.Time, MFS.Status) ~ ., data = uniSigExp)  #MFS.Time, MFS.Status
#multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)
multi_var_coefs <- multiCoxSum$coefficients
multi_sign_gene <- rownames(multi_var_coefs)
multi_var_coefs <- multi_var_coefs[,'coef']
length(multi_sign_gene)
multi_sign_gene 
multi_sign_gene=gsub("`","",multi_sign_gene)

#load the test dataset
load("Dvst_TCGA.SARC_n259.rdata")
clin = read.table("/SARC_n261_210914.txt",
                  header = T,sep = "\t",check.names = F)
clin.os.expr = clin[!is.na(clin$OS.Status),]
colnames(clin.os.expr)[1]="ID"

test_expr = vsd
test_expr = test_expr[multi_sign_gene,]
test_expr = as.data.frame(t(test_expr))
test_expr$ID = rownames(test_expr)
test_expr = merge(test_expr,clin.os.expr,by.x = "ID",by.y = "ID")
rownames(test_expr) = test_expr$ID
test_expr = test_expr[,c("OS.Time","OS.Status",multi_sign_gene)]
test_expr$EILncSig_Score <- predict(multiCox,type="risk",newdata=test_expr)


#
library(survivalROC)
#1-3-5-median AUC
n.span = 0.05*259^(-.2)
roc = survivalROC(Stime=test_expr$OS.Time, status=test_expr$OS.Status, marker = test_expr$EILncSig_Score, 
                  predict.time = 365*5, method="NNE",span = n.span)
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)

#
#ROC
pdf(file="ROC.pdf",width=5.75,height=5.75)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc1=survivalROC(Stime=test_expr$OS.Time, status=test_expr$OS.Status, marker = test_expr$EILncSig_Score, 
                 predict.time = 365*1, method="NNE",span = n.span)
roc2=survivalROC(Stime=test_expr$OS.Time, status=test_expr$OS.Status, marker = test_expr$EILncSig_Score, 
                 predict.time = 365*3, method="NNE",span = n.span)
roc3=survivalROC(Stime=test_expr$OS.Time, status=test_expr$OS.Status, marker = test_expr$EILncSig_Score, 
                 predict.time = 365*5, method="NNE",span = n.span)
plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='forestgreen', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("Time-dependent ROC curve for OS prediction"), 
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1)
abline(0,1)
lines(roc2$FP,roc2$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="steelblue",lwd=2)
lines(roc3$FP,roc3$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="salmon",lwd=2)
legend("bottomright", 
       c("1-year AUC: 0.606","3-year AUC: 0.604","5-year AUC: 0.632"),
       lwd=2,
       bty = "n",
       col=c("forestgreen","steelblue","salmon"))
dev.off()

#surv_cutpoint
library(survminer)
surv_cutpoint(
  test_expr,
  time = "OS.Time",
  event = "OS.Status",
  "EILncSig_Score",
  minprop = 0.3, #min30%
  progressbar = TRUE
)

#cutoff 
test_expr = as.data.frame(test_expr)
test_expr$Risk_level = as.vector(ifelse(test_expr$EILncSig_Score > 0.4668533,"High","Low"))
colnames(test_expr)[10] = "EILncSig_Score"
colnames(test_expr)[11] = "EILncSig_Level"
write.table(test_expr,"Risktable_tcga.sarc_0.4668533.txt",sep="\t",quote=F)

#KM
library(survival)
library(survminer)
risktable=read.table("Risktable_tcga.sarc_0.4668533.txt",header=T,sep="\t",check.names=F)

risktable$OS.Time = risktable$OS.Time/365   

diff=survdiff(Surv(OS.Time, OS.Status) ~EILncSig_Level,data = risktable)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = T)

fit <- survfit(Surv(OS.Time, OS.Status) ~EILncSig_Level, data = risktable)

pdf(file="survival os tcga.sarc.pdf",width = 5.5,height =5.75)
ggsurvplot(fit, 
           data=risktable,
           size = 1.28,
           fontsize = 2.5,
           title = "Overall Survival \nTCGA-SARC",
           conf.int= F,
           pval=paste0("p= ",pValue),
           pval.coord = c(11.5,1),
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

risktable=read.table("Risktable_tcga.sarc_0.4668533.txt",sep="\t",header=T,row.names = 1,check.names=F)
risktable=risktable[order(risktable$EILncSig_Score),]
riskClass=risktable[,"EILncSig_Level"]
lowLength=length(riskClass[riskClass=="Low"])
highLength=length(riskClass[riskClass=="High"])
line=risktable[,"EILncSig_Score"]
line[line>10]=10

cutoff = 0.4668533

#RiskScore
#pdf 6.5*3.0
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
color=as.vector(risktable$OS.Status)
color[color==1]="salmon" #red
color[color==0]="steelblue" #lighblue
#pdf(file="SurvStat.pdf",width = 12,height = 5)
plot(risktable$OS.Time/365,
     pch=19,
     cex=0.7,
     xlab="Patients ordered by EILncSig Score (Red: Deceased, Blue: Alive)",
     ylab="OS time (years)",
     col=color)
abline(v=lowLength,lty=2,lwd=2,col="brown")
#dev.off()

#Heatmap
risktable1=risktable[,3:9]
risktable1=t(risktable1)
risktable2=as.matrix(risktable1)

annotation=data.frame(EILncSig_Level=risktable[,ncol(risktable)])
rownames(annotation)=rownames(risktable)

my_colour = list(EILncSig_Level = c(High = "salmon", Low = "steelblue"))

bk = unique(c(seq(-1.75,1.75, length= 50)))
library(pheatmap)
pdf(file="heatmap 7 gene.pdf",height = 3.5,width = 7.5) 
pheatmap(risktable2,
         annotation_col = annotation,
         breaks = bk,
         #color = colorRampPalette(c("blue", "white", "red"))(100),
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_cols = F,
         scale = "row",
         main = "TCGA-SARC",
         #gaps_col = 46,
         fontsize_row=11,
         fontsize_col=3,
         annotation_colors = my_colour,
         show_colnames= F) 
dev.off()


################################################################################
#C-indec

library(survcomp)
rt=read.table("Risktable_tcga.sarc_0.4668533.txt",header=T,sep="\t",check.names=F,row.names=1)
cindex <- concordance.index(x=rt$EILncSig_Score,
                            surv.time = rt$OS.Time, 
                            surv.event = rt$OS.Status,
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
risktable = read.table("Risktable_tcga.sarc_0.4668533.txt",header = T,sep = "\t",check.names = F)
risktable$id = rownames(risktable)

clin = read.table("SARC_n261_210914.txt",header = T,sep = "\t",check.names = F)
clin = clin[,c(1,2,3,9,11,13,14)]

risktable =merge(risktable,clin,by.x = "id",by.y="sample_id")
risktable = na.omit(risktable)

#
risktable$EILncSig_Level = gsub("Low","1.Low",risktable$EILncSig_Level)
risktable$EILncSig_Level = gsub("High","2.High",risktable$EILncSig_Level)

risktable$Residual_tumor = gsub("R0","1.R0",risktable$Residual_tumor)
risktable[!risktable$Residual_tumor == "1.R0","Residual_tumor"] = "2.R1/R2"

risktable$Gender = gsub("Female","1.Female",risktable$Gender)
risktable$Gender = gsub("Male","2.Male",risktable$Gender)

risktable$Local_Recurrence = gsub("No","1.No",risktable$Local_Recurrence)
risktable$Local_Recurrence = gsub("Yes","2.Yes",risktable$Local_Recurrence)

risktable$Tumor_depth = gsub("Superficial","01.Superficial",risktable$Tumor_depth)
risktable$Tumor_depth = gsub("Deep","2.Deep",risktable$Tumor_depth)

risktable$Metastasis_at_diagnosis = gsub("No","1.No",risktable$Metastasis_at_diagnosis)
risktable$Metastasis_at_diagnosis = gsub("Yes","2.Yes",risktable$Metastasis_at_diagnosis)

colnames(risktable)
#step1£¬unicox
f <- as.formula(Surv(OS.Time,OS.Status) ~ Local_Recurrence)
cox <- coxph(f,data=risktable)
coxsum = summary(cox)
outtable = cbind(HR = coxsum$conf.int[,"exp(coef)"],
                 LCI = coxsum$conf.int[,"lower .95"],
                 UCI = coxsum$conf.int[,"upper .95"],
                 P = coxsum$coefficients[,"Pr(>|z|)"])
outtable

#step2£¬multicox 
f1 <- as.formula(Surv(OS.Time,OS.Status) ~ Age_at_diagnosis + Gender
                 + Residual_tumor + Local_Recurrence + Tumor_depth
                 + Metastasis_at_diagnosis + EILncSig_Level)

coxmulti <- coxph(f1,data=risktable)
coxsum = summary(coxmulti)
outtable = cbind(HR = coxsum$conf.int[,"exp(coef)"],
                 LCI = coxsum$conf.int[,"lower .95"],
                 UCI = coxsum$conf.int[,"upper .95"],
                 P = coxsum$coefficients[,"Pr(>|z|)"])
write.table(outtable,"out.txt",sep = "\t")

#coxinput
library(survminer)
#9.75*5.65
ggforest(model = coxmulti, data = risktable, main = 'Multivariate Cox Regression Analysis\n TCGA-SARC Overall Survival',
         fontsize = 1)



################################################################################

#COX-Regression with clin.feature
#load clin and risktable
#load clin and risktable
risktable = read.table("Risktable_tcga.sarc_0.4668533.txt",header = T,sep = "\t",check.names = F)
risktable$id = rownames(risktable)
clin = read.table("SARC_n261_210914.txt",header = T,sep = "\t",check.names = F)
clin = clin[,c(1,12,13,14)]
risktable =merge(risktable,clin,by.x = "id",by.y="sample_id")
#risktable = na.omit(risktable)

risktable = read.table("Risktable_tcga.sarc_0.4668533.txt",header = T,sep = "\t",check.names = F)
risktable$id = rownames(risktable)
clin1 = read.table("SARC.subtype.txt",header = T,sep = "\t")
clin1 = clin1[,c(1,19,40,46,47)]
risktable =merge(risktable,clin1,by.x = "id",by.y="patient")
#risktable = na.omit(risktable)

#--------------------------------------------------------------------------------

colnames(risktable)[13] = "Tumor_Necrosis"
risktable = risktable[!is.na(risktable$Tumor_Necrosis),]
library(ggplot2)
library(ggpubr)
risktable[risktable$Tumor_Necrosis %in% c("Extensive Necrosis",
                                          "Moderate Necrosis"),"binary"] = "EM"
risktable[risktable$Tumor_Necrosis %in% c("Focal necrosis",
                                          "No necrosis"),"binary"] = "FN"
my_comparisons <- list(c("EM","FN"))
hline <- data.frame(binary=c("EM","FN"),
                    Tumor_Necrosis=c("Extensive Necrosis","Moderate Necrosis",
                                     "Focal necrosis","No necrosis"),
                    v=c(1.3, 1.3))

ggplot(risktable, 
       aes(x = binary,
           y = EILncSig_Score,
           fill = Tumor_Necrosis)) +
  geom_boxplot(aes(group=interaction(Tumor_Necrosis,binary)), 
               width=0.65, position=position_dodge(width = 1),
               outlier.shape=NA) +
  scale_fill_manual(values = c("steelblue","darkorange","forestgreen","salmon")) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     label.y = 1.25)+
  geom_errorbar(data = hline,size=0.7, width=0.45,aes(y=v, ymax=v, ymin=v))+
  theme_light() + xlab(NULL) +  ylab("EILncSig Score")+
  ggtitle("Association of Tumor Necrosis and EILncSig Score \nTCGA-SARC") +
  theme(axis.text.x = element_blank()) + coord_cartesian(ylim=c(0, 1.42))

#6*3.75
#--------------------------------------------------------------------------------

colnames(risktable)[15] = "Metastatic_Disease"
risktable = risktable[!is.na(risktable$Metastatic_Disease),]
library(ggplot2)
library(ggpubr)

ggplot(risktable, 
       aes(x = Metastatic_Disease,
           y = EILncSig_Score,
           fill = Metastatic_Disease)) +
  geom_boxplot(aes(group=Metastatic_Disease), 
               width=0.65, position=position_dodge(width = 1),
               outlier.shape=NA) +
  scale_fill_manual(values = c("forestgreen","brown3")) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format",
                     label.y = 1.29)+
  theme_light() + xlab(NULL) +  ylab("EILncSig Score")+
  ggtitle("Association of Metastasis and EILncSig Score \nTCGA-SARC") +
  theme(axis.text.x = element_blank()) + coord_cartesian(ylim=c(0, 1.32))
#5.5*3.75
#--------------------------------------------------------------------------------

colnames(risktable)[13] = "Relapse"
risktable = risktable[!is.na(risktable$Relapse),]
library(ggplot2)
library(ggpubr)

ggplot(risktable, 
       aes(x = Relapse,
           y = EILncSig_Score,
           fill = Relapse)) +
  geom_boxplot(aes(group=Relapse), 
               width=0.65, position=position_dodge(width = 1),
               outlier.shape=NA) +
  scale_fill_manual(values = c("forestgreen","brown3")) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format",
                     label.y = 1.45)+
  theme_light() + xlab(NULL) +  ylab("EILncSig Score")+
  ggtitle("Association of Relapse and EILncSig Score \nTCGA-SARC") +
  theme(axis.text.x = element_blank()) + coord_cartesian(ylim=c(0, 1.52))
#5.5*3.75


#--------------------------------------------------------------------------------

colnames(risktable)[16] = "iCluster"
risktable = risktable[!is.na(risktable$iCluster),]
risktable$iCluster = paste0("C",risktable$iCluster)

ggplot(risktable, 
       aes(x = iCluster,
           y = EILncSig_Score,
           fill = iCluster)) +
  geom_boxplot(aes(group=iCluster), 
               width=0.65, position=position_dodge(width = 1),
               outlier.shape=NA) +
  scale_fill_manual(values = c("brown3","salmon","darkorange","forestgreen","steelblue")) +
  stat_compare_means(#comparisons = my_comparisons,
                     #method = "wilcox.test",
                     label = "p.format",
                     label.y = 1.60)+
  theme_light() + xlab(NULL) +  ylab("EILncSig Score")+
  ggtitle("Association of Molecular Subtype and EILncSig Score \nTCGA-SARC") +
  theme(axis.text.x = element_blank()) + coord_cartesian(ylim=c(0, 1.65))


library(ggstatsplot)
library(ggpubr)
library(pals)

pdf("icluster2.pdf",width = 8.2,height = 3.85)
ggbarstats(risktable,
           EILncSig_Level, 
           iCluster,
           bar.proptest = T, 
           label = "none",
           palette = 'Set2',#
           results.subtitle = T,
           type = "nonparametric",
           proportion.test = F
) + labs(x = '',y = '',title="Molecular Subtype (Integrative clustering)") + theme_light()+
  theme(axis.text.x = element_text(size = 11.5)) +
  coord_flip() #
dev.off()


#--------------------------------------------------------------------------------
#RFS

risktable = read.table("Risktable_tcga.sarc_0.4668533.txt",header = T,sep = "\t",check.names = F)
risktable$id = rownames(risktable)
clin1 = read.table("SARC.subtype.txt",header = T,sep = "\t")
clin1 = clin1[,c(1,19,20)]
risktable =merge(risktable,clin1,by.x = "id",by.y="patient")
risktable = na.omit(risktable)
colnames(risktable)[13] = "RFS.Status"
colnames(risktable)[14] = "RFS.Time"

risktable[risktable$RFS.Status == "Relapse","RFS.Status"] = 1
risktable[risktable$RFS.Status == "Non-Relapse","RFS.Status"] = 0
risktable$RFS.Status = as.numeric(risktable$RFS.Status)

#KM
library(survival)
library(survminer)

risktable$RFS.Time = risktable$RFS.Time/365  

diff=survdiff(Surv(RFS.Time, RFS.Status) ~EILncSig_Level,data = risktable)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = T)

fit <- survfit(Surv(RFS.Time, RFS.Status) ~EILncSig_Level, data = risktable)

pdf(file="survival rfs tcga.sarc.pdf",width = 5.5,height =5.75)
ggsurvplot(fit, 
           data=risktable,
           size = 1.28,
           fontsize = 2.5,
           title = "Relapse-free Survival \nTCGA-SARC",
           conf.int= F,
           pval=paste0("p= ",pValue),
           pval.coord = c(8.5,1),
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



#
cindex <- concordance.index(x=risktable$EILncSig_Score,
                            surv.time = risktable$RFS.Time, 
                            surv.event = risktable$RFS.Status,
                            method = "noether")
cindex$c.index
cindex$se
cindex$lower
cindex$upper
cindex$p.value




#
risktable$RFS.Status = as.numeric(risktable$RFS.Status)
library(survivalROC)
#1-3-5-median AUC
n.span = 0.05*206^(-.2)
roc = survivalROC(Stime=risktable$RFS.Time, status=risktable$RFS.Status, marker = risktable$EILncSig_Score, 
                  predict.time = 365*1, method="NNE",span = n.span)
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)

#
#ROC
pdf(file="ROC rfs.pdf",width=5.75,height=5.75)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc1=survivalROC(Stime=risktable$RFS.Time, status=risktable$RFS.Status, marker = risktable$EILncSig_Score, 
                 predict.time = 365*1, method="NNE",span = n.span)
roc2=survivalROC(Stime=risktable$RFS.Time, status=risktable$RFS.Status, marker = risktable$EILncSig_Score, 
                 predict.time = 365*3, method="NNE",span = n.span)
roc3=survivalROC(Stime=risktable$RFS.Time, status=risktable$RFS.Status, marker = risktable$EILncSig_Score, 
                 predict.time = 365*5, method="NNE",span = n.span)
plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='forestgreen', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("Time-dependent ROC curve for RFS prediction"), 
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1)
abline(0,1)
lines(roc2$FP,roc2$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="steelblue",lwd=2)
lines(roc3$FP,roc3$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="salmon",lwd=2)
legend("bottomright", 
       c("1-year AUC: 0.612","3-year AUC: 0.609","5-year AUC: 0.594"),
       lwd=2,
       bty = "n",
       col=c("forestgreen","steelblue","salmon"))
dev.off()
