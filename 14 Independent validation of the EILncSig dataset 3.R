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

#load expr matrix
load("pan.sarcoma.array.symbol.n1085.sva.rdata")
#load clin
load("sample.batch.n1085.rdata")
#
expr = pan.sarcoma.array.symbol.n1085.sva
clin = sample.batch.n1085
clin.os = clin[!is.na(clin$OS.Status),]
clin.os = clin.os[clin.os$DataSet_ID == "E-TABM-1202",]
expr.os = as.data.frame(t(expr))
expr.os$ID = rownames(expr.os)
data.os = merge(clin.os[,c(1,9,10)],expr.os,by.x = "ID",by.y="ID")
rownames(data.os)=data.os$ID
#
test_expr = data.os[,c("OS.Time","OS.Status",multi_sign_gene)]
test_expr$EILncSig_Score <- predict(multiCox,type="risk",newdata=test_expr)

#
library(survivalROC)
#1-3-5-median AUC
n.span = 0.05*101^(-.2)
roc = survivalROC(Stime=test_expr$OS.Time, status=test_expr$OS.Status, marker = test_expr$EILncSig_Score, 
                  predict.time = 365*1, method="NNE",span = n.span)
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
       c("1-year AUC: 0.617","3-year AUC: 0.571","5-year AUC: 0.605"),
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
test_expr$Risk_level = as.vector(ifelse(test_expr$EILncSig_Score > 1.221731,"High","Low"))
colnames(test_expr)[10] = "EILncSig_Score"
colnames(test_expr)[11] = "EILncSig_Level"
write.table(test_expr,"Risktable_E.TABM.1202_1.221731.txt",sep="\t",quote=F)


#KM
library(survival)
library(survminer)
risktable=read.table("Risktable_E.TABM.1202_1.221731.txt",header=T,sep="\t",check.names=F)

risktable$OS.Time = risktable$OS.Time/365   

diff=survdiff(Surv(OS.Time, OS.Status) ~EILncSig_Level,data = risktable)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = F)

fit <- survfit(Surv(OS.Time, OS.Status) ~EILncSig_Level, data = risktable)

pdf(file="survival os E-TABM-1202.pdf",width = 5.5,height =5.75)
ggsurvplot(fit, 
           data=risktable,
           size = 1.28,
           fontsize = 2.5,
           title = "Overall Survival \nWilliamson et al.'s RMS cohort (E-TABM-1202)",
           conf.int= F,
           pval=paste0("p= ",pValue),
           pval.coord = c(8.75,1),
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

risktable=read.table("Risktable_E.TABM.1202_1.221731.txt",sep="\t",header=T,row.names = 1,check.names=F)
risktable=risktable[order(risktable$EILncSig_Score),]
riskClass=risktable[,"EILncSig_Level"]
lowLength=length(riskClass[riskClass=="Low"])
highLength=length(riskClass[riskClass=="High"])
line=risktable[,"EILncSig_Score"]
line[line>10]=10

cutoff = 1.221731

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
         main = "Williamson et al.'s RMS cohort (E-TABM-1202)",
         #gaps_col = 46,
         fontsize_row=11,
         fontsize_col=3,
         annotation_colors = my_colour,
         show_colnames= F) 
dev.off()


################################################################################
#C-indec

library(survcomp)
rt=read.table("Risktable_E.TABM.1202_1.221731.txt",header=T,sep="\t",check.names=F,row.names=1)
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
risktable = read.table("Risktable_E.TABM.1202_1.221731.txt",header = T,sep = "\t",check.names = F)
risktable$id = rownames(risktable)

clin = read.table("E-TABM-1202_info_n101.txt",header = T,sep = "\t",check.names = F)
clin = clin[,c(1:6)]

risktable =merge(risktable,clin,by.x = "id",by.y="id")
risktable = na.omit(risktable)

#
risktable$EILncSig_Level = gsub("Low","1.Low",risktable$EILncSig_Level)
risktable$EILncSig_Level = gsub("High","2.High",risktable$EILncSig_Level)


colnames(risktable)
#unicox
f <- as.formula(Surv(OS.Time,OS.Status) ~ EILncSig_Level)
cox <- coxph(f,data=risktable)
coxsum = summary(cox)
outtable = cbind(HR = coxsum$conf.int[,"exp(coef)"],
                 LCI = coxsum$conf.int[,"lower .95"],
                 UCI = coxsum$conf.int[,"upper .95"],
                 P = coxsum$coefficients[,"Pr(>|z|)"])
outtable



