#ICI genes------------------------------------------------------------------------------
load("vst_TCGA.SARC_n259.rdata")
expr = as.data.frame(t(vsd))
ICI = read.table("ICI gene symbol.txt",header = T)
ICI$ID %in% colnames(expr) #check
expr = expr[,colnames(expr) %in% ICI$ID]
expr$ID = rownames(expr)
clin = read.table("SARC_n261_210914.txt",
                  header = T,sep = "\t",check.names = F)
clin = clin[!is.na(clin$OS.Status),]
colnames(clin)[1]="ID"
expr = merge(expr,clin,by.x = "ID",by.y = "ID")
expr = expr[,c(1:45)] #useless cols
risklevel = read.table("Risktable_tcga.sarc_0.4668533.txt",header = T,check.names = F,sep = "\t")
risklevel$ID = row.names(risklevel)
risklevel = risklevel[,10:12]
expr.risk = merge(expr,risklevel,by.x = "ID",by.y = "ID")
rm(clin,expr,vsd,risklevel,ICI)
#
mean.val = expr.risk[,c(2:45,47)]
mean.val.low = mean.val[mean.val$EILncSig_Level == "Low",]
mean.val.high = mean.val[mean.val$EILncSig_Level == "High",]

mean.val.low.mean = apply(mean.val.low[,1:44], 2, mean) # 2 - column
mean.val.high.mean = apply(mean.val.high[,1:44], 2, mean) # 2 - column

mean_val = cbind(ID = colnames(mean.val)[1:44],
                 Low.mean = mean.val.low.mean,
                 High.mena = mean.val.high.mean)
mean_val = as.data.frame(mean_val)
mean_val$diff = as.numeric(mean_val$High.mena) - as.numeric(mean_val$Low.mean) 
rm(mean.val,mean.val.high,mean.val.low)

# take VTCN1 as an example
library(ggplot2)
library(ggpubr)
compare_means(VTCN1 ~ EILncSig_Level,data = expr.risk)
colnames(expr.risk)
ggplot(expr.risk, 
       aes(x = EILncSig_Level,
           y = VTCN1)) +
  geom_violin(aes(fill=EILncSig_Level), trim = FALSE) + 
  geom_boxplot(width = 0.1)+
  scale_fill_manual(values = c("salmon","steelblue"))+
  stat_compare_means(
    label = "p.format",
    method = "wilcox.test",
    label.y = 15,
    hide.ns = T
    )+
  theme_light() + xlab(NULL) + 
  ylab("Normalized expression") + 
  ggtitle("VTCN1") +
  theme(axis.text.x = element_text(size = 0))  + coord_cartesian(ylim=c(2.0, 15.25))

#correlation
ggplot(expr.risk, aes(x=EILncSig_Score, y=VTCN1)) + 
  geom_point(color="black",size = 1.0) + 
  geom_smooth(method="lm", se=T)+
  # geom_smooth()+
  theme_light()+
  stat_cor(data=expr.risk, method = "spearman",
           digits = 3,
           label.y= 11.35)+
  labs(x="EILncSig Score", 
       y="Normalized expression", 
       title="VTCN1")
#3.0-3.0

#Immunotherapy cohort------------------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(ggstatsplot)
#download GSE176307 processed file from GEO
expr <- read.table("GSE176307_baci_rsem_RS_BACI_headers_tab.txt",header = T,sep = "\t",check.names = F)
library(DESeq2)
genenames = expr$ID
rawdata = as.matrix(expr[,-1])
rawdata = round(rawdata)
dim(rawdata)
dim(rawdata)
group <- c(rep('Tumor',90))
condition = factor(group)
sample <- data.frame(row.names = colnames(rawdata), condition)

ddsCountTable <- DESeqDataSetFromMatrix(countData = rawdata,
                                        colData = sample,
                                        design= ~ 1) #~ condition; ~ 1
ddsCountTable$condition<- relevel(ddsCountTable$condition, 
                                  ref = "Tumor") # 
ddsCountTable <- estimateSizeFactors(ddsCountTable)
dim(ddsCountTable)
dds <- DESeq(ddsCountTable)

normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)
vsd <- vst(dds,blind=T)
vsd <- assay(vsd)
vsd = as.data.frame(vsd)
normalized_counts = cbind(ID= genenames,
                          normalized_counts)
vsd = cbind(ID= genenames,
            vsd)
save(normalized_counts, file="DESeq2.normalized_counts_n90.rdata")
save(vsd, file="vst_n90.rdata")
#load EILncSig model - recheck
load("uniSigExp.rdata")
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
load("vst_n90.rdata")
expr = vsd
#check gene names and keep official symbols
"MIR22HG" %in% expr$ID
"MIR155HG" %in% expr$ID
"LINC01140" %in% expr$ID
"LBX2-AS1" %in% expr$ID
"MCM3AP-AS1" %in% expr$ID
"CTD-2284J15.1"%in% expr$ID
expr[expr$ID == "CTD-2284J15.1","ID"] = "ENSG00000254231" # information on the ensembl Human (GRCh37.p13)
"RP11-568N6.1"%in% expr$ID
expr[expr$ID == "RP11-568N6.1","ID"] = "ENSG00000260101" # information on the ensembl Human (GRCh37.p13)
#extract
expr = expr[expr$ID %in% multi_sign_gene,]
rownames(expr) = expr$ID
expr = expr[,-1]
expr = as.data.frame(t(expr))
expr$EILncSig_Score <- predict(multiCox,type="risk",newdata=expr)
#load clinical info
clin = read.table("GSE176307_info.txt",header = T,sep = "\t",check.names = F)
clin2 = read.table("GSE176307_BACI_Omniseq_Sample_Name_Key_submitted_GEO_v2.csv",header = T,sep = ",",check.names = F)
colnames(clin2)[1]="ID"
clin=clin[-12,] #11 / 12
clin[clin$ID == "BACI165_1","ID"] = "BACI165" #BACI165_1 / BACI165_2
clin = merge(clin,clin2,by.x="ID",by.y="ID")
clin[clin$ID=="BACI165",8]="RS-03238964" #RS-03238964 / RS-03239001
colnames(clin)[8] = "sampleid"
expr$ID= rownames(expr)
expr.clin = merge(expr,clin,by.x="ID",by.y="sampleid")

library("survminer")
surv_cutpoint(
  expr.clin,
  time = "PFS.Time",
  event = "PFS.Status",
  "EILncSig_Score",
  minprop = 0.3,
  progressbar = TRUE
)
surv_cutpoint(
  expr.clin,
  time = "OS.Time",
  event = "OS.Status",
  "EILncSig_Score",
  minprop = 0.3,
  progressbar = TRUE
)
#KM
expr.clin$EILncSig_level = as.vector(ifelse(expr.clin$EILncSig_Score > 1.207528,"High","Low"))
write.table(expr.clin,"risktable_1.207528.txt",quote = F,row.names = F)
library(survival)
rt = expr.clin
rt$PFS.Time = rt$PFS.Time/30
diff=survdiff(Surv(PFS.Time, PFS.Status) ~EILncSig_level,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = F)
fit <- survfit(Surv(PFS.Time, PFS.Status) ~EILncSig_level, data = rt)
pdf(file="pfs.pdf",width = 5.5,height =5.75)
ggsurvplot(fit, 
           data=rt,
           size = 1.28,
           fontsize = 2.5,
           title = "Progression-free Survival \nKim et al.'s ICB cohort (GSE176307)",
           conf.int= F,
           pval=paste0("p= ",pValue),
           pval.coord = c(23.5,1),
           pval.size=5.0,
           surv.median.line = "hv",
           risk.table= T,
           palette = c("salmon","steelblue"),
           legend.labs=c("High", "Low"),
           legend.title="EILncSig Risk level",
           xlab="Follow up time (months)",
           #break.time.by = 365,
           ggtheme = theme_light(),
           risk.table.y.text.col = T,
           risk.table.height = 0.18,
           risk.table.y.text = F,
           #ncensor.plot = T,
           #ncensor.plot.height = 0.18,
           conf.int.style = "ribbon")
dev.off()

rt = expr.clin
rt$OS.Time = rt$OS.Time/30
diff=survdiff(Surv(OS.Time, OS.Status) ~EILncSig_level,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = F)
fit <- survfit(Surv(OS.Time, OS.Status) ~EILncSig_level, data = rt)
pdf(file="os.pdf",width = 5.5,height =5.75)
ggsurvplot(fit, 
           data=rt,
           size = 1.28,
           fontsize = 2.5,
           title = "Overall Survival \nKim et al.'s ICB cohort (GSE176307)",
           conf.int= F,
           pval=paste0("p= ",pValue),
           pval.coord = c(23.5,1),
           pval.size=5.0,
           surv.median.line = "hv",
           risk.table= T,
           palette = c("salmon","steelblue"),
           legend.labs=c("High", "Low"),
           legend.title="EILncSig Risk level",
           xlab="Follow up time (months)",
           #break.time.by = 365,
           ggtheme = theme_light(),
           risk.table.y.text.col = T,
           risk.table.height = 0.18,
           risk.table.y.text = F,
           #ncensor.plot = T,
           #ncensor.plot.height = 0.18,
           conf.int.style = "ribbon")
dev.off()

expr.clin$ICB_Response = NA
expr.clin[expr.clin$Response %in% c("CR","PR"),"ICB_Response"] = "CR/PR"
expr.clin[expr.clin$Response %in% c("SD","PD"),"ICB_Response"] = "SD/PD"
ggbarstats(expr.clin,
           ICB_Response,
           EILncSig_level, 
           bar.proptest = T, 
           label = "none",
           palette = 'Set2',#
           results.subtitle = T,
           type = "nonparametric",
           proportion.test = F
) + labs(x = '',y = '',title="Association of EILncSig Level and ICB Response\nKim et al.'s ICB cohort (GSE176307)") + 
  theme_light()+
  theme(axis.text.x = element_text(size = 11.5)) +
  # scale_fill_manual(values = c("steelblue","forestgreen","darkorange","salmon","brown3"))+
  scale_fill_manual(values = c("darkorange1","forestgreen")) +
  coord_flip() #
#10.2*3.5
