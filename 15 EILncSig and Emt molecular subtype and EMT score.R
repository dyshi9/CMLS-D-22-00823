#re-check model - gene and coefs value
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
load("pan.sarcoma.array.symbol.n1085.sva.rdata")
load("sample.batch.n1085.rdata")
expr = pan.sarcoma.array.symbol.n1085.sva
clin = sample.batch.n1085
emt.c =read.table("EMT_cluster.txt",
                  header = T,sep = "\t")
emt.c[emt.c$cluster == "1","cluster"] = "C1"
emt.c[emt.c$cluster == "2","cluster"] = "C2"
colnames(emt.c)[1] = "id"
colnames(emt.c)[2] = "EMT_Cluster"
test_expr = as.data.frame(t(expr))
test_expr = test_expr[,multi_sign_gene]
test_expr$riskscore <- predict(multiCox,type="risk",newdata=test_expr)
test_expr$id = rownames(test_expr)
test_expr = merge(test_expr,emt.c,by.x="id",by.y="id")

ggplot(test_expr, 
       aes(x = EMT_Cluster,
           y = riskscore,
           fill = EMT_Cluster)) +
  geom_violin(aes(fill=EMT_Cluster), trim = FALSE) + 
  geom_boxplot(width = 0.1,fill="white")+
  scale_fill_manual(values = c("steelblue","salmon")) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    size = 4.75, 
    label.y = 7.25)+
  theme_light() + xlab(NULL) +  ylab("EILncSig Score")+
  ggtitle("Association of EMT Clusters and EILncSig Score \nMerged Pan-Sarcoma dataset") +
  theme(axis.text.x = element_blank()) 

#5.5*3.75
#------------------------------------------------------------------------------
risktable = cbind(id=rownames(test_expr),score = test_expr$riskscore)
colnames(risktable)[2] = "EILncSig_Score"
risktable = as.data.frame(risktable)
risktable$EILncSig_Score = as.numeric(risktable$EILncSig_Score)
# from the step - EMT clusters analysis - EMT signature scores
# Hollern and Tuan - EMT score
emt.score = read.table("Tuan emt score.txt",sep = "\t",header = T) 

cortable = merge(risktable,emt.score,by.x="id",by.y="ID")


library(ggplot2)
library(ggpubr)
#
cortable = as.data.frame(cortable)

ggplot(cortable, aes(x=emt.score, y=EILncSig_Score)) + 
  geom_point(color="black",size = 1) + 
  geom_smooth(method="lm", se=T)+
  #geom_smooth()+
  theme_light()+
  stat_cor(data=cortable, method = "pearson",
           digits = 3,
           label.y= 8)+
  labs(x="EILncSig Score", 
       y="EMT score", 
       title="Correlation of EMT score and EILncSig Score\nTuan et al.'s EMT signature") #Tuan et al.'s EMT signature
#4.5*4.5
