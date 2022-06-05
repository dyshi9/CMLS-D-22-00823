load("pan.sarcoma.array.symbol.n1085.sva.rdata")
load("EMTgenes.rdata")
#-----------------------------------------------------------------
#prepare a phenotype txt file, which describe the phentype of an EMT related gene - Epithelial or mesenchymal
pheno = read.table("PMID25214461 pheno.txt",
                    header = T,sep = "\t")
pheno = read.table("PMID29346386 pheno.txt",
                    header = T,sep = "\t")
EMTsig1 = EMTsig[EMTsig$symbol %in% c(PMID29346386$V1),] # PMID25214461
EMTsig1 = merge(EMTsig1,pheno,by.x="symbol",by.y="id")
table(pheno$pheno)
#-----------------------------------------------------------------

expr = pan.sarcoma.array.symbol.n1085.sva[EMTsig1$true.symbol,]
expr = as.data.frame(t(expr))
expr$ID = row.names(expr)

rt=read.table("EMT_cluster.txt",header=T,sep="\t",check.names=F)
colnames(rt)[1] = "ID"
rownames(rt) = 1:nrow(rt)

expr = merge(rt,expr,by.x = "ID",by.y="ID")
expr$cluster = paste0("C",expr$cluster)
colnames(expr)[2] = "EMT_Cluster"

gene.e = as.data.frame(colnames(expr))
gene.e$p = colnames(expr)%in%EMTsig1[EMTsig1$pheno == "Epithelial","true.symbol"]
gene.e = as.vector(gene.e[gene.e$p==T,1])
expr.e = expr[,c("ID","EMT_Cluster",gene.e)]

gene.m = as.data.frame(colnames(expr))
gene.m$p = colnames(expr)%in%EMTsig1[EMTsig1$pheno == "Mesenchymal","true.symbol"]
gene.m = as.vector(gene.m[gene.m$p==T,1])
expr.m = expr[,c("ID","EMT_Cluster",gene.m)]

expr.e$score.e = apply(expr.e[,3:ncol(expr.e)],1,mean)
expr.m$score.m = apply(expr.m[,3:ncol(expr.m)],1,mean)

expr.e = expr.e[,c("ID","EMT_Cluster","score.e")]
expr.m = expr.m[,c("ID","EMT_Cluster","score.m")]

expr.score = cbind(expr.e,expr.m)
expr.score$emt.score = expr.score$score.m - expr.score$score.e
expr.score = expr.score[,c(1,2,3,6,7)]


library(ggplot2)
library(ggpubr)
ggplot(expr.score, 
       aes(x = EMT_Cluster,
           y = emt.score)) +
  geom_violin(aes(fill=EMT_Cluster), trim = FALSE) + 
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = c("royalblue","red3"))+
  stat_compare_means(
    label = "p.signif",
    method = "wilcox.test",
    hide.ns = T,
    size = 5.5, 
    label.y = 3.59)  + 
  theme_light() + xlab(NULL) + 
  ylab("EMT score") + ggtitle("Hollern et al.'s EMT signature") +
  theme(axis.text.x = element_text(size = 0))

library(ggplot2)
library(ggpubr)
ggplot(expr.score, 
       aes(x = EMT_Cluster,
           y = score.m)) +
  geom_violin(aes(fill=EMT_Cluster), trim = FALSE) + 
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = c("steelblue","salmon"))+
  stat_compare_means(
    label = "p.signif",
    method = "wilcox.test",
    hide.ns = T,
    size = 5.5, 
    label.y = 8.85)  + 
  theme_light() + xlab(NULL) + 
  ylab("Average expression of Mesenchymal markers") + ggtitle("Hollern et al.'s EMT signature") +
  theme(axis.text.x = element_text(size = 0))

library(ggplot2)
library(ggpubr)
ggplot(expr.score, 
       aes(x = EMT_Cluster,
           y = score.e)) +
  geom_violin(aes(fill=EMT_Cluster), trim = FALSE) + 
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = c("steelblue","salmon"))+
  stat_compare_means(
    label = "p.format",
    method = "wilcox.test",
    hide.ns = T,
    size = 4.0, 
    label.y = 7.58)  + 
  theme_light() + xlab(NULL) + 
  ylab("Average expression of Epithelial markers") + ggtitle("Hollern et al.'s EMT signature") +
  theme(axis.text.x = element_text(size = 0))



