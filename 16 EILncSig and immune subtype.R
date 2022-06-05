library(ggplot2)
library(ggpubr)
riskscore = read.table("Risktable_tcga.sarc_0.4668533.txt",header = T,sep = "\t",check.names = F)
riskscore$ID = rownames(riskscore)
# TCGA-SARC immune subtype - obtianed from https://doi.org/10.1016/j.immuni.2018.03.023
# prepare a txt file TCGAsarc immune type.txt
imtype = read.table("TCGAsarc immune type.txt",header = T,sep = "\t")
risktable = merge(riskscore,imtype,by="ID")
colnames(risktable)[14] = "Immune Subtype" 
risktable = risktable[!is.na(risktable$`Immune Subtype`),]

ggplot(risktable, 
       aes(x = `Immune Subtype`,
           y = EILncSig_Score,
           fill = `Immune Subtype`)) +
  geom_boxplot(aes(group= `Immune Subtype`), 
               width=0.65, position=position_dodge(width = 1),
               outlier.shape=NA) +
  scale_fill_manual(values = c("brown3","salmon","darkorange","forestgreen","steelblue")) +
  stat_compare_means(#comparisons = my_comparisons,
    #method = "wilcox.test",
    label = "p.format",
    label.y = 1.45)+
  theme_light() + xlab(NULL) +  ylab("EILncSig Score")+
  ggtitle("Association of Immune Subtype and EILncSig Score \nTCGA-SARC") +
  theme(axis.text.x = element_blank()) + coord_cartesian(ylim=c(0, 1.50))
#6.8*3.75

library(ggstatsplot)
library(pals)

pdf("immune sub2.pdf",width = 9,height = 3.85)
ggbarstats(risktable,
           EILncSig_Level, 
           `Immune Subtype`,
           bar.proptest = T, 
           label = "none",
           palette = 'Set2',#
           results.subtitle = T,
           type = "nonparametric",
           proportion.test = F
) + labs(x = '',y = '',title="Association of EILncSig level and Immune Subtype") + theme_light()+
  theme(axis.text.x = element_text(size = 11.5)) +
  # scale_fill_manual(values = c("skyblue","purple","brown","darkorange"))+
  scale_fill_manual(values = c("royalblue","red2")) +
  coord_flip() #
dev.off()


pdf("immune sub3.pdf",width = 9.5,height = 3.5)
ggbarstats(risktable,
           `Immune Subtype`,
           EILncSig_Level, 
           bar.proptest = T, 
           label = "none",
           palette = 'Set2',#
           results.subtitle = T,
           type = "nonparametric",
           proportion.test = F
) + labs(x = '',y = '',title="Association of EILncSig level and Immune Subtype") + theme_light()+
  theme(axis.text.x = element_text(size = 11.5)) +
  scale_fill_manual(values = c("steelblue","forestgreen","darkorange","salmon","brown3"))+
  # scale_fill_manual(values = c("royalblue","red2")) +
  coord_flip() #
dev.off()
#------------

colnames(risktable)
colnames(risktable)[8]="WWP1-AS1" #official symbol
colnames(risktable)[10]="AFTPH-DT" #official symbol


ggplot(risktable, 
       aes(x = `Immune Subtype`,
           y = `MIR22HG`)) +
  geom_violin(aes(fill=`Immune Subtype`), trim = FALSE) + 
  geom_boxplot(width=0.1) +
  scale_fill_manual(values = c("brown3","salmon","darkorange","forestgreen","steelblue")) +
  stat_compare_means(#comparisons = my_comparisons,
    method = "kruskal.test",
    label = "p.format",
    label.y = 13.5
    )+
  theme_light() + xlab(NULL) +  ylab("Normalized expression")+
  ggtitle("Association of Immune Subtype and EILncSig expression \nMIR22HG") +
  theme(axis.text.x = element_blank()) + coord_cartesian(ylim=c(4.8, 13.72))
#6.8*3.5