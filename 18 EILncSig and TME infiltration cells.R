################################################################################
library(pals)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(pheatmap)
library(ggcorrplot)

# Read TME file - txt file output via CIBERSORTx online tool
TME.results <-read.table("CIBERSORTx_Job2_tpm.Adjusted.txt",sep = "\t",
                         header = T,row.names = 1,check.names = F)
TME.results = as.data.frame((TME.results))
TME.results = TME.results[,1:22]
rownames(TME.results) = substr(rownames(TME.results),1,12)

# Read EILncSig level phenotype file of TCGA
pheno.data <- read.table("Risktable_tcga.sarc_0.4668533.txt", header = T, sep = "\t")
pheno.data = as.data.frame(pheno.data[,c(10,11)])
pheno.data$ID = rownames(pheno.data)

#
pheno.data = pheno.data[row.names(pheno.data) %in% as.vector(row.names(TME.results)),]
TME.results = TME.results[row.names(TME.results) %in% as.vector(row.names(pheno.data)),]

# Merge TME file with phenotype data
sample.list <- as.character(unique(pheno.data$ID))
pheno.data <-
  pheno.data[match(sample.list, pheno.data$ID), ]
TME.results <-
  TME.results[match(sample.list, rownames(TME.results)), ]
TME.data <- cbind(pheno.data, TME.results)

# Show TME cells
TME.cells <- colnames(TME.results)
TME.cells

#Cell componment boxplot
plot.info <- NULL
for (i in 1:length(TME.cells)) {
  idx.sub <- which(colnames(TME.data) == TME.cells[i])
  sub <- data.frame(
    PATIENT_ID = TME.data$ID,
    Cell_Type = TME.cells[i],
    EILncSig_Level = TME.data$EILncSig_Level,
    Composition = TME.data[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}

#plot stack
ggbarplot(
  plot.info,
  x = "EILncSig_Level",
  y = "Composition",
  size = 0,
  fill = "Cell_Type",
  color = "Cell_Type",
  ylab = "CIBERSORTx abs.score",
  title = "TME cells infiltration between EILncSig levels",
  xlab = ""
) +
  theme_light() + theme(axis.text.x = element_text(size = 11)) +
  scale_fill_manual(values=as.vector(polychrome(26)[1:22]))+
  scale_colour_manual(values=as.vector(polychrome(26)[1:22]))+
  theme(legend.position = NULL)
# pdf 7.5-5.5


#plot bar
table(plot.info$Cell_Type)
cellskeep1 = c("Macrophages M0","Macrophages M1","Macrophages M2")
cellskeep2 = c("B cells memory","B cells naive","NK cells activated",
               "NK cells resting","Dendritic cells activated","Dendritic cells resting",
               "Eosinophils","Monocytes","Neutrophils","Plasma cells",
               "Mast cells activated","Mast cells resting",
               "T cells CD4 memory activated","T cells CD4 memory resting","T cells CD4 naive",
               "T cells CD8","T cells follicular helper","T cells gamma delta",
               "T cells regulatory (Tregs)")

plot.info1 = plot.info[plot.info$Cell_Type %in% cellskeep1,]
plot.info2 = plot.info[plot.info$Cell_Type %in% cellskeep2,]


#1 pdf 4*4
ggboxplot(
  plot.info1,
  x = "Cell_Type",
  y = "Composition",
  color = "black",
  fill = "EILncSig_Level",
  ylab = "CIBERSORTx abs.score",
  title = "TME cells infiltration",
  xlab = "",
  outlier.shape = NA,
  main = ""
) +
  stat_compare_means(
    aes(group=EILncSig_Level),
    method = "wilcox.test", #kruskal.test wilcox.test
    size = 5,  ##
    label = "p.signif", ###p
    #label = "p.format", 
    label.y = 1.39,
    #ref.group = "cluster 2",# Pairwise comparison against all
    hide.ns = T
  )+scale_fill_manual(values=c("red3","royalblue"))+
  scale_colour_manual(values=c("red3","royalblue"))+
  theme_light() + theme(axis.text.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(
    angle = 30,
    hjust = 1,
    vjust = 1
  )) + coord_cartesian(ylim=c(0, 1.5))

#2 pdf 12.5*4
ggboxplot(
  plot.info2,
  x = "Cell_Type",
  y = "Composition",
  color = "black",
  fill = "EILncSig_Level",
  ylab = "CIBERSORTx abs.score",
  title = "TME cells infiltration",
  xlab = "",
  outlier.shape = NA,
  main = ""
) +
  stat_compare_means(
    aes(group=EILncSig_Level),
    method = "wilcox.test", #kruskal.test wilcox.test
    size = 5,  ##
    label = "p.signif", ###p
    #label = "p.format", 
    label.y = 0.5,
    #ref.group = "cluster 2",# Pairwise comparison against all
    hide.ns = T
  )+scale_fill_manual(values=c("red3","royalblue"))+
  scale_colour_manual(values=c("red3","royalblue"))+
  theme_light() + theme(axis.text.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(
    angle = 30,
    hjust = 1,
    vjust = 1
  )) + coord_cartesian(ylim=c(0, 0.515))


#-----------------------------------------------------------------------------------------
#correlation

colnames(TME.data)
# take NK cells activated as an example
ggplot(TME.data, aes(x=EILncSig_Score, y=`NK cells activated`)) + 
  geom_point(color="black",size = 1.0) + 
  geom_smooth(method="lm", se=T)+
  # geom_smooth()+
  theme_light()+
  stat_cor(data=TME.data, 
           method = "spearman",
           digits = 3,
           label.y= 0.375)+
  labs(x="EILncSig Score", 
       y="CIBERSORTx abs.score", 
       title="TME cells infiltration \nNK cells activated")
#3.5-3.5

