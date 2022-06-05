#load matrix row=gene col=sample
load("pan.sarcoma.array.symbol.n1085.sva.rdata")
expr = pan.sarcoma.array.symbol.n1085.sva
library("org.Hs.eg.db") 
library("AnnotationDbi")
genes=as.vector(row.names(expr))
entrezIDs = mapIds(org.Hs.eg.db,
                   keys= genes, #Column containing keytype
                   column="ENTREZID",
                   keytype="SYMBOL",
                   multiVals="first")
entrezIDs <- as.character(entrezIDs)
out = cbind(expr,entrezID=entrezIDs)
out = as.data.frame(out)
out = out[!is.na(out$entrezID),]
row.names(out) = out[,ncol(out)]
out=out[,-ncol(out)]
#prepare the normalized counts matrix (rownames are ENTREZIDs of genes) 
expr_for_GSVA = out
save(expr_for_GSVA,file = "expr_for_GSVA.rdata")

library(GSVA) 
library(limma) 
library(pheatmap)
load("expr_for_GSVA.rdata")
#prepare the database imported from .gmt file. 

#take GO_BP as an example
original_gmt_GSVA <- readLines('c5.go.bp.v7.3.entrez.gmt') 
strsplit_no_name <- function(gmt.list_layer){ 
  as.character(unlist(strsplit(gmt.list_layer, split = '\t', fixed = T)))[-2]} 
database_list_GSVA <- lapply(original_gmt_GSVA, strsplit_no_name) 

for (layers in 1:length(database_list_GSVA)) { 
  names(database_list_GSVA)[layers] <- database_list_GSVA[layers][[1]][1] 
  database_list_GSVA[layers][[1]] <- database_list_GSVA[layers][[1]][-1]}
rm(layers,original_gmt_GSVA,strsplit_no_name)
#the input data must be in matrix format. 
#We should set the distribution as 'Poisson' when using RNA-seq integer counts data rather than the default 'Gaussian'. 
#If the matrix is derived from microarray or the RNA-seq integer counts have been transformed to log2 form or..., 
#we can directly set kcdf as 'Gaussian'.
es <- gsva(as.matrix(expr_for_GSVA), 
           database_list_GSVA, 
           mx.diff=TRUE, verbose=TRUE, 
           method='gsva', kcdf='Gaussian', parallel.sz=8)
write.table(es,"GSVA_go.bp.txt",row.names = T,sep="\t")

es.1 = read.table("GSVA_go.bp.txt",header = T,sep = "\t",check.names = F)
es.all = es.1
save(es.all,file = "es.all.rdata")

library(effectsize)
# make sample in row, variable in column
es_t = as.data.frame(t(es.all))
es_t_z <- standardize(es_t, robust = TRUE) 
es.all.stad = as.data.frame(t(es_t_z))
save(es.all.stad,file = "es.all.stad.rdata")


library(limma)
condition_table = read.table("EMT_cluster.txt",header = T,check.names = F,sep = "\t")
colnames(condition_table)[1] = "ID"
condition_table$cluster = paste0("C",condition_table$cluster)
condition_table_for_limma <- model.matrix(~0 + condition_table$cluster) 
colnames(condition_table_for_limma) <- levels(as.factor(condition_table$cluster))
rownames(condition_table_for_limma) <- condition_table$ID 
contrast_matrix <- makeContrasts('C2-C1',levels = as.factor(condition_table$cluster)) 
#generate the contrast matrix.
#It means that genes in normal group will be at denominator level. 
es.limma = es.all
lmfit <- lmFit(es.limma, condition_table_for_limma) 
lmfit <- contrasts.fit(lmfit, contrast_matrix) 
lmfit <- eBayes(lmfit) 
es.pathway = topTable(lmfit,coef=T, n=Inf,adjust="BH", sort.by = "p")
es.pathway.sig = es.pathway[(es.pathway$adj.P.Val <0.05)&(abs(es.pathway$logFC) >= 0.2),] #non-zscore
write.table(es.pathway.sig,"es.gobp.pathway.txt",sep = "\t",row.names = T,quote = F)

#plot heatmap
#choose pathway to display - select from es.XXXXX.txt and re-save as needed pathways.txt
path = read.table("needed pathways.txt",sep = "\t",header = T)
#read all es
es.1 = read.table("./gsva scores/GSVA_go.bp.txt",header = T,sep = "\t",check.names = F)
es.2 = read.table("./gsva scores/GSVA_wiki.txt",header = T,sep = "\t",check.names = F)
es.3 = read.table("./gsva scores/GSVA_kegg.txt",header = T,sep = "\t",check.names = F)
es.4 = read.table("./gsva scores/GSVA_react.txt",header = T,sep = "\t",check.names = F)

#sort
es.all = rbind(es.1,es.2,es.3,es.4)
rm(es.1,es.2,es.3,es.4)
es.all = es.all[path$ID,]
rownames(es.all) = path$name

#make annotation - cluster
emtc = read.table("EMT_cluster.txt",header = T,check.names = F,sep = "\t")
emtc = emtc[order(emtc$cluster),]
emt.c = as.data.frame(emtc$cluster)
rownames(emt.c) = rownames(emtc)
colnames(emt.c) = "EMT_Cluster"
emt.c$EMT_Cluster = gsub("1","C1",emt.c$EMT_Cluster)
emt.c$EMT_Cluster = gsub("2","C2",emt.c$EMT_Cluster)
rm(emtc)

# annotation - sample datasets
load("sample.batch.n1085.rdata")
sample.batch.n1085 = sample.batch.n1085[,c(1,2)]
rownames(sample.batch.n1085) = sample.batch.n1085$ID
sample.batch.n1085 = sample.batch.n1085[row.names(emt.c),] #- sort by emc.c
sample.batch = as.data.frame(sample.batch.n1085$DataSet_ID)
rownames(sample.batch) = sample.batch.n1085$ID
colnames(sample.batch) = "Datasets"
rm(sample.batch.n1085)

emt.c_datasets = cbind(emt.c,sample.batch)

#make annotation - gene sets
gene.annot = path[,c(3,2)]
gene_annot = as.data.frame(gene.annot$Gene_Sets)
rownames(gene_annot) = gene.annot$name
colnames(gene_annot) = "Genesets"
rm(gene.annot)

#sort es.all
es.all = es.all[,row.names(emt.c)]


#plot
library(pheatmap)

datasetscol = as.character(polychrome(17))

my_colour = list(
  EMT_Cluster = c(C1 = "steelblue", C2 = "salmon"),
  Datasets = c(`E-MEXP-1922 `= datasetscol[1],`E-MEXP-3628`= datasetscol[2],`E-MEXP-964`= "darkred",
               `E-TABM-1202`= "brown1",GSE13433  = "tan4",GSE142162 = datasetscol[6],
               GSE14827  = datasetscol[7],GSE17618 = datasetscol[8],GSE20196 = "maroon",
               GSE20559 = datasetscol[10],GSE23980 = datasetscol[11],GSE34620 = datasetscol[12],
               GSE34800 = datasetscol[13],GSE37371 = datasetscol[14],GSE66533 = datasetscol[15],
               GSE71118 = datasetscol[16],GSE87437 = datasetscol[17]),
  Genesets =c(`GO Biological Process` = "#3182BD",KEGG="#E6550D",
              REACTOME="#31A354",WikiPathways="#756BB1"))


pdf(file="heatmap GSVA.pdf",height=8,width=11)
bk = unique(c(seq(-1.75,1.75, length= 100)))
pheatmap(es.all, 
         annotation_col = emt.c_datasets, 
         annotation_row = gene_annot,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         breaks = bk,
         cluster_cols = F,
         cluster_rows = T,
         fontsize=9,
         fontsize_row=7,
         fontsize_col=6.5,
         scale="row",
         gaps_col = 636,
         #cutree_cols = 2,
         show_colnames=F,
         show_rownames=T,
         main = "GSVA Enrichment score",
         annotation_colors = my_colour
)
dev.off()
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

################################################################################

# Evaluation of TME infiltrating cells

# Cibersortx results were obtained by using the CibersortX online tool.
# CIBERSORTx_pan1085_Adjusted_origninal.txt

library(pals)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(pheatmap)
library(ggcorrplot)

# Read TME file
TME.results <-read.table("CIBERSORTx_pan1085_Adjusted_origninal.txt",sep = "\t",
                         header = T,row.names = 1,check.names = F)
TME.results = as.data.frame((TME.results))
TME.results = TME.results[,1:22]
# Read phenotype file of TCGA
pheno.data <- read.table("EMT_cluster.txt", header = T, sep = "\t")
colnames(pheno.data)[1] = "ID"
colnames(pheno.data)[2] = "EMT_Cluster"
pheno.data$EMT_Cluster = paste0("C",pheno.data$EMT_Cluster)

#取交集
pheno.data = pheno.data[row.names(pheno.data) %in% as.vector(row.names(TME.results)),]
TME.results = TME.results[row.names(TME.results) %in% as.vector(row.names(pheno.data)),]
pheno.data$PATIENT_ID = row.names(pheno.data)

# Merge TME file with phenotype data
sample.list <- as.character(unique(pheno.data$PATIENT_ID))
pheno.data <-
  pheno.data[match(sample.list, pheno.data$PATIENT_ID), ]
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
    PATIENT_ID = TME.data$PATIENT_ID,
    Cell_Type = TME.cells[i],
    EMT_Cluster = TME.data$EMT_Cluster,
    Composition = TME.data[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}

#plot stack
ggbarplot(
  plot.info,
  x = "EMT_Cluster",
  y = "Composition",
  size = 0,
  fill = "Cell_Type",
  color = "Cell_Type",
  ylab = "CIBERSORTx abs.score",
  title = "TME cells infiltration",
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
               "NK cells resting","Dendritic cells activated","Dendritic cells resting")
cellskeep3 = c("Eosinophils","Monocytes","Neutrophils","Plasma cells",
               "Mast cells activated","Mast cells resting")
cellskeep4 = c("T cells CD4 memory activated","T cells CD4 memory resting","T cells CD4 naive",
               "T cells CD8","T cells follicular helper","T cells gamma delta",
               "T cells regulatory (Tregs)")

plot.info1 = plot.info[plot.info$Cell_Type %in% cellskeep1,]
plot.info2 = plot.info[plot.info$Cell_Type %in% cellskeep2,]
plot.info3 = plot.info[plot.info$Cell_Type %in% cellskeep3,]
plot.info4 = plot.info[plot.info$Cell_Type %in% cellskeep4,]

#1 pdf 5.5*3.0
ggboxplot(
  plot.info1,
  x = "Cell_Type",
  y = "Composition",
  color = "black",
  fill = "EMT_Cluster",
  ylab = "CIBERSORTx abs.score",
  title = "TME cells infiltration",
  xlab = "",
  outlier.shape = NA,
  main = ""
) +
  stat_compare_means(
    aes(group=EMT_Cluster),
    method = "wilcox.test", #kruskal.test wilcox.test
    size = 5,  ##字体大小
    label = "p.signif", ###p值格式
    #label = "p.format", 
    #label.y = 19.05,
    #ref.group = "cluster 2",# Pairwise comparison against all
    hide.ns = T
  )+scale_fill_manual(values=c("steelblue","salmon"))+
  scale_colour_manual(values=c("steelblue","salmon"))+
  theme_light() + theme(axis.text.x = element_text(size = 10.5)) +
  theme(axis.text.x = element_text(
    angle = 0,
    hjust = 1,
    vjust = 1
  )) + coord_flip(ylim=c(0, 1.05)) 

#2 pdf 6.5*4.25
ggboxplot(
  plot.info2,
  x = "Cell_Type",
  y = "Composition",
  color = "black",
  fill = "EMT_Cluster",
  ylab = "CIBERSORTx abs.score",
  title = "TME cells infiltration",
  xlab = "",
  outlier.shape = NA,
  main = ""
) +
  stat_compare_means(
    aes(group=EMT_Cluster),
    method = "wilcox.test", #kruskal.test wilcox.test
    size = 5,  ##字体大小
    label = "p.signif", ###p值格式
    #label = "p.format", 
    label.y = 0.125,
    #ref.group = "cluster 2",# Pairwise comparison against all
    hide.ns = T
  )+scale_fill_manual(values=c("steelblue","salmon"))+
  scale_colour_manual(values=c("steelblue","salmon"))+
  theme_light() + theme(axis.text.x = element_text(size = 10.5)) +
  theme(axis.text.x = element_text(
    angle = 0,
    hjust = 1,
    vjust = 1
  )) + coord_flip(ylim = c(0,0.135))

#3 pdf 6.5*4.25
ggboxplot(
  plot.info3,
  x = "Cell_Type",
  y = "Composition",
  color = "black",
  fill = "EMT_Cluster",
  ylab = "CIBERSORTx abs.score",
  title = "TME cells infiltration",
  xlab = "",
  outlier.shape = NA,
  main = ""
) +
  stat_compare_means(
    aes(group=EMT_Cluster),
    method = "wilcox.test", #kruskal.test wilcox.test
    size = 5,  ##字体大小
    label = "p.signif", ###p值格式
    #label = "p.format", 
    label.y = 0.2,
    #ref.group = "cluster 2",# Pairwise comparison against all
    hide.ns = T
  )+scale_fill_manual(values=c("steelblue","salmon"))+
  scale_colour_manual(values=c("steelblue","salmon"))+
  theme_light() + theme(axis.text.x = element_text(size = 10.5)) +
  theme(axis.text.x = element_text(
    angle = 0,
    hjust = 1,
    vjust = 1
  )) + coord_flip(ylim = c(0,0.205))


#4 pdf 6.5*4.75
ggboxplot(
  plot.info4,
  x = "Cell_Type",
  y = "Composition",
  color = "black",
  fill = "EMT_Cluster",
  ylab = "CIBERSORTx abs.score",
  title = "TME cells infiltration",
  xlab = "",
  outlier.shape = NA,
  main = ""
) +
  stat_compare_means(
    aes(group=EMT_Cluster),
    method = "wilcox.test", #kruskal.test wilcox.test
    size = 5,  ##字体大小
    label = "p.signif", ###p值格式
    #label = "p.format", 
    label.y = 0.28,
    #ref.group = "cluster 2",# Pairwise comparison against all
    hide.ns = T
  )+scale_fill_manual(values=c("steelblue","salmon"))+
  scale_colour_manual(values=c("steelblue","salmon"))+
  theme_light() + theme(axis.text.x = element_text(size = 10.5)) +
  theme(axis.text.x = element_text(
    angle = 0,
    hjust = 1,
    vjust = 1
  )) + coord_flip(ylim = c(0,0.2950))

# Read TME file
TME.results <-read.table("CIBERSORTx_pan1085_Adjusted_origninal.txt",sep = "\t",
                         header = T,row.names = 1,check.names = F)
TME.results = as.data.frame((TME.results))
TME.results = TME.results[,1:22]
# Read phenotype file of TCGA
pheno.data <- read.table("EMT_cluster.txt", header = T, sep = "\t")
colnames(pheno.data)[1] = "ID"
colnames(pheno.data)[2] = "EMT_Cluster"
pheno.data$EMT_Cluster = paste0("C",pheno.data$EMT_Cluster)

#取交集
pheno.data = pheno.data[row.names(pheno.data) %in% as.vector(row.names(TME.results)),]
TME.results = TME.results[row.names(TME.results) %in% as.vector(row.names(pheno.data)),]
pheno.data$PATIENT_ID = row.names(pheno.data)

# Merge TME file with phenotype data
sample.list <- as.character(unique(pheno.data$PATIENT_ID))
pheno.data <-
  pheno.data[match(sample.list, pheno.data$PATIENT_ID), ]
TME.results <-
  TME.results[match(sample.list, rownames(TME.results)), ]
TME.data <- cbind(pheno.data, TME.results)


colnames(TME.data) = gsub(" ","_",colnames(TME.data))
colnames(TME.data)[12] = "T_cells_regulatory_Tregs"
colnames(TME.data)

cellname = as.vector(colnames(TME.data))[4:25]

# pval
wilcoxoutall= data.frame()
for (i in cellname)
{
  wilcoxout = wilcox.test(x=TME.data[TME.data$EMT_Cluster == "C1",i],
                          y=TME.data[TME.data$EMT_Cluster == "C2",i],
                          alternative = "two.sided")
  wilcoxoutall = rbind(wilcoxoutall,cbind(id=i,pval=wilcoxout$p.value))
}

# difference
meanoutall= data.frame()
for (i in cellname)
{
  meanout = mean(TME.data[TME.data$EMT_Cluster == "C1",i]) - 
    mean(TME.data[TME.data$EMT_Cluster == "C2",i])
  meanoutall = rbind(meanoutall,cbind(id=i,meanval=meanout))
}

wilcox.res = merge(wilcoxoutall,meanoutall,by="id")
colnames(wilcox.res)

wilcox.res$pval = as.numeric(wilcox.res$pval)
wilcox.res$meanval = as.numeric(wilcox.res$meanval)
wilcox.res$`-log10(pval)` = -log10(wilcox.res$pval)
wilcox.res1 = wilcox.res[!wilcox.res$id %in% c("Macrophages_M2",
                                               "Macrophages_M0",
                                               "Macrophages_M1"),]
wilcox.res1 = wilcox.res1[wilcox.res1$pval <0.05,]
wilcox.res2 = wilcox.res[wilcox.res$id %in% c("Macrophages_M2",
                                              "Macrophages_M0",
                                              "Macrophages_M1"),]
wilcox.res2 = wilcox.res2[wilcox.res2$pval <0.05,]


#plot 6*4.0
library(ggplot2)
ggplot(data = wilcox.res1, 
       aes(x = reorder(id, -meanval), 
           y = meanval)) +
  geom_bar(stat = "identity", aes(fill = `-log10(pval)`)) +
  theme_light() + coord_flip() +
  ylab("Difference of CIBERSORTx scores") + 
  xlab("")+
  ggtitle("EMT Cluster: C1 vs C2") +
  scale_fill_gradient2(
    low ="royalblue",
    mid = "purple",
    high ="red3",
    midpoint = 12,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  )


################################################################################

# Evaluation of TME infiltrating cells

# xCell

library(immunedeconv)

#a matrix with genes in rows and samples in columns
load("pan.sarcoma.array.symbol.n1085.sva.rdata")

res.xcell = deconvolute_xcell(pan.sarcoma.array.symbol.n1085.sva,
                              arrays = T,
                              rnaseq = F,
                              cell.types.use = NULL,
                              scale = TRUE)
res.xcell = t(res.xcell)
write.table(res.xcell,"res.xcell.txt",sep = "\t",row.names = T,quote = F)

#----------------------
#plot
library(pals)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(pheatmap)
library(ggcorrplot)

# Read TME file
TME.results <-read.table("res.xcell.txt",sep = "\t",
                         header = T,row.names = 1,check.names = F)
TME.results = as.data.frame((TME.results))
TME.results = TME.results[,65:67]
# Read phenotype file of TCGA
pheno.data <- read.table("EMT_cluster.txt", header = T, sep = "\t")
colnames(pheno.data)[1] = "ID"
colnames(pheno.data)[2] = "EMT_Cluster"
pheno.data$EMT_Cluster = paste0("C",pheno.data$EMT_Cluster)


TME.results$ID = rownames(TME.results)
TME.data = merge(TME.results,pheno.data,by="ID")

colnames(TME.data)


# 1-immuneScore 2-stromaScore 3-MicroenvironmentScore
ggplot(TME.data, 
       aes(x = EMT_Cluster,
           y = MicroenvironmentScore)) +
  geom_violin(aes(fill=EMT_Cluster), trim = FALSE) + 
  geom_boxplot(width = 0.1)+
  scale_fill_manual(values = c("salmon","steelblue"))+
  stat_compare_means(
    label = "p.signif",
    method = "wilcox.test",
    hide.ns = T,
    size = 4.75, 
    label.y =1.855)+
  theme_light() + xlab(NULL) + 
  ylab("xCell score") + ggtitle("xCell - Microenvironment") +
  theme(axis.text.x = element_text(size = 0))

