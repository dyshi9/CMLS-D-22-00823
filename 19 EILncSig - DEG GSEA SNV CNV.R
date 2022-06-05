# DEG-----------------------------------------------------------------------------------------
Type.sample=read.table("Risktable_tcga.sarc_0.4668533.txt",sep="\t",header=T,check.names=F)
Type.sample$ID = rownames(Type.sample)
Type.sample=Type.sample[,c(12,11,10)]
Type.sample=Type.sample[order(Type.sample$EILncSig_Score),]
#load RNA-seq expr 
load("TCGA.SARC_n259_count.symbol.rdata")
data = TCGA.SARC_n259_count.symbol
colnames(data) <- substr(colnames(data), 1, 12)
data = as.matrix(data)
data = data[,colnames(data) %in% as.vector(Type.sample$ID)]
expr = data[,Type.sample$ID]
library(DESeq2)
expr = round(expr)
expr = expr[rowMeans(expr)>0,]
head(expr)
dim(expr)
table(Type.sample$EILncSig_Level)
group <- c(rep('Low_Risk',152),rep('High_Risk',107))
condition = factor(group)
sample <- data.frame(row.names = colnames(expr), condition)
ddsCountTable <- DESeqDataSetFromMatrix(countData = expr,
                                        colData = sample,
                                        design= ~ condition)
ddsCountTable$condition<- relevel(ddsCountTable$condition, 
                                  ref = "Low_Risk") # 
ddsCountTable <- estimateSizeFactors(ddsCountTable)
keep <- rowSums(counts(ddsCountTable) >= 5) >= 3
ddsCountTable <- ddsCountTable[keep,]
dim(ddsCountTable)
dds <- DESeq(ddsCountTable)
normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)
save(normalized_counts, file="normalized_counts.sarc259.risklevel.rdata")
vsd1 <- vst(dds,blind=FALSE)  
vsd <- assay(vsd1)
save(vsd, file="VST.sarc259.risklevel.rdata")

contrast_Group <- c("condition","High_Risk","Low_Risk")
res <- results(dds,  contrast=contrast_Group)

baseHigh_Risk <- counts(dds, normalized=TRUE)[,colData(dds)$condition == "High_Risk"]
baseMeanHigh_Risk <- as.data.frame(rowMeans(baseHigh_Risk))
colnames(baseMeanHigh_Risk) <- "High_Risk"
head(baseMeanHigh_Risk)

baseLow_Risk <- counts(dds, normalized=TRUE)[,colData(dds)$condition == "Low_Risk"]
baseMeanLow_Risk <- as.data.frame(rowMeans(baseLow_Risk))
colnames(baseMeanLow_Risk) <- "Low_Risk"
head(baseMeanLow_Risk)

res <- cbind(baseMeanHigh_Risk, baseMeanLow_Risk, as.data.frame(res))
head(res)

res <- cbind(ID=rownames(res), as.data.frame(res))
res$padj[is.na(res$padj)] <- 1
res <- res[order(res$pvalue),]
head(res)
res$regulation = as.factor(ifelse(res$padj < 0.05 & abs(res$log2FoldChange) >= 1,
                                  ifelse(res$log2FoldChange >= 1,'Up','Down'),'No Significance'))
table(res$regulation)
write.table(as.data.frame(res), file="DESeq2.diff_sarc_EILncSig.txt", sep="\t", quote=F, row.names=F)

library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 1.50,
                labSize = 2.0,
                xlim = c(min(res$log2FoldChange), max(res$log2FoldChange)),
                ylim = c(0, max(-log10(res$padj))),
                title = 'DEG: EILncSig High-Risk vs Low-Risk',
                subtitle = "",
                caption = '',
                col = c('grey30', 'forestgreen', 'royalblue', 'firebrick1'),
                colAlpha = 0.6
)
#8.0-7.5
load("EMTgenes.rdata")
EMTsig1 = EMTsig
resemt = res[res$ID %in% as.vector(EMTsig1$true.symbol),]
table(resemt$regulation)
write.table(as.data.frame(resemt), file="DESeq2.diff_sarc.emt_EILncSig.txt", sep="\t", quote=F, row.names=F)
EnhancedVolcano(resemt,
                lab = rownames(resemt),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 2.5,
                labSize = 3.2,
                xlim = c(min(resemt$log2FoldChange), max(resemt$log2FoldChange)),
                ylim = c(0, max(-log10(resemt$padj))),
                title = 'DEG: EILncSig High-Risk vs Low-Risk',
                subtitle = "EMT Signature",
                caption = '',
                col = c('grey30', 'forestgreen', 'royalblue', 'firebrick1'),
                colAlpha = 0.6
)

#plot heatmap
# load expr
expr.heatmap = log2(normalized_counts + 1)
res = read.table("DESeq2.diff_sarc_EILncSig.txt",header = T,sep = "\t",check.names = F)
#choose genes
load("EMTgenes.rdata")
EMTsig1 = EMTsig
resemt = res[res$ID %in% as.vector(EMTsig1$true.symbol),]
table(resemt$regulation)
resemt = resemt[resemt$regulation != "No Significance",]
#read EIlncSig
clin1 = read.table("Risktable_tcga.sarc_0.4668533.txt",sep = "\t",header = T,check.names = F)
clin1$sample_id = rownames(clin1)
clin2 = read.table("SARC_n261_210914.txt",sep = "\t",header = T,check.names = F)
clin= merge(clin1,clin2,by="sample_id")
colnames(clin)[1] = "Tumor_Sample_Barcode" # 
clin = clin[,c(1,3,11,12,14,25)]
colnames(clin)[2] = "Virtual_Status"
colnames(clin)[6] = "Metastatic_Disease" 
clin$Virtual_Status = gsub("0","Alive",clin$Virtual_Status)
clin$Virtual_Status = gsub("1","Dead",clin$Virtual_Status)
clin[is.na(clin$Metastatic_Disease),"Metastatic_Disease"] = "N/A"
clin=clin[order(clin$EILncSig_Score),]
immune.type =read.table("TCGAsarc immune type.txt",header = T,sep = "\t")
clin$immune.type  = NA
for (i in clin$Tumor_Sample_Barcode)
{
  if (i %in% immune.type$ID)
  {
    clin[clin$Tumor_Sample_Barcode == i,7] = immune.type[immune.type$ID == i,3]
}}
clin[is.na(clin$immune.type),"immune.type"] = "N/A"
tme.type =read.table("tME_cluster4_annotated.txt",header = T,sep = "\t")
clin$tme.type  = NA
for (i in clin$Tumor_Sample_Barcode)
{
  if (i %in% tme.type$ID)
  {
    clin[clin$Tumor_Sample_Barcode == i,8] = tme.type[tme.type$ID == i,2]
  }}
clin[is.na(clin$tme.type),"tme.type"] = "N/A"
#sort
expr.heatmap = expr.heatmap[,clin$Tumor_Sample_Barcode]
expr.heatmap = expr.heatmap[resemt$ID,]
#make annotation - sample
clin = as.data.frame(clin)
rownames(clin) = clin$Tumor_Sample_Barcode
clin = clin[,c(4,6,7,8)]
colnames(clin)[3] = "Immune_Type"
colnames(clin)[4] = "TME_Pattern"
#make annotation - gene sets
pheno1 = read.table("PMID25214461 pheno.txt",header = T,sep = "\t")
#pheno2 = read.table("PMID26420858 pheno.txt",header = T,sep = "\t") # no detailed information
pheno3 = read.table("PMID28680090 pheno.txt",header = T,sep = "\t")
pheno4 = read.table("PMID29346386 pheno.txt",header = T,sep = "\t")
pheno5 = read.table("PMID31069256 pheno.txt",header = T,sep = "\t")
pheno = rbind(pheno1,pheno3,pheno4,pheno5)
pheno = pheno[!duplicated(pheno$id),]
pheno[!pheno$id %in% EMTsig$true.symbol,"id"] #check
table(resemt$ID %in% pheno$id) #check
resemt[!resemt$ID %in% pheno$id,"ID"]
gene.ann = resemt[,c(1,10)]
gene.ann$pheno  = NA
for (i in gene.ann$ID)
{
  if (i %in% pheno$id)
  {
    gene.ann[gene.ann$ID == i,3] = pheno[pheno$id == i,2]
  }}
gene.ann[is.na(gene.ann$pheno),"pheno"] = "N/A"
table(gene.ann$pheno)
gene.ann.name = gene.ann$ID
gene.ann = gene.ann[,-c(1,2)]
gene.ann = as.data.frame(gene.ann)
rownames(gene.ann) = gene.ann.name
colnames(gene.ann) = "Phenotype"
#plot
library(pheatmap)
library(pals)
colnames(clin)
table(gene.ann$Phenotype)
datasetscol = as.character(polychrome(17))
my_colour = list(
  EILncSig_Level = c(Low  = "royalblue", High  = "red3"),
  Metastatic_Disease = c(No  = "forestgreen", Yes   = "brown", `N/A` = "grey60"),
  Immune_Type = c(`IFN-gamma Dominant`= "orange",`Inflammatory`= "darkgreen",
                  `Lymphocyte Depleted`= "steelblue",`Wound Healing`= "red3",
                  `TGF-beta Dominant`  = "pink",`N/A` = "grey60"),
  TME_Pattern = c(`Immune-Enriched Fibrotic` = "purple", `Depleted` = "darkorange",
                    `Immune-Enriched Non-Fibrotic` = "skyblue",`Fibrotic` = "brown"),
  Phenotype = c(Epithelial = "#3182BD", Mesenchymal = "#E6550D", `N/A` = "grey60")
  )
gene.notNA = cbind(rownames(gene.ann),gene.ann$Phenotype)
gene.notNA = gene.notNA[gene.notNA[,2] != "N/A",]
expr.heatmap.notNA = expr.heatmap[gene.notNA[,1],]
pdf(file="DEG emt.pdf",height=4.25,width=6)
bk = unique(c(seq(-1.75,1.75, length= 100)))
pheatmap(expr.heatmap, 
         annotation_col = clin, 
         annotation_row = gene.ann,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         breaks = bk,
         cluster_cols = F,
         cluster_rows = T,
         fontsize=7,
         scale="row",
         gaps_col = 152,
         show_colnames=F,
         show_rownames=F,
         main = "EILncSig_Level High vs Low (EMT Signature)",
         annotation_colors = my_colour
)
dev.off()

#all genes
res = read.table("DESeq2.diff_sarc_EILncSig.txt",header = T,sep = "\t")
res.sig.up = res[res$regulation == "Up",]
res.sig.down = res[res$regulation == "Down",]
res.sig = rbind(res.sig.up,res.sig.down)
expr.heatmap = log2(normalized_counts + 1)
expr.heatmap = expr.heatmap[,row.names(clin)]
expr.heatmap = expr.heatmap[res.sig$ID,]
my_colour = list(
  EILncSig_Level = c(Low  = "royalblue", High  = "red3"),
  Metastatic_Disease = c(No  = "forestgreen", Yes   = "brown", `N/A` = "grey60"),
  Immune_Type = c(`IFN-gamma Dominant`= "orange",`Inflammatory`= "darkgreen",
                  `Lymphocyte Depleted`= "steelblue",`Wound Healing`= "red3",
                  `TGF-beta Dominant`  = "pink",`N/A` = "grey60"),
  TME_Pattern = c(`Immune-Enriched Fibrotic` = "purple", `Depleted` = "darkorange",
                  `Immune-Enriched Non-Fibrotic` = "skyblue",`Fibrotic` = "brown")
)
clin = as.data.frame(clin)
rownames(clin) = clin$Tumor_Sample_Barcode
clin = clin[,c(4,6,7,8)]
colnames(clin)[3] = "Immune_Type"
colnames(clin)[4] = "TME_Pattern"
pdf(file="DEG all.pdf",height=4.25,width=7.0)
bk = unique(c(seq(-1.75,1.75, length= 100)))
pheatmap(expr.heatmap, 
         annotation_col = clin,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         breaks = bk,
         cluster_cols = F,
         cluster_rows = F,
         fontsize=8,
         scale="row",
         gaps_col = 152,
         show_colnames=F,
         show_rownames=F,
         main = "EILncSig_Level High vs Low",
         annotation_colors = my_colour
)
dev.off()

# GSEA-----------------------------------------------------------------------------------------
#use the DESeq2 result file
#1st colounm-gene symbolㄛ2nd column-logFC/t=statistics
#To entrezID
library("org.Hs.eg.db") 
library("AnnotationDbi")
rt=read.table("DESeq2.diff_sarc_EILncSig.txt",sep="\t",check.names=F,header=T) 
genes=as.vector(rt[,1])
entrezIDs = mapIds(org.Hs.eg.db,
                   keys= genes, #Column containing keytype
                   column="ENTREZID",
                   keytype="SYMBOL",
                   multiVals="first")
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
out = out[!is.na(out$entrezID),]
write.table(out,file="entrezID.txt",sep="\t",quote=F,row.names=F)

library("clusterProfiler")
library("enrichplot")
library("ggplot2")
rt=read.table("entrezID.txt",sep="\t",header=T,check.names=F)   
rt=rt[is.na(rt[,"entrezID"])==F,] 
genelist_input<-rt$stat
names(genelist_input) = as.character(rt$entrezID)
genelist_input = sort(genelist_input, decreasing = TRUE)
head(genelist_input)
#gmt
gmtfile <- "c5.go.bp.v7.3.entrez.gmt" 
gmtlist <- read.gmt(gmtfile)
gmtlist$term = gsub("GOBP_","",gmtlist$term)
#GSEA
GSEA_input = genelist_input
GSEA_input <- sort(GSEA_input, decreasing = T)
gsea_result <- GSEA(GSEA_input,TERM2GENE = gmtlist, pvalueCutoff  = 0.05, pAdjustMethod = "BH", verbose = T)
gsea_resultplus <- setReadable(gsea_result, 'org.Hs.eg.db', keyType="ENTREZID")
write.table(gsea_resultplus,file="gesagobp.txt",sep="\t",quote=F,row.names = F)
edo = gsea_result
# select specific items to be displayed
edo@result = edo@result[c(24,38,185,230,297,339,371,403,472,502,20,3,5,11,18,25,34,39,44,115),]
pdf(file="GSEA GOBP.pdf",width = 10,height = 7)
ridgeplot(edo,
          showCategory = 20,
          fill = "p.adjust",
          core_enrichment = TRUE,
          label_format = 50,
          orderBy = "NES",
          decreasing = FALSE)
dev.off()

gmtfile <- "c2.cp.reactome.v7.3.entrez.gmt"
gmtlist <- read.gmt(gmtfile)
gmtlist$term = gsub("REACTOME_","",gmtlist$term)
GSEA_input = genelist_input
GSEA_input <- sort(GSEA_input, decreasing = T)
gsea_result <- GSEA(GSEA_input,TERM2GENE = gmtlist, pvalueCutoff  = 0.05, pAdjustMethod = "BH", verbose = T)
gsea_resultplus <- setReadable(gsea_result, 'org.Hs.eg.db', keyType="ENTREZID")
write.table(gsea_resultplus,file="gesareact.txt",sep="\t",quote=F,row.names = F)
edo = gsea_result
edo@result = edo@result[c(58,69,93,106,122,127,138,147,201,208,2,6,8,9,11,17,21,32,35,40),]
pdf(file="GSEA Reactome.pdf",width = 9.5,height = 7)
ridgeplot(edo,
          showCategory = 20,
          fill = "p.adjust",
          core_enrichment = TRUE,
          label_format = 50,
          orderBy = "NES",
          decreasing = FALSE)
dev.off()


# SNV-----------------------------------------------------------------------------------------
library(SummarizedExperiment)
library(TCGAbiolinks)
library(readr)
#tcga_project <- "TCGA-SARC"
tcga_project <- "SARC"
#Mutation data (hg38)
maf <- GDCquery_Maf(tcga_project, pipelines = "mutect2") #mutect2
table(maf$Variant_Classification)
library(maftools)
#keep all mut
include.silent = c("3'Flank", "3'UTR","5'Flank", "5'UTR", "IGR","Intron",#"Targeted_Region"
                   "RNA","Silent",
                   "Frame_Shift_Del", "Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
                   "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
                   "Splice_Site","Splice_Region","Translation_Start_Site")
maffile = "TCGA.SARC.mutect.somatic.maf"
maf_all = read.maf(maffile,
                   isTCGA = T,
                   useAll = T,
                   vc_nonSyn = include.silent) #maf file must read by maftool can't read.table
maf_all = maf_all@data
nrow(maf_all[!duplicated(maf_all$Tumor_Sample_Barcode),])
write.table(maf_all,"SARC_n237.allmut_simID.maf",sep = "\t",row.names = F,quote = F)

clin1 = read.table("Risktable_tcga.sarc_0.4668533.txt",sep = "\t",header = T,check.names = F)
clin1$sample_id = rownames(clin1)
clin2 = read.table("SARC_n261_210914.txt",sep = "\t",header = T,check.names = F)
clin= merge(clin1,clin2,by="sample_id")
colnames(clin)[1] = "Tumor_Sample_Barcode" 
clin = clin[,c(1,3,11,12,14,25)]
colnames(clin)[2] = "Virtual_Status"
colnames(clin)[6] = "Metastatic_Disease" 
clin$Virtual_Status = gsub("0","Alive",clin$Virtual_Status)
clin$Virtual_Status = gsub("1","Dead",clin$Virtual_Status)
clin[is.na(clin$Metastatic_Disease),"Metastatic_Disease"] = "N/A"
write.table(clin,"clinformaf.txt",row.names = F,sep = "\t",quote = F)
maffile = "SARC_n237.allmut_simID.maf"
#keep all mut
include.silent = c("3'Flank", "3'UTR","5'Flank", "5'UTR", "IGR","Intron",#"Targeted_Region"
                   "RNA","Silent",
                   "Frame_Shift_Del", "Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
                   "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
                   "Splice_Site","Splice_Region","Translation_Start_Site")
maf.sarc = read.maf(maffile,
                    clinicalData = "clinformaf.txt",
                    isTCGA = T,
                    useAll = T,
                    vc_nonSyn = include.silent)
write.mafSummary(maf=maf.sarc, basename="SARC_n237.allmut_simID") #
library(pals)
vc_cols =as.vector(polychrome(25))[6:25]
names(vc_cols) = c("Frame_Shift_Del","3'Flank", "3'UTR", "5'UTR", "Multi_Hit","IGR","Intron",
                   "RNA","Silent",
                   "Missense_Mutation", "Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
                   "Nonsense_Mutation", "Nonstop_Mutation",
                   "Splice_Site","Splice_Region","Translation_Start_Site","5'Flank")
pdf(file="waterfall.pdf",width=8,height=7.5)
oncoplot(maf = maf.sarc, top=20, colors = vc_cols,gene_mar=6,logColBar=F,draw_titv=TRUE,)
dev.off()
pdf(file="summary.pdf",width=6.5,height=5)
plotmafSummary(maf=maf.sarc, color = vc_cols,
               rmOutlier=TRUE, addStat="median", dashboard=TRUE, titvRaw = T)
dev.off()#
maftmb = tmb(maf = maf.sarc, 
             captureSize = 50, 
             logScale = TRUE)   
write.table(maftmb,"tmb.txt",sep = "\t",row.names = F,quote = F)
# 
clin <- read.table("clinformaf.txt", header=T, check.names = F,sep="\t")
clin.high <- subset(clin, EILncSig_Level=="High")$Tumor_Sample_Barcode
clin.low <- subset(clin, EILncSig_Level=="Low")$Tumor_Sample_Barcode
maftest.high <- subsetMaf(maf=maf.sarc, tsb=clin.high, isTCGA=TRUE)
maftest.low <- subsetMaf(maf=maf.sarc, tsb=clin.low, isTCGA=TRUE)
#
fabcolors1 = c("steelblue","salmon")
names(fabcolors1) = c("Alive", "Dead")
fabcolors2 = c("royalblue","red2")
names(fabcolors2) = c("Low", "High")
fabcolors3 = c("forestgreen","brown","grey60")
names(fabcolors3) = c("No", "Yes","N/A")
fabcolors = list(Virtual_Status = fabcolors1, EILncSig_Level = fabcolors2, Metastatic_Disease = fabcolors3)
pdf(file="waterfall_high.pdf",width=5.5,height=4.5)
oncoplot(maf=maftest.high, top=10, draw_titv=F, gene_mar = 5,barcode_mar = 4,colors = vc_cols,
         clinicalFeatures = c("Virtual_Status","Metastatic_Disease","EILncSig_Level"),
         annotationColor = fabcolors,
         annotationFontSize = 1.7,
         legendFontSize = 1.7)
dev.off()
pdf(file="waterfall_low.pdf",width=5.5,height=4.5)
oncoplot(maf=maftest.low, top=10, draw_titv=F, gene_mar = 5,barcode_mar = 4,colors = vc_cols,
         clinicalFeatures = c("Virtual_Status","Metastatic_Disease","EILncSig_Level"),
         annotationColor = fabcolors,
         annotationFontSize = 1.7,
         legendFontSize = 1.7)
dev.off()
pdf(file="waterfall_low.legend.pdf",width=10,height=10)
oncoplot(maf=maftest.low, top=10, draw_titv=F, gene_mar = 5,barcode_mar = 4,colors = vc_cols,
         clinicalFeatures = c("Virtual_Status","Metastatic_Disease","EILncSig_Level"),
         annotationColor = fabcolors,
         annotationFontSize = 1.7,
         legendFontSize = 1.7)
dev.off()
pdf(file="waterfall_high.legend.pdf",width=10,height=10)
oncoplot(maf=maftest.high, top=10, draw_titv=F, gene_mar = 5,barcode_mar = 4,colors = vc_cols,
         clinicalFeatures = c("Virtual_Status","Metastatic_Disease","EILncSig_Level"),
         annotationColor = fabcolors,
         annotationFontSize = 1.7,
         legendFontSize = 1.7)
dev.off()

# compare specific genes  - take MUC4 as an exmaple
pdf("MUC4.pdf",width = 6.75,height = 3.75) # 6.75*3.75
lollipopPlot2(m1=maftest.high, m2=maftest.low, m1_name="EILncSig_Level High", m2_name="EILncSig_Level Low", 
              gene="MUC4",
              AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short",
              showDomainLabel = FALSE, colors = vc_cols,
              domainLabelSize = 1)
dev.off()

#TMB
#plot
library(dplyr) 
library(ggplot2)
library(ggpubr)
library(ggExtra)
#load
clin = read.table("clinformaf.txt",header = T,sep = "\t")
tmb = read.table("tmb.txt",sep = "\t",header = T)
TMBtest = merge(tmb,clin,by.x = "Tumor_Sample_Barcode",by.y = "Tumor_Sample_Barcode")
colnames(TMBtest)
ggplot(TMBtest,aes(x = EILncSig_Level, y = total_perMB_log)) + 
  geom_violin(aes(fill=EILncSig_Level), trim = FALSE) + 
  geom_boxplot(width = 0.2)+
  scale_fill_manual(values = c("salmon","steelblue"))+
  stat_compare_means(label = "p.format",
                     method = "wilcox.test",
                     hide.ns = T,
                     size = 3.5, 
                     #label.y =2.5,
                     label.y.npc = "top")+
  theme_light() + xlab(NULL) + 
  ylab("Tumor Mutation Burden (per MB)")+
  theme(axis.text.x = element_text(size = 10))
#pdf 5.5*5.5
ggplot(TMBtest, aes(x=EILncSig_Score, y=total_perMB_log)) + 
  geom_point(color="black") + 
  geom_smooth(method="lm", se=T)+
  theme_light()+
  stat_cor(data=TMBtest, method = "pearson")+
  labs(x="", 
       y="", 
       title="")

# CNV-----------------------------------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
query <- GDCquery(project = "TCGA-SARC",
                  data.category = "Copy Number Variation",
                  data.type = "Masked Copy Number Segment") #Copy Number Segment #Masked Copy Number Segment
GDCdownload(query, method = "api")
CNV_download <- GDCprepare(query = query,
                           save = T,
                           save.filename = "TCGA-SARC.CNV_download.rdata")
#
rm(list = ls())
load("TCGA-SARC.CNV_download.rdata")
tumorCNV <- data
tumorCNV <- tumorCNV[,2:7]
tumorCNV <- tumorCNV[,c("Sample", "Chromosome",
                        "Start", "End", "Num_Probes", "Segment_Mean")]
head(tumorCNV)
sampletable = as.data.frame(tumorCNV[!duplicated(tumorCNV$Sample),"Sample"])
sampletable$id1 = substr(sampletable$Sample,14,16)  
table(sampletable$id1)

tum.sam <- substr(tumorCNV$Sample,14,16) %in% c("01A","01B")
table(tum.sam)
tumorCNV <- tumorCNV[tum.sam,]
nrow(tumorCNV[!duplicated(tumorCNV$Sample),])
#check duplicated samples
tumorCNV$sample1 = substr(tumorCNV$Sample,1,16)
tumorCNV$sample2 = substr(tumorCNV$sample1,1,12)
tumorCNV$sample3 = substr(tumorCNV$sample1,14,16)
table(tumorCNV[!duplicated(tumorCNV$sample1),"sample3"])
nrow(tumorCNV[!duplicated(tumorCNV$sample2),])
nrow(tumorCNV[!duplicated(tumorCNV$sample1),])
#check duplicated samples and delete duplicated samples
sampletable$id2 = substr(sampletable$Sample,1,12)
sampletable=sampletable[sampletable$id1 %in% c("01A","01B"),]
sampletable = sampletable[!duplicated(sampletable$id2),]
tumorCNV = tumorCNV[tumorCNV$Sample %in% sampletable$Sample,]
nrow(tumorCNV[!duplicated(tumorCNV$Sample),])
tumorCNV$Sample = tumorCNV$sample2
tumorCNV = tumorCNV[,1:6]
write.table(tumorCNV, file = "segment_sarc.n260.txt", sep = "\t", row.names = F, quote = F)

tumorCNV = read.table("segment_sarc.n260.txt",sep = "\t",header = T)
clin = read.table("Risktable_tcga.sarc_0.4668533.txt",header = T,sep = "\t")
clin.high = rownames(clin[clin$EILncSig_Level == "High",])
clin.low = rownames(clin[clin$EILncSig_Level == "Low",])
tumorCNV.high = tumorCNV[tumorCNV$Sample %in% clin.high,]
tumorCNV.low = tumorCNV[tumorCNV$Sample %in% clin.low,]

nrow(tumorCNV.high[!duplicated(tumorCNV.high$Sample),])
nrow(tumorCNV.low[!duplicated(tumorCNV.low$Sample),])

write.table(tumorCNV.high, file = "segment_sarc.high.n107.txt", sep = "\t", row.names = F, quote = F)
write.table(tumorCNV.low, file = "segment_sarc.low.n151.txt", sep = "\t", row.names = F, quote = F)

# go to GISTIC2.0 - online tools
# download all results (Low/High EILncSig groups) from GISTIC2.0 

#annotate gene name
library(GenomicRanges)
# get a dataframe with c("GeneSymbol","Chr","Start","End") h19/hg38
# load gtf
options(stringsAsFactors = FALSE)
load("Homo_sapiens.GRCh38.104.chr.gtf.Rdata")
gtf_df[is.na(gtf_df$gene_name),"gene_name"] = gtf_df[is.na(gtf_df$gene_name),"gene_id"]
genes = gtf_df
genes = genes[genes$type == "gene",]
genes = genes[,c(1,2,3,12)]
genes = genes[genes[,4]!="" & genes[,1]%in%c(1:22,"X","Y"),]
genes[,1] = sapply(genes[,1],as.character)
genes[genes$seqnames == "X","seqnames"] = "23"
genes[genes$seqnames == "Y","seqnames"] = "24"
genes[,1] = sapply(genes[,1],as.integer)
genes = genes[order(genes[,2]),] #order start
genes = genes[order(genes[,1]),] #order chr
colnames(genes) = c("Chr","Start","End","GeneSymbol")
save(genes,file = "hg38.ens104_chr.gene.position_ordered.rdata")
load("hg38.ens104_chr.gene.position_ordered.rdata")
library(GenomicRanges)
genes_GR = makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)
#load CNV list
sCNV = read.table("./High/High_Risk.all_lesions.conf_90.txt",header = T,check.names = F,sep = "\t")
sCNV = sCNV[,c(1,2,3,6,7)]
sCNV = sCNV[sCNV[,4] <0.05,]
sCNV = sCNV[sCNV[,5] <0.05,]
sCNV = sCNV[!duplicated(sCNV$Descriptor),]
write.table(sCNV,"High_Risk.widepeakrange.txt",sep = "\t",row.names = F)
# use excel to check the txt file and re-save
sCNV = read.table("High_Risk.widepeakrange.txt",header = T,check.names = F,sep = "\t")
sCNV$Chr = gsub("chr","",sCNV$Chr)
sCNV$Chr = as.integer(sCNV$Chr)
sCNV$Start = as.integer(sCNV$Start)
sCNV$End = as.integer(sCNV$End)
sCNV = sCNV[order(sCNV[,4]),] #order start
sCNV = sCNV[order(sCNV[,3]),] #order chr
sCNV_GR = makeGRangesFromDataFrame(sCNV,keep.extra.columns = TRUE)
hits = findOverlaps(genes_GR, sCNV_GR, type="within")
sCNV_ann = cbind(sCNV[subjectHits(hits),],genes[queryHits(hits),])
sCNV_ann = sCNV_ann[!duplicated(sCNV_ann$GeneSymbol),]
sCNV_ann$Aberrant_Region = paste0(sCNV_ann[,3],":",sCNV_ann[,4],"每",sCNV_ann[,5])
sCNV_ann$Gene_Region = paste0(sCNV_ann[,8],":",sCNV_ann[,9],"每",sCNV_ann[,10])
sCNV_ann = sCNV_ann[,c(1,2,6,7,11,12,13)]
sCNV_ann = sCNV_ann[,c("Unique Name","Descriptor","Aberrant_Region","q-values",
                       "Residual q values after removing segments shared with higher peaks",
                       "Gene_Region","GeneSymbol")]
write.table(sCNV_ann,"sCNV_genename_high.txt",quote = F,sep = "\t",row.names = F)
#
load("hg38.ens104_chr.gene.position_ordered.rdata")
library(GenomicRanges)
genes_GR = makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)
#load CNV list
sCNV = read.table("./Low/Low_Risk.all_lesions.conf_90.txt",header = T,check.names = F,sep = "\t")
sCNV = sCNV[,c(1,2,3,6,7)]
sCNV = sCNV[sCNV[,4] <0.05,]
sCNV = sCNV[sCNV[,5] <0.05,]
sCNV = sCNV[!duplicated(sCNV$Descriptor),]
write.table(sCNV,"Low_Risk.widepeakrange.txt",sep = "\t",row.names = F)
# use excel to check the txt file and re-save
sCNV = read.table("Low_Risk.widepeakrange.txt",header = T,check.names = F,sep = "\t")
sCNV$Chr = gsub("chr","",sCNV$Chr)
sCNV$Chr = as.integer(sCNV$Chr)
sCNV$Start = as.integer(sCNV$Start)
sCNV$End = as.integer(sCNV$End)
sCNV = sCNV[order(sCNV[,4]),] #order start
sCNV = sCNV[order(sCNV[,3]),] #order chr
sCNV_GR = makeGRangesFromDataFrame(sCNV,keep.extra.columns = TRUE)
hits = findOverlaps(genes_GR, sCNV_GR, type="within")
sCNV_ann = cbind(sCNV[subjectHits(hits),],genes[queryHits(hits),])
sCNV_ann = sCNV_ann[!duplicated(sCNV_ann$GeneSymbol),]
sCNV_ann$Aberrant_Region = paste0(sCNV_ann[,3],":",sCNV_ann[,4],"每",sCNV_ann[,5])
sCNV_ann$Gene_Region = paste0(sCNV_ann[,8],":",sCNV_ann[,9],"每",sCNV_ann[,10])
sCNV_ann = sCNV_ann[,c(1,2,6,7,11,12,13)]
sCNV_ann = sCNV_ann[,c("Unique Name","Descriptor","Aberrant_Region","q values",
                       "Residual q values after removing segments shared with higher peaks",
                       "Gene_Region","GeneSymbol")]
write.table(sCNV_ann,"sCNV_genename_low.txt",quote = F,sep = "\t",row.names = F)
#

#maftools
options(stringsAsFactors = FALSE)
library(maftools)
test.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
test.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 
test = read.maf(maf = test.maf)
# here maf is used for the gene position
# load h38/19 position file
colnames(genes) = c("Chromosome","Start_Position","End_Position","Hugo_Symbol")
ncol(laml@data)
emptyDF <- t(data.frame(row1 = 1:17))
colnames(emptyDF) = colnames(laml@data)
emptyDF =emptyDF[rep(1,60568),]
rownames(emptyDF) = 1:60568
emptyDF=as.data.frame(emptyDF)
emptyDF$Hugo_Symbol = genes$Hugo_Symbol
emptyDF$Chromosome = genes$Chromosome
emptyDF$Start_Position = genes$Start_Position
emptyDF$End_Position = genes$End_Position
class(emptyDF) = c( "data.table" ,"data.frame")
laml@data = emptyDF
#Done for background
sarc_high.gistic <- readGistic(gisticAllLesionsFile = './High/High_Risk.all_lesions.conf_90.txt',
                               gisticAmpGenesFile = './High/High_Risk.amp_genes.conf_90.txt',
                               gisticDelGenesFile = './High/High_Risk.del_genes.conf_90.txt',
                               gisticScoresFile = './High/High_Risk.scores.gistic',
                               isTCGA = F)
pdf("GISTIC2_HighRisk.pdf",width = 9,height = 3.0)
gisticChromPlot(gistic = sarc_high.gistic, 
                ref.build = 'hg38',#hg38/hg19
                fdrCutOff = 0.1, #0.25"
                markBands = "none", # NULL / "all"
                mutGenes = NULL,
                txtSize = 0.5,
                cytobandTxtSize = 0.7,
                cytobandOffset = 0.07,
                y_lims = c(2,-1)
) 
dev.off()

pdf("GISTIC2_HighRisk_ref.pdf",width = 9,height = 3.0)
gisticChromPlot(gistic = sarc_high.gistic, 
                ref.build = 'hg38',#hg38/hg19
                fdrCutOff = 0.05, #0.25"
                markBands = "all", # NULL / "all"
                mutGenes = NULL,
                txtSize = 0.5,
                cytobandTxtSize = 0.7,
                cytobandOffset = 0.1,
                y_lims = c(2,-1)
) 
dev.off()

sarc_low.gistic <- readGistic(gisticAllLesionsFile = './Low/Low_Risk.all_lesions.conf_90.txt',
                              gisticAmpGenesFile = './Low/Low_Risk.amp_genes.conf_90.txt',
                              gisticDelGenesFile = './Low/Low_Risk.del_genes.conf_90.txt',
                              gisticScoresFile = './Low/Low_Risk.scores.gistic',
                              isTCGA = T)

pdf("GISTIC2_LowRisk.pdf",width = 9,height = 3)
gisticChromPlot(gistic = sarc_low.gistic, 
                ref.build = 'hg38',#hg38/hg19
                fdrCutOff = 0.05, #0.1
                markBands = "none", # NULL / "all"
                txtSize = 0.5,
                cytobandTxtSize = 0.7,
                cytobandOffset = 0.07,
                y_lims = c(2,-1))
dev.off()

pdf("GISTIC2_LowRisk_ref.pdf",width = 9,height = 3)
gisticChromPlot(gistic = sarc_low.gistic, 
                ref.build = 'hg38',#hg38/hg19
                fdrCutOff = 0.05, #0.1
                markBands = "all", # NULL / "all"
                txtSize = 0.5,
                cytobandTxtSize = 0.7,
                cytobandOffset = 0.1,
                y_lims = c(2,-1)) 
dev.off()

# use the gain loss cytoband count from the GIStic2 results (thresholded.by_genes.txt but no duplicated rows)
cnv.h = read.table("CNV gain loss cytoband_High.txt",header = T,sep = "\t",check.names = F)
cnv.l = read.table("CNV gain loss cytoband_Low.txt",header = T,sep = "\t",check.names = F)
rownames(cnv.h) = cnv.h$Cytoband
rownames(cnv.l) = cnv.l$Cytoband
cnv.h=cnv.h[,-1]
cnv.l=cnv.l[,-1]
cnv.h = as.data.frame(t(cnv.h))
cnv.l = as.data.frame(t(cnv.l))
cnv.h.Gain = apply(cnv.h,1,function(x) sum(x[x > 0]))
cnv.h.Loss = apply(cnv.h,1,function(x) sum(x[x < 0]))
cnv.l.Gain = apply(cnv.l,1,function(x) sum(x[x > 0]))
cnv.l.Loss = apply(cnv.l,1,function(x) sum(x[x < 0]))
cnv_h = cbind(ID = row.names(cnv.h),
              Gain = cnv.h.Gain,
              Loss = cnv.h.Loss,
              EILncSig_Level = "High")
cnv_l = cbind(ID = row.names(cnv.l),
              Gain = cnv.l.Gain,
              Loss = cnv.l.Loss,
              EILncSig_Level = "Low")
CNV_events = rbind(cnv_h,cnv_l)
#
risk = read.table("Risktable_tcga.sarc_0.4668533.txt",header = T,sep = "\t")
risk$ID = rownames(risk)
CNV_events = merge(CNV_events,risk,by="ID")
CNV_events = CNV_events[,c(1,2,3,4,14)]
colnames(CNV_events)[4] = "EILncSig_Level"
CNV_events$Gain = as.integer(CNV_events$Gain)
CNV_events$Loss = 0 - as.integer(CNV_events$Loss)
CNV_events$Total.Events = CNV_events$Gain + CNV_events$Loss 
write.table(CNV_events,"CNV_events.txt",sep = "\t",quote = F,row.names = F)
CNV_events = read.table("CNV_events.txt",header = T,sep = "\t")
#
library(ggplot2)
library(ggpubr)
rt <- NULL
for (i in c("Gain","Loss")) {
  idx.sub <- which(colnames(CNV_events) == i)
  sub <- data.frame(
    ID = CNV_events$ID,
    CNV = i,
    EILncSig_Level = CNV_events$EILncSig_Level,
    EILncSig_Score = CNV_events$EILncSig_Score,
    Total.Events = CNV_events$Total.Events,
    CNV_Events = CNV_events[,idx.sub]
  )
  rt <- rbind(rt, sub)
}
ggplot(data = rt, 
       aes(x = reorder(ID, -Total.Events), 
           y = CNV_Events)) +  
  geom_bar(stat = "identity", aes(fill = CNV)) +
  theme_light() +
  ylab("CNV Events") + 
  xlab(NULL)+
  ggtitle("") + theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c("salmon","steelblue")) 
colnames(CNV_events)
ggplot(CNV_events, 
       aes(x = EILncSig_Level,
           y = Gain)) +
  geom_violin(aes(fill=EILncSig_Level), trim = FALSE) + 
  geom_boxplot(width = 0.1)+
  scale_fill_manual(values = c("salmon","steelblue"))+
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y = 550 
  )+
  theme_light() + xlab(NULL) +  ylab("CNV Events")+
  ggtitle("CNV amplication") +
  theme(axis.text.x = element_blank()) + coord_cartesian(ylim=c(-120, 570))
ggplot(CNV_events, 
       aes(x = EILncSig_Level,
           y = Loss)) +
  geom_violin(aes(fill=EILncSig_Level), trim = FALSE) + 
  geom_boxplot(width = 0.1)+
  scale_fill_manual(values = c("salmon","steelblue"))+
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y = 550 
  )+
  theme_light() + xlab(NULL) +  ylab("CNV Events")+
  ggtitle("CNV deletion") +
  theme(axis.text.x = element_blank()) + coord_cartesian(ylim=c(-120, 570))
ggplot(CNV_events, aes(x=EILncSig_Score, y=Gain)) + 
  geom_point(color="black",size = 1.25) + 
  geom_smooth(method="lm", se=T)+
  #geom_smooth()+
  theme_light()+
  stat_cor(data=CNV_events, method = "spearman",
           digits = 3,
           label.y= 450)+
  labs(x="EILncSig Score", 
       y="CNV Events", 
       title="CNV amplication")