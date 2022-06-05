

library(GSVA) 
library(limma) 
library(pheatmap)

#pan-cancer TME pattern gene signature obtained from https://doi.org/10.1016/j.ccell.2021.04.014 
#prepare the TME pattern gene signature as txt file and save it as gmt file. 

original_gmt_GSVA <- readLines('TME.gmt') 
strsplit_no_name <- function(gmt.list_layer){ 
  as.character(unlist(strsplit(gmt.list_layer, split = '\t', fixed = T)))} 
database_list_GSVA <- lapply(original_gmt_GSVA, strsplit_no_name) 

for (layers in 1:length(database_list_GSVA)) { 
  names(database_list_GSVA)[layers] <- database_list_GSVA[layers][[1]][1] 
  database_list_GSVA[layers][[1]] <- database_list_GSVA[layers][[1]][-1]
  }
rm(layers,original_gmt_GSVA,strsplit_no_name)

#--GSE71118--------------------------------------------------------------------------
#load expr matrix
load("pan.sarcoma.array.symbol.n1085.sva.rdata")
#load clin
load("sample.batch.n1085.rdata")
sample.batch = sample.batch.n1085[!is.na(sample.batch.n1085$MFS.Status),]
sample.batch = sample.batch[sample.batch$DataSet_ID=="GSE71118",]
expr = pan.sarcoma.array.symbol.n1085.sva[,sample.batch$ID]
rm(pan.sarcoma.array.symbol.n1085.sva,sample.batch.n1085)
#load matrix row=gene col=sample
expr_norm_for_GSVA = expr
#rm(out,expr,genes,entrezIDs)
es <- gsva(as.matrix(expr_norm_for_GSVA), 
           database_list_GSVA, 
           mx.diff=TRUE, verbose=TRUE, 
           method='gsva', kcdf='Gaussian', parallel.sz=8)
write.table(es,"GSVA_TME.GSE71118.txt",row.names = T,sep="\t")


#--E-TABM-1202--------------------------------------------------------------------------
#load expr matrix
load("pan.sarcoma.array.symbol.n1085.sva.rdata")
#load clin
load("sample.batch.n1085.rdata")
sample.batch = sample.batch.n1085[!is.na(sample.batch.n1085$OS.Status),]
sample.batch = sample.batch[sample.batch$DataSet_ID=="E-TABM-1202",]
expr = pan.sarcoma.array.symbol.n1085.sva[,sample.batch$ID]
rm(pan.sarcoma.array.symbol.n1085.sva,sample.batch.n1085)
#load matrix row=gene col=sample
expr_norm_for_GSVA = expr
#rm(out,expr,genes,entrezIDs)
es <- gsva(as.matrix(expr_norm_for_GSVA), 
           database_list_GSVA, 
           mx.diff=TRUE, verbose=TRUE, 
           method='gsva', kcdf='Gaussian', parallel.sz=8)
write.table(es,"GSVA_TME.E-TABM-1202.txt",row.names = T,sep="\t")


#--TCGA-SARC--------------------------------------------------------------------------
#load expr matrix
load("vst_TCGA.SARC_n259.rdata")
#load matrix row=gene col=sample
expr_norm_for_GSVA = vsd
#rm(out,expr,genes,entrezIDs)
es <- gsva(as.matrix(expr_norm_for_GSVA), 
           database_list_GSVA, 
           mx.diff=TRUE, verbose=TRUE, 
           method='gsva', kcdf='Gaussian', parallel.sz=8)
write.table(es,"GSVA_TME.TCGA-SARC.txt",row.names = T,sep="\t")

#--TARGET-OS--------------------------------------------------------------------------
#load expr matrix
load("vst_TARGET.OS_n98.rdata")
#load matrix row=gene col=sample
expr_norm_for_GSVA = vsd
#rm(out,expr,genes,entrezIDs)
es <- gsva(as.matrix(expr_norm_for_GSVA), 
           database_list_GSVA, 
           mx.diff=TRUE, verbose=TRUE, 
           method='gsva', kcdf='Gaussian', parallel.sz=8)
write.table(es,"GSVA_TME.TARGET-OS.txt",row.names = T,sep="\t")


#######################################################################################
#read GSVA scores
#median-transform
#center and scale by MAD

library(effectsize)
es = read.table("GSVA_TME.E-TABM-1202.txt",header = T,sep = "\t",check.names = F)
# sample in row, variable in column
es_t = as.data.frame(t(es))
es_t_z <- standardize(es_t, robust = TRUE) 
es = as.data.frame(t(es_t_z))
rm(es_t,es_t_z)
write.table(es,"GSVA_TME.E-TABM-1202.standardized.txt",row.names = T,sep="\t")

es = read.table("GSVA_TME.GSE71118.txt",header = T,sep = "\t",check.names = F)
# sample in row, variable in column
es_t = as.data.frame(t(es))
es_t_z <- standardize(es_t, robust = TRUE) 
es = as.data.frame(t(es_t_z))
rm(es_t,es_t_z)
write.table(es,"GSVA_TME.GSE71118.standardized.txt",row.names = T,sep="\t")

es = read.table("GSVA_TME.TARGET-OS.txt",header = T,sep = "\t",check.names = F)
# sample in row, variable in column
es_t = as.data.frame(t(es))
es_t_z <- standardize(es_t, robust = TRUE) 
es = as.data.frame(t(es_t_z))
rm(es_t,es_t_z)
write.table(es,"GSVA_TME.TARGET-OS.standardized.txt",row.names = T,sep="\t")

es = read.table("GSVA_TME.TCGA-SARC.txt",header = T,sep = "\t",check.names = F)
# sample in row, variable in column
es_t = as.data.frame(t(es))
es_t_z <- standardize(es_t, robust = TRUE) 
es = as.data.frame(t(es_t_z))
rm(es_t,es_t_z)
write.table(es,"GSVA_TME.TCGA-SARC.standardized.txt",row.names = T,sep="\t")

#combine
es1 = read.table("GSVA_TME.E-TABM-1202.standardized.txt",header = T,sep = "\t",check.names = F)
es2 = read.table("GSVA_TME.GSE71118.standardized.txt",header = T,sep = "\t",check.names = F)
es3 = read.table("GSVA_TME.TARGET-OS.standardized.txt",header = T,sep = "\t",check.names = F)
es4 = read.table("GSVA_TME.TCGA-SARC.standardized.txt",header = T,sep = "\t",check.names = F)
es = cbind(es1,es2,es3,es4)
rm(es1,es2,es3,es4)
#######################################################################################

library(ConsensusClusterPlus)
workDir="cluster/"

#¾ÛÀà
GSVAscores = as.matrix(es)

maxK=6
results = ConsensusClusterPlus(GSVAscores,
                               maxK=maxK,
                               reps=1000,
                               pItem=0.8,
                               pFeature=1,
                               title=workDir,
                               clusterAlg="km",
                               distance="euclidean", 
                               seed=520,
                               verbose=T,
                               plot= "pdf")

#
clusterNum=4                 #

cluster=results[[clusterNum]][["consensusClass"]]
outTab=cbind(colnames(GSVAscores),cluster)
outTab=as.data.frame(outTab)
colnames(outTab)[1] = "ID"

write.table(outTab,file="TME_4cohorts_4clusters.txt",sep="\t",row.names = F,quote=F)

#####################################################################################################################

#plot heatmap
library(pheatmap)
#combine
es1 = read.table("GSVA_TME.E-TABM-1202.standardized.txt",header = T,sep = "\t",check.names = F)
es2 = read.table("GSVA_TME.GSE71118.standardized.txt",header = T,sep = "\t",check.names = F)
es3 = read.table("GSVA_TME.TARGET-OS.standardized.txt",header = T,sep = "\t",check.names = F)
es4 = read.table("GSVA_TME.TCGA-SARC.standardized.txt",header = T,sep = "\t",check.names = F)
es = cbind(es1,es2,es3,es4)
rm(es1,es2,es3,es4)
#
es.heatmap = es
clusterfile = read.table("TME_4cohorts_4clusters.txt",header = T,check.names = F,sep = "\t")

#
clusterfile = clusterfile[order(clusterfile$cluster),]
table(clusterfile$cluster)
clusterfile[clusterfile$cluster == 1,2] = "Fibrotic"
clusterfile[clusterfile$cluster == 2,2] = "Immune-Enriched Fibrotic"
clusterfile[clusterfile$cluster == 3,2] = "Depleted"
clusterfile[clusterfile$cluster == 4,2] = "Immune-Enriched Non-Fibrotic"

write.table(clusterfile,file="TME_cluster4_annotated.txt",sep="\t",row.names = F,quote=F)


#make annotation - sample datasets
load("sample.batch.n1085.rdata")
sample.batch.n1085 = sample.batch.n1085[,c(1,2)]
rownames(sample.batch.n1085) = sample.batch.n1085$ID

sample.ann = as.data.frame(colnames(es))
colnames(sample.ann)[1] = "ID"
sample.ann$Datasets = NA
for (i in sample.ann$ID)
{
  if (i %in% sample.batch.n1085$ID)
  {
    sample.ann[sample.ann$ID == i,"Datasets"] = sample.batch.n1085[sample.batch.n1085$ID == i,2]
}}
sample.ann[substr(sample.ann$ID,1,4) == "TCGA","Datasets"] = "TCGA-SARC"
sample.ann[substr(sample.ann$ID,1,6) == "TARGET","Datasets"] = "TARGET-OS"
table(sample.ann$Datasets)

#make annotation - sample EILncSig
EILncSig1 = read.table("EMTILncSig/Risktable_training_1.245959.txt",
                       header = T,sep = "\t")
EILncSig2 = read.table("EMTILncSig validation/validation E-TABM-1202/Risktable_E.TABM.1202_1.221731.txt",
                       header = T,sep = "\t")
EILncSig3 = read.table("EMTILncSig validation/validation tcga sarc/Risktable_tcga.sarc_0.4668533.txt",
                       header = T,sep = "\t")
EILncSig4 = read.table("EMTILncSig validation/validation target os/Risktable_target.os_0.4010929.txt",
                       header = T,sep = "\t")
rownames(EILncSig1) = EILncSig1$id
EILncSig1 = EILncSig1[,c("EILncSig_Score","EILncSig_Level")]
EILncSig2 = EILncSig2[,c("EILncSig_Score","EILncSig_Level")]
EILncSig3 = EILncSig3[,c("EILncSig_Score","EILncSig_Level")]
EILncSig4 = EILncSig4[,c("EILncSig_Score","EILncSig_Level")]

EILncSig1 = EILncSig1[order(EILncSig1$EILncSig_Score),]
EILncSig2 = EILncSig2[order(EILncSig2$EILncSig_Score),]
EILncSig3 = EILncSig3[order(EILncSig3$EILncSig_Score),]
EILncSig4 = EILncSig4[order(EILncSig4$EILncSig_Score),]

EILncSig = rbind(EILncSig1,EILncSig2,EILncSig3,EILncSig4)
EILncSig$ID = rownames(EILncSig)

sample.ann = merge(sample.ann,EILncSig,by="ID")
rownames(sample.ann) = sample.ann$ID
sample.ann = sample.ann[,-1] 
sample.ann = sample.ann[,-2] # delete EILncSig_Score

sample.ann = sample.ann[order(sample.ann$EILncSig_Level),] #- sort by EILncSig_Level
sample.ann$ID = rownames(sample.ann)

annsample = merge(sample.ann,clusterfile,by.x = "ID",by.y = "ID")

rownames(annsample) = annsample$ID
annsample=annsample[,-1]
colnames(annsample)[3] = "TME Pattern"
colnames(annsample)[1] = "Datasets"
colnames(annsample)[2] = "EILncSig Level"

colnames(annsample)
save(annsample,file = "annsample.rdata")

#
annsig = read.table("TME-feature.txt",header = T,sep = "\t",check.names = F)
rownames_annsig = annsig[,2]
annsig = as.data.frame(annsig[,1])
rownames(annsig) = rownames_annsig
colnames(annsig) = "TME Functional Signatures"
#
rm(clusterfile,EILncSig,EILncSig1,EILncSig2,EILncSig3,EILncSig4,sample.ann,sample.batch.n1085)


#
library(pals)
my_colour = list(
   `EILncSig Level` =c(High = 'red2', Low = 'royalblue'),
   Datasets = c(`E-TABM-1202`= "forestgreen",`GSE71118`= "brown1", 
                `TCGA-SARC`  = "purple", `TARGET-OS` = "darkorange"),
   `TME Pattern` = c(`Immune-Enriched Fibrotic` = "purple", `Depleted` = "darkorange",
                     `Immune-Enriched Non-Fibrotic` = "skyblue",`Fibrotic` = "brown"),
  `TME Functional Signatures`= c(`Angiogenesis` = "#875692", `Anti tumor microenvironment` = "#F38400",
                       `Antigen presentation` = "#A1CAF1",`B cells` = "#BE0032",
                       `CAF` = "#C2B280", `Checkpoint inhibition` = "#848482",
                       `Cytotoxic T and NK cells` = "#008856",`Granulocytes` = "#E68FAC",
                       `MDSC` = "#0067A5", `Treg` = "#604E97",
                       `Tumor features` = "#F6A600",`Tumor promotive immune infiltrate` = "#B3446C")
)

colnames(annsample)
annsample = annsample[order(annsample$Datasets),]
annsample = annsample[order(annsample$`EILncSig Level`),]
annsample = annsample[order(annsample$`TME Pattern`),]

es.heatmap = es.heatmap[rownames(annsig),] #TME pattern name
es.heatmap1 = es.heatmap[,rownames(annsample)] #TME pattern name


pdf("TME.subtype.heatmap.pdf",height = 5.25,width =8.25)
bk = unique(c(seq(-1.75,1.75, length= 50)))
pheatmap(es.heatmap1, 
         breaks = bk,
         annotation_col = annsample,
         annotation_row = annsig,
         show_colnames = F, 
         show_rownames = T, 
         scale ="row",
         main = "Standardized GSVA Enrichment Score of TME Functional Signatures",
         cluster_cols = F, 
         cluster_rows = F, 
         fontsize = 7,
         gaps_col = c(242,449,632), # /50/ 8,31,36,50,60,69,83
         gaps_row = c(7,15,27), #18,50,69
         annotation_colors = my_colour,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
)
dev.off()

##############################################################################################
##
load("annsample.rdata")
colnames(annsample)
library(ggstatsplot)
library(ggpubr)
library(pals)

EILncSig1 = read.table("D:/CRAN/Working_Files/04 panSARC/8 EMTILncSig/Risktable_training_1.245959.txt",
                       header = T,sep = "\t")
EILncSig2 = read.table("D:/CRAN/Working_Files/04 panSARC/9 EMTILncSig validation/validation E-TABM-1202/Risktable_E.TABM.1202_1.221731.txt",
                       header = T,sep = "\t")
EILncSig3 = read.table("D:/CRAN/Working_Files/04 panSARC/9 EMTILncSig validation/validation tcga sarc/Risktable_tcga.sarc_0.4668533.txt",
                       header = T,sep = "\t")
EILncSig4 = read.table("D:/CRAN/Working_Files/04 panSARC/9 EMTILncSig validation/validation target os/Risktable_target.os_0.4010929.txt",
                       header = T,sep = "\t")
rownames(EILncSig1) = EILncSig1$id
EILncSig1 = EILncSig1[,c("EILncSig_Score","EILncSig_Level")]
EILncSig2 = EILncSig2[,c("EILncSig_Score","EILncSig_Level")]
EILncSig3 = EILncSig3[,c("EILncSig_Score","EILncSig_Level")]
EILncSig4 = EILncSig4[,c("EILncSig_Score","EILncSig_Level")]
EILncSig = rbind(EILncSig1,EILncSig2,EILncSig3,EILncSig4)
EILncSig$ID = rownames(EILncSig)

annsample1 = merge(EILncSig,annsample,by = "row.names")
colnames(annsample1)


##############################################################################################

annsample2 = annsample1[annsample1$`TME Pattern` != "N/A",]
table(annsample2$Datasets)
colnames(annsample2)

ggbarstats(annsample2,
           `TME Pattern`,
           `EILncSig Level`, 
           bar.proptest = T, 
           #package = "pals",
           palette = 'Set1',#
           results.subtitle = T,
           type = "nonparametric",
           proportion.test = F
) + labs(x = '',y = '',title="Association of EILncSig level and TME pattern") + theme_light()+
  theme(axis.text.x = element_text(size = 11.5)) +
  scale_fill_manual(values = c("skyblue","purple","brown","darkorange"))+
  coord_flip() #
#10.2*3.5


#################################################################################
library(Rtsne)
library(scatterplot3d)

data <- es.heatmap
datat = t(data)
datat = as.data.frame(datat)
annsampleRtsne = annsample[order(annsample$`TME Pattern`),]
datat = datat[rownames(annsampleRtsne),]
datat$label = annsampleRtsne$`TME Pattern`
datat$label<-as.factor(datat$label)
datatm = datat[,-ncol(datat)]
datatm = as.matrix(datatm)
table(datat$label)
set.seed(47)
tsne<- Rtsne(datatm, dims = 3, perplexity=30,  verbose=TRUE,max_iter = 1000)
scatterplot3d(tsne$Y[,1:3],
              pch=c(rep(20,242),rep(20,207),rep(20,183),rep(20,135)),
              color=c(rep("darkorange",242),rep("brown",207),rep("purple",183),rep("skyblue",135)),
              angle=60,
              main='t-SNE visualization of TME patterns',
              xlab = "t-SNE_dim1",
              ylab = "t-SNE_dim2",
              zlab = "t-SNE_dim3",
              cex.symbols=1.25)
