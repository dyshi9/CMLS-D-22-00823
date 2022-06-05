load("pan.sarcoma.array.symbol.n1085.sva.rdata")
load("sample.batch.n1085.rdata")
emt.c = read.table("EMT_cluster.txt",header = T,sep = "\t")
colnames(emt.c)[1] = "ID"
sample.batch.n1085 = merge(sample.batch.n1085,emt.c,by="ID")
# EMT cluster to factor
table(sample.batch.n1085$cluster)
sample.batch.n1085[sample.batch.n1085$cluster == "1","cluster"] = 1
sample.batch.n1085[sample.batch.n1085$cluster == "2","cluster"] = 0
table(sample.batch.n1085$cluster)
class(sample.batch.n1085$cluster)
sample.batch.n1085$cluster = as.factor(sample.batch.n1085$cluster) 
#delete useless col
sample.batch.n1085 = sample.batch.n1085[,c(1,3,4,7,8,14)]
#change name
table(sample.batch.n1085$Sarcoma_Type2)
sample.batch.n1085[sample.batch.n1085$Sarcoma_Type2 == "Sarcoma (not specified)",
                   "Sarcoma_Type2"] = "Sarcoma_ns"
#fix age
class(sample.batch.n1085$Age)
sample.batch.n1085$Age= as.numeric(sample.batch.n1085$Age)
#fix gender
table(sample.batch.n1085$Gender)
sample.batch.n1085$Gender = as.factor(sample.batch.n1085$Gender)
sample.batch.n1085$Gender = gsub("female","0",sample.batch.n1085$Gender)
sample.batch.n1085$Gender = gsub("Female","0",sample.batch.n1085$Gender)
sample.batch.n1085$Gender = gsub("male","1",sample.batch.n1085$Gender)
sample.batch.n1085$Gender = gsub("Male","1",sample.batch.n1085$Gender)
#fix metastasis
table(sample.batch.n1085$Metastasis)
sample.batch.n1085$Metastasis = gsub("No","0",sample.batch.n1085$Metastasis)
sample.batch.n1085$Metastasis = gsub("Yes","1",sample.batch.n1085$Metastasis)
#to a new table
sample.trait = sample.batch.n1085
colnames(sample.trait)[4] = "Sarcoma_Type"
colnames(sample.trait)[6] = "EMT_Cluster"
sample.trait$Metastasis = as.factor(sample.trait$Metastasis)
sample.trait$Sarcoma_Type = as.factor(sample.trait$Sarcoma_Type)
# check class
class(sample.trait$Age)
class(sample.trait$Gender)
class(sample.trait$Sarcoma_Type)
class(sample.trait$Metastasis)
class(sample.trait$EMT_Cluster)

save(sample.trait,file = "sample_trait.rdata")
#------------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
#
expr = pan.sarcoma.array.symbol.n1085.sva
expr = as.matrix(expr)
save(expr,file="expr.rdata")

load("sample_trait.rdata")
load("expr.rdata")

datExpr0=as.data.frame(t(expr)) #

library(WGCNA)
#评估矩阵信息是否合格
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#
sampleTree = hclust(dist(datExpr0), method = "average")
#
sizeGrWindow(12,9)
#
par(cex = 0.45);
par(mar = c(1,4,2,1)) #0,4,2,0
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub="", xlab="",cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 185, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 185, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#pdf 7*4

#=========================================================
#load
traitData = sample.trait

dim(traitData)
names(traitData)

# 
tumorSamples = rownames(datExpr);
#
traitRows = match(tumorSamples, traitData$ID);
#
traitRows <- na.omit(traitRows)
#
datTraits = traitData[traitRows, -1];
#
rownames(datTraits) = traitData[traitRows, 1];
#
collectGarbage()


matchnum <- match(row.names(datTraits),row.names(datExpr))
datExpr <- datExpr[matchnum,]
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
#
# Convert traits to a color representation: white means low, red means high, grey means missing entry

### delete Sarcoma_Type & metastasis
datTraits = datTraits[,c(-3,-4)]
###
#factor to numeric
datTraits$Age = as.numeric(datTraits$Age)
datTraits$Gender = as.numeric(datTraits$Gender)
#datTraits$Metastasis = as.numeric(datTraits$Metastasis)
datTraits$EMT_Cluster = as.numeric(datTraits$EMT_Cluster)
#factor to numeric

datTraitsColor <- numbers2colors(datTraits, signed = FALSE);
#
# Plot the sample dendrogram and the colors underneath.
sizeGrWindow(12,9)
plotDendroAndColors(sampleTree2, datTraitsColor,
                    groupLabels = names(datTraits),
                    colorHeight = 0.2,
                    colorHeightBase = 0.2,
                    colorHeightMax = 0.4,
                    rowWidths = NULL,
                    dendroLabels = NULL,
                    addGuide = FALSE, guideAll = FALSE,
                    guideCount = 50, guideHang = 0.2,
                    addTextGuide = FALSE,
                    cex.colorLabels = 0.8,
                    cex.dendroLabels = 0.7, 
                    cex.rowText = 0.8,
                    marAll = c(1, 8.5, 3, 1), saveMar = TRUE, #c(1, 5, 3, 1)
                    main = "Sample dendrogram and trait heatmap")
save(datExpr, datTraits,file = "G-01-dataInput.RData")



# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
rm(list = ls())
lnames = load(file = "G-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
#
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",
     type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,cex=cex1,col="red")


#
# When you have relatively few genes (<5000) use the following code
# here we define the adjacency matrix using soft thresholding with beta=
ADJ1=abs(cor(datExpr,use="p"))^7
k=as.vector(apply(ADJ1,2,sum, na.rm=T))

# When you have a lot of genes use the following code
k=softConnectivity(datE=datExpr,power=7)

# Plot a histogram of k and a scale free topology plot
# check the quality of the scale free topology
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")


#=========================================================
softPower = 7;
## adjacency = adjacency(datExpr, power = softPower)


dissTOM = 1-TOMsimilarity(adjacency(datExpr, power = softPower))

save(dissTOM,file = "dissTOM.rdata")

geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (denodrogram)
sizeGrWindow(12,12)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


minModuleSize = 30; 
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro =FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# 
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# 
sizeGrWindow(8,12)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#7.0*3.5 pdf


MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

MEDiss = 1-cor(MEs);
#
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
#
MEDissThres = 0.30
# 
abline(h=MEDissThres, col = "red")
# 
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
#===================================================================


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleColors = mergedColors
# 
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "G-02-networkConstruction-StepByStep.RData")

library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "G-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "G-02-networkConstruction-StepByStep.RData");
rm(adjacency,dissTOM)

nGenes = ncol(datExpr);
#
nSamples = nrow(datExpr);
# 
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
#
MEs = orderMEs(MEs0)
#
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

row.names(MEs0)=row.names(datExpr)
write.table(MEs0,file="MEs0.txt",sep='\t',quote=F,row.names=T)

sizeGrWindow(10,8)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 9.5, 3, 3)); 
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               #naColor = "grey",
               cex.lab.x = 0.75,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

EMT_Cluster = as.data.frame(datTraits$EMT_Cluster);
names(EMT_Cluster) = "EMT_Cluster"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
#
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue =  as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
#
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
#
geneTraitSignificance = as.data.frame(cor(datExpr, EMT_Cluster, use ="p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples));
names(geneTraitSignificance) = paste("GS.", names(EMT_Cluster), sep="");
names(GSPvalue) = paste("p.GS.", names(EMT_Cluster), sep="");
#
#choose specific color/module
module = "royalblue"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for EMT_Cluster",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module) #6*4


geneInfo0 = data.frame(geneSymbol = rownames(geneTraitSignificance),moduleColor = moduleColors,geneTraitSignificance,GSPvalue)

# Order modules by their significance for ER
modOrder = order(-abs(cor(MEs, EMT_Cluster, use = "p")));

# Add module membership information in the chosen order #在上面的表格中加入MM的信息

for (mod in 1:ncol(geneModuleMembership))  
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,modOrder[mod]],MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]],sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.EMT_Cluster));
geneInfo = geneInfo0[geneOrder, ]
#
write.csv(geneInfo, file = "geneInfo.csv")

#
setwd();
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "G-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "G-02-networkConstruction-StepByStep.RData");
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 7);

#
plotTOM = dissTOM^7;

#
diag(plotTOM) = NA;

#
# Call the plot function
sizeGrWindow(9,9)
#
pdf(file = "Rplot1.pdf", width = 10, height =10)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot",
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon') )
dev.off()


nSelect = 1000
# For reproducibility, we set the random seed
set.seed(520);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;

TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

#
pdf(file = "Rplot1.pdf", width = 8, height =8)

TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, randomly selected genes (n = 1000)",
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon') )

dev.off()
#===================================================================

MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# 
EMT_Cluster = as.data.frame(datTraits$EMT_Cluster);
names(EMT_Cluster) = "EMT_Cluster"
# 
MET = orderMEs(cbind(MEs, EMT_Cluster))

# 
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.8)
#pdf(file = "Rplot.pdf", width = 10, height =10)
plotEigengeneNetworks(MET, "", 
                      marDendro = c(0,4,1,5), marHeatmap = c(7,7.5,1,2), 
                      cex.lab = 0.8,
                      xLabelsAngle = 90)


pdf(file = "11.1.pdf", width = 7.5, height =4)
par(cex =1)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro =c(8.5,4,2,0),plotHeatmaps = FALSE)
dev.off()

# 单独画热图
pdf(file = "11.2.pdf", width = 7.5, height =7.5)
par(cex =0.8)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(9,9,2,2),
                      plotDendrograms = FALSE, plotAdjacency = FALSE,printAdjacency = F, 
                      colorLabels = T,excludeGrey = F, 
                      xLabels = colnames(MET),
                      yLabels = colnames(MET),
                      ySymbols = colnames(MET),
                      xLabelsAngle = 90)
dev.off()
#===================================================================
#mark lncRNA
#===================================================================
load("Homo_sapiens.GRCh38.104.chr.gtf.Rdata")
gtf_df[is.na(gtf_df$gene_name),"gene_name"] = gtf_df[is.na(gtf_df$gene_name),"gene_id"]
gtf_df = gtf_df[gtf_df$type == "gene",]
gtf_df = gtf_df[!duplicated(gtf_df$gene_id),]
table(gtf_df$gene_biotype)
gtf_df_lnc = gtf_df[gtf_df$gene_biotype == "lncRNA",]

geneinfo = read.csv("geneInfo.csv",header = T,row.names = 1)
geneinfo.lnc = geneinfo[geneinfo$geneSymbol %in% gtf_df_lnc$gene_name,]
write.csv(geneinfo.lnc,"geneinfo.lnc.csv")
