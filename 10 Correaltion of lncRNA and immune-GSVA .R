#GSVA
load("pan.sarcoma.array.symbol.n1085.sva.rdata")
library(GSVA) 
#prepare a txt file with immune gene Category, gene id, gene symbol from the Immport database and save as .gmt file. 
original_gmt_GSVA <- readLines('immport_Symbol.ensembl104.gmt') 
strsplit_no_name <- function(gmt.list_layer){ 
  as.character(unlist(strsplit(gmt.list_layer, split = '\t', fixed = T)))}
database_list_GSVA <- lapply(original_gmt_GSVA, strsplit_no_name) 
for (layers in 1:length(database_list_GSVA)) { 
  names(database_list_GSVA)[layers] <- database_list_GSVA[layers][[1]][1] 
  database_list_GSVA[layers][[1]] <- database_list_GSVA[layers][[1]][-1]  # [-2] 去除第二列
}
rm(layers,original_gmt_GSVA,strsplit_no_name)

#load matrix row=gene col=sample
expr = pan.sarcoma.array.symbol.n1085.sva
es <- gsva(as.matrix(expr), 
           database_list_GSVA, 
           mx.diff=TRUE, verbose=TRUE, 
           method='gsva', kcdf='Gaussian', parallel.sz=8)
write.table(es,"GSVA_immport.txt",row.names = T,sep="\t",quote = F)

load("Homo_sapiens.GRCh38.104.chr.gtf.Rdata")
gtf_df[is.na(gtf_df$gene_name),"gene_name"] = gtf_df[is.na(gtf_df$gene_name),"gene_id"]
gtf_df = gtf_df[gtf_df$type == "gene",]
gtf_df = gtf_df[!duplicated(gtf_df$gene_id),]
table(gtf_df$gene_biotype)
gtf_df_lnc = gtf_df[gtf_df$gene_biotype == "lncRNA",]
load("pan.sarcoma.array.symbol.n1085.sva.rdata")
expr=pan.sarcoma.array.symbol.n1085.sva
expr.lnc = expr[rownames(expr) %in% gtf_df_lnc$gene_name,]
GSVA = read.table("GSVA_immport.txt",header = T,sep = "\t",check.names = F)

#cor test
#row-sample; column-variable
library(psych)
expr.lnc = t(expr.lnc)
GSVA = t(GSVA)
GSVA = as.data.frame(GSVA)
expr.lnc = as.data.frame(expr.lnc)
#calculate
Totalsigout=data.frame()
for(i in colnames(GSVA)){
  allOutTab=data.frame()
  GSVA.i = as.data.frame(GSVA[,i])
  colnames(GSVA.i) = i
  corr_test.res = corr.test(x = GSVA.i, 
                            y = expr.lnc,
                            use = "pairwise", #complete
                            method="pearson",
                            adjust="BH", #none
                            alpha=0.05,
                            ci=FALSE,
                            minlength=5)
  cor.r = corr_test.res$r
  cor.p = corr_test.res$p
  cor.adjp = corr_test.res$p.adj
  allOutTab = t(rbind(cor.r,cor.p,cor.adjp))
  colnames(allOutTab) = c("cor","p","fdr")
  allOutTab = as.data.frame(allOutTab)
  Totalsigout =rbind(Totalsigout,
                     cbind(Gene = row.names(allOutTab),
                           Immport_term = i,
                           Cor = allOutTab[,1],
                           Pval = allOutTab[,2],
                           FDR = allOutTab[,3]))
  } #end cycle
Totalsigout$FDR = as.numeric(Totalsigout$FDR)
Totalsigout$Pval = as.numeric(Totalsigout$Pval)
Totalsigout$Cor = as.numeric(Totalsigout$Cor)
Totalsigout.table = Totalsigout[abs(Totalsigout$Cor)>= 0.3 & Totalsigout$FDR < 0.05,]
#end
table(duplicated(Totalsigout.table$Gene))
table(Totalsigout.table$Immport_term)
write.table(Totalsigout.table,"cor_lnc&immport_0.3.txt",row.names = F,sep = "\t",quote = F)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
load("Homo_sapiens.GRCh38.104.chr.gtf.Rdata")
gtf_df[is.na(gtf_df$gene_name),"gene_name"] = gtf_df[is.na(gtf_df$gene_name),"gene_id"]
gtf_df = gtf_df[gtf_df$type == "gene",]
gtf_df = gtf_df[!duplicated(gtf_df$gene_id),]
table(gtf_df$gene_biotype)
gtf_df_lnc = gtf_df[gtf_df$gene_biotype == "lncRNA",]
load("pan.sarcoma.array.symbol.n1085.sva.rdata")
expr=pan.sarcoma.array.symbol.n1085.sva
expr.lnc = expr[rownames(expr) %in% gtf_df_lnc$gene_name,]
#prepare a file with immune gene symbol curated by the immport database, and save as rdata
load("immport_Symbol.ensembl104.rdata")
expr.immune = expr[rownames(expr) %in% immport_Symbol.ensembl104$Symbol_updated,]
#cor test
#row-sample; column-variable
library(psych)
expr.lnc = t(expr.lnc)
expr.immune = t(expr.immune)
expr.immune = as.data.frame(expr.immune)
expr.lnc = as.data.frame(expr.lnc)

#calculate
Totalsigout=data.frame()
for(i in colnames(expr.immune)){
  allOutTab=data.frame()
  expr.immune.i = as.data.frame(expr.immune[,i])
  colnames(expr.immune.i) = i
  corr_test.res = corr.test(x = expr.immune.i, 
                            y = expr.lnc,
                            use = "pairwise", #complete
                            method="pearson",
                            adjust="BH", #none
                            alpha=0.05,
                            ci=FALSE,
                            minlength=5)
  cor.r = corr_test.res$r
  cor.p = corr_test.res$p
  cor.adjp = corr_test.res$p.adj
  allOutTab = t(rbind(cor.r,cor.p,cor.adjp))
  colnames(allOutTab) = c("cor","p","fdr")
  allOutTab = as.data.frame(allOutTab)
  Totalsigout =rbind(Totalsigout,
                     cbind(Gene = row.names(allOutTab),
                           Immport_term = i,
                           Cor = allOutTab[,1],
                           Pval = allOutTab[,2],
                           FDR = allOutTab[,3]))
} #end cycle
Totalsigout$FDR = as.numeric(Totalsigout$FDR)
Totalsigout$Pval = as.numeric(Totalsigout$Pval)
Totalsigout$Cor = as.numeric(Totalsigout$Cor)
Totalsigout.table = Totalsigout[abs(Totalsigout$Cor)>= 0.3 & Totalsigout$FDR < 0.05,]
#end

table(duplicated(Totalsigout.table$Gene))
write.table(Totalsigout.table,"cor_lnc&immune.gene_0.3.txt",row.names = F,sep = "\t",quote = F)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
load("Homo_sapiens.GRCh38.104.chr.gtf.Rdata")
gtf_df[is.na(gtf_df$gene_name),"gene_name"] = gtf_df[is.na(gtf_df$gene_name),"gene_id"]
gtf_df = gtf_df[gtf_df$type == "gene",]
gtf_df = gtf_df[!duplicated(gtf_df$gene_id),]
table(gtf_df$gene_biotype)
gtf_df_lnc = gtf_df[gtf_df$gene_biotype == "lncRNA",]
load("pan.sarcoma.array.symbol.n1085.sva.rdata")
expr=pan.sarcoma.array.symbol.n1085.sva
expr.lnc = expr[rownames(expr) %in% gtf_df_lnc$gene_name,]
tme.cell =  read.table("CIBERSORTx_pan1085_Adjusted.txt",check.names = F,header = T,sep = "\t")
row.names(tme.cell) = tme.cell[,1]
tme.cell = tme.cell[,2:23]
tme.cell =as.data.frame(t(tme.cell))
#cor test
#row-sample; column-variable
library(psych)
expr.lnc = t(expr.lnc)
expr.lnc = as.data.frame(expr.lnc)
tme.cell = t(tme.cell)
tme.cell = as.data.frame(tme.cell)
#calculate
Totalsigout=data.frame()
for(i in colnames(tme.cell)){
  allOutTab=data.frame()
  tme.cell.i = as.data.frame(tme.cell[,i])
  colnames(tme.cell.i) = i
  corr_test.res = corr.test(x = tme.cell.i, 
                            y = expr.lnc,
                            use = "pairwise", 
                            method="pearson", 
                            adjust="BH",
                            alpha=0.05,
                            ci=FALSE,
                            minlength=5)
  cor.r = corr_test.res$r
  cor.p = corr_test.res$p
  cor.adjp = corr_test.res$p.adj
  allOutTab = t(rbind(cor.r,cor.p,cor.adjp))
  colnames(allOutTab) = c("cor","p","fdr")
  allOutTab = as.data.frame(allOutTab)
  Totalsigout =rbind(Totalsigout,
                     cbind(Gene = row.names(allOutTab),
                           Immport_term = i,
                           Cor = allOutTab[,1],
                           Pval = allOutTab[,2],
                           FDR = allOutTab[,3]))
} #end cycle
Totalsigout$FDR = as.numeric(Totalsigout$FDR)
Totalsigout$Pval = as.numeric(Totalsigout$Pval)
Totalsigout$Cor = as.numeric(Totalsigout$Cor)
Totalsigout.table = Totalsigout
#end
table(duplicated(Totalsigout.table$Gene))
write.table(Totalsigout.table,"cor_lnc&tme.cell_origninal_p",row.names = F,sep = "\t",quote = F)
