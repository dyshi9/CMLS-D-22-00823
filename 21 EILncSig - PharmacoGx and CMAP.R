library(PharmacoGx)
#download-------------------------------------------------------------
availablePSets(canonical = F)
PSet.CMAP_2016 <- downloadPSet("CMAP_2016")
PSet.CMAP.sig_2016 <- downloadPertSig("CMAP_2016")
save(PSet.CMAP_2016,file = "CMAP_2016.rdata")
save(PSet.CMAP.sig_2016,file = "CMAP_2016.sig.rdata")

load("Homo_sapiens.GRCh38.104.chr.gtf.Rdata")
gtf_df[is.na(gtf_df$gene_name),"gene_name"] = gtf_df[is.na(gtf_df$gene_name),"gene_id"]
gtf_df = gtf_df[,c("gene_id","gene_name")]
gtf_df = gtf_df[!duplicated(gtf_df$gene_id),]
DEG = read.table("DESeq2.diff_sarc_EILncSig.txt",header = T,sep = "\t")
DEG = DEG[DEG$regulation != "No Significance",]
DEG = merge(DEG,gtf_df,by.x="ID",by.y="gene_name")
DEG$abs.stat = abs(DEG$stat)
DEG = DEG[order(DEG$abs.stat,decreasing = T),]
table(DEG$regulation)
load("D:/CRAN/00 code/05 pharmacogx drug/drug.perturbation_cmap2016.rdata")
# 
EILncSig_genes = DEG[,c(11,10)]
rownames(EILncSig_genes) = EILncSig_genes$gene_id
colnames(EILncSig_genes) = c("feature","direction")
EILncSig_genes$direction = gsub("Up","-1",EILncSig_genes$direction)
EILncSig_genes$direction = gsub("Down","1",EILncSig_genes$direction)
EILncSig_genes = EILncSig_genes[1:1000,]

res <- apply(drug.perturbation[,,c("tstat", "fdr")],
             2, function(x, EILncSig){
               return(connectivityScore(x=x,
                                        y=EILncSig[,2,drop=FALSE],
                                        method="fgsea", nperm=100))
             }, EILncSig=EILncSig_genes)

rownames(res) <- c("Connectivity", "P Value")
res <- t(res)
res <- res[order(res[,1], decreasing=TRUE),]
res = as.data.frame(res)

res.sig1 = res[res$`P Value` <0.05,]

write.table(res,"cmap.res_DEG.sign_1000.txt",sep = "\t",quote = F)
write.table(res.sig1,"cmap.res_DEG.sign_1000_sig.txt",sep = "\t",quote = F)

library(ggplot2)
cmap = read.table("cmap.res_DEG.sign_1000_sig.txt",header = T,sep = "\t",check.names = F)
cmap$ID = rownames(cmap)
cmap$direction = NA
cmap[cmap$Connectivity >0,"direction"] = "pos"
cmap[cmap$Connectivity <0,"direction"] = "neg"
colnames(cmap)[1] = "Connectivity Score"
colnames(cmap)
cmap1 = cmap[abs(cmap$`Connectivity Score`) >= 0.25,]
cmap2 = cmap[c(1:10,53:62),]
ggplot(cmap2, 
       aes(x = reorder(ID,`Connectivity Score`), 
           y = `Connectivity Score`, 
           size = `P Value`,
           colour = `Connectivity Score`)
) + 
  geom_point() + 
  scale_colour_gradient2(
    low ="red",
    mid = "white",
    midpoint = 0,
    high ="blue")  + 
  theme_bw() + 
  scale_size_continuous(range = c(4.5, 2))+
  xlab("") + 
  ylab("")+
  ggtitle("") + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0,size=13.5, hjust = 0.55,vjust = 0),
        axis.text.y = element_text(size=15))

