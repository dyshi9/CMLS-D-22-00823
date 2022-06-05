#load symbol for refname
load("Homo_sapiens.GRCh38.104.chr.gtf.Rdata")
gtf_df[is.na(gtf_df$gene_name),"gene_name"] = gtf_df[is.na(gtf_df$gene_name),"gene_id"]
gtf_df = gtf_df[,c("gene_id","gene_name")]
gtf_df = gtf_df[!duplicated(gtf_df$gene_id),]

genesymbols = gtf_df
save(genesymbols,file="gene_name.rdata")

#load txt that was downloaded from EMTome website
PMID25214461 = read.table("PMID25214461.txt",skip = 2,header = F)
PMID28680090 = read.table("PMID28680090.txt",skip = 2,header = F)
PMID31069256 = read.table("PMID31069256.txt",skip = 2,header = F)
PMID29293502 = read.table("PMID29293502.txt",skip = 2,header = F)
PMID29346386 = read.table("PMID29346386.txt",skip = 2,header = F)

EMTsig = rbind(PMID25214461,PMID28680090,PMID31069256,PMID29293502,PMID29346386)


EMTsig = as.data.frame(EMTsig[!duplicated(EMTsig$V1),])
colnames(EMTsig) = "symbol"
EMTsig$ref = EMTsig$symbol %in% genesymbols$gene_name
EMTsig$true.symbol = EMTsig$symbol

#update offical gene symbols
"C1orf54" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "C1ORF54","true.symbol"] = "C1orf54" 
"C4orf19" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "C4ORF19","true.symbol"] = "C4orf19" 
"C11orf24" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "C11ORF24","true.symbol"] = "C11orf24" 
"C1orf210" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "C1ORF210","true.symbol"] = "C1orf210" 
"C19orf33" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "C19ORF33","true.symbol"] = "C19orf33" 
"C1orf116" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "C1ORF116","true.symbol"] = "C1orf116" 
"C2orf15" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "C2ORF15","true.symbol"] = "C2orf15" 
"C6orf132" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "C6ORF132","true.symbol"] = "C6orf132" 
"C5orf66" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "C5ORF66","true.symbol"] = "C5orf66" 
"LOC100505880" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "LOC100505880","true.symbol"] = "N/A" 
"DENND2B" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "ST5","true.symbol"] = "DENND2B" 
"DYNC2I2" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "WDR34","true.symbol"] = "DYNC2I2" 
"CRACDL" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "KIAA1211L","true.symbol"] = "CRACDL" 
"PIK3IP1" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "MGC17330","true.symbol"] = "PIK3IP1" 
"ASAP1-IT1" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "ASAP1-IT1","true.symbol"] = "N/A" 
"LINC00312" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "LINC00312","true.symbol"] = "N/A" 
"CCL9" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "CCL9","true.symbol"] = "N/A" 
"EXPI" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "EXPI","true.symbol"] = "N/A" 
"FCRLS" %in% genesymbols$gene_name
EMTsig[EMTsig$symbol == "FCRLS","true.symbol"] = "N/A" 

#
EMTsig$true.ref = EMTsig$true.symbol %in% genesymbols$gene_name
save(list = ls(all.names = TRUE),file = "EMTgenes.rdata")
