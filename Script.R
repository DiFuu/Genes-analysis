#Read the data and convert it to a .txt
library(readxl)

file <- read_excel("Genes_all.xlsx", 1)
write.table(file, "Genes_all.txt",sep="\t")

#We read the file
geneList <- read.table("Genes_all.txt", row.names=NULL, quote="\"", comment.char="", header=TRUE)
geneList <- geneList[,2]

head(geneList)

####################
#GENOMIC ANNOTATION#
####################

library(biomaRt)
library(org.At.tair.db)

mart<- useDataset("athaliana_eg_gene", useMart(biomart="plants_mart",host="https://plants.ensembl.org"))

Genes_annotation <- getBM(
  filters= "ensembl_gene_id",
  attributes= c("ensembl_gene_id", "tair_symbol", "description"),
  values= geneList,
  mart= mart
)

dim(Genes_annotation)

head(Genes_annotation)

write.table(Genes_annotation,file="Genes_annotation.txt",sep="\t",row.names=F)

###############
#GO ENRICHMENT#
###############

library(clusterProfiler) # citation needed (GO enrichment, KEGG enrichment and plots)
library(org.At.tair.db) # AGI codes and other symbols for A. thaliana genes
# keytypes(org.At.tair.db) # To see the different keys (symbol, entry, etc.) available in the database
# head(keys(org.At.tair.db, keytype = "TAIR"))

## https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
library(enrichplot) # plots using enrichgo objects (check the link)

library(DOSE)
library(GOSemSim) #citation needed (simplify internally calls this package)

# Warning: if you want to save your files you can use write.xlsx / write.table
# BUT you will loose the enrichgo format (if you load them again in R you'll have dataframes
# to preserve the format use: save(your_enrichresult, file = "your_enrichresult.RData") and load("your_enrichresult.RData") 

# Basic enrichment function:
GOBP <- 
  enrichGO(geneList, # list of genes (AGI vector - genes in TAIR format -)
           OrgDb = org.At.tair.db, # database / organism
           keyType = "TAIR", # key that you are using (symbol, entry, TAIR format, etc.)
           ont = "BP", # BP, CC or MF
           pvalueCutoff = 0.05,
           pAdjustMethod = "BY", # I used Bonferroni-Yuketielly (?) but you can change to whatever method
           minGSSize = 1, # min number of genes per ontology in the result
           readable = FALSE # if TRUE the program changes ATxGxxxx by names (i.e., ERF2), but if the name is not recorded you'll obtain NA
  )

barplot(GOBP)

#Data
dataBP <- data.frame(GOBP$ID, GOBP$Description, GOBP$GeneRatio, GOBP$BgRatio, GOBP$pvalue, GOBP$p.adjust, GOBP$qvalue, GOBP$geneID, GOBP$Count, GOBP@ontology)
dataBP$GOBP.geneID <- gsub("/",", ", dataBP$GOBP.geneID)
colnames(dataBP) <- c("ID","term","GeneRatio","BgRatio","pvalue","adj_pval","qvalue","genes","count", "category")


write.table(dataBP,file="GOBP.txt",sep="\t",row.names=F)


#KEGG
options(clusterProfiler.download.method = "wininet")
KEGG <-
  enrichKEGG(geneList,
             organism = "ath",
             keyType = "kegg",
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             minGSSize = 10,
             maxGSSize = 500,
             qvalueCutoff = 0.2,
             use_internal_data = FALSE)

barplot(KEGG)

KEGGrd <- setReadable(KEGG, 'org.At.tair.db', 'TAIR')
cnetplot(KEGGrd, showCategory = 20, foldChange = NULL, layout = "kk", colorEdge = T, circular = F, node_label = "all", cex_category = 1, cex_gene = 1, node_label_size = NULL, cex_label_category = 1, cex_label_gene = 1)

######################
#GOPLOT VISUALIZATION#
######################

library(GOplot)

datosgenes <- read.table("datosgenes.txt", header=TRUE)
datosGO <- read.table("GOBP_readable.txt", header=TRUE)
#process <- unique(datosGO[,"term"])
process <- read.table("process.txt", header=TRUE)
process <- as.character(process[,1 ])


#Data (datosGO must contain ID, term, adj_pval, genes, count and category, datosgenes must contain: ID y logFC)
circ <- circle_dat(datosGO,datosgenes)

#We select the data from the category of Biological Process
dataGO <- circ[circ$category == "BP",]
genes <- unique(dataGO[,c("genes","logFC")])

chord <- chord_dat(circ, genes, process)

head(chord) #1 is assigned to the process, 0 means not

chord <- chord_dat(data = circ, genes = genes, process = process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 9, lfc.col = c("black", "white"))

