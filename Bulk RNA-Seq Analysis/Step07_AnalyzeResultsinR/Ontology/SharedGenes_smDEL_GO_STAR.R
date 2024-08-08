library("dplyr")
library("tidyr")
library("DOSE")
library("ggplot2")
library("clusterProfiler")
library("org.Hs.eg.db")
library("biomaRt")
library("data.table")
library("enrichplot")
library("ggnewscale")
library("tibble")

#Set working directory
directory <- "../Cotney_Lab/PWS_RNASeq/STAR/smDEL"
setwd(directory)

#Read in files of shared genes
geneUP <- read.delim("VennDiagram_IntersectUpResultsENSEMBL.txt", sep="\t", head=F)
geneDOWN <- read.delim("VennDiagram_IntersectDownResultsENSEMBL.txt", sep="\t", head=F)

#Load dds object
load("PWS_DESeqSTAR.RData")
#Make list of genes for universe (all genes from alignment)
univ <- row.names(dds)

#Break up into up & down regulated gene lists
geneUP <- geneUP$V1
geneDOWN <- geneDOWN$V1
#Create combined list for both up & down regulated lists
geneALL <- append(geneUP, geneDOWN)

#Read in combined results object
ResGeneID_ALL <- read.csv(file = 'CombinedH9CT2_Results_DESeqSTAR.csv', header = TRUE, row.names = "X")
#Calculate average log2FoldChange
ResGeneID_ALL$log2FoldChange <- rowMeans(ResGeneID_ALL[,c(3,48)])

#Generate object that has DESeq2 results + counts for shared genes
geneALLdf <- as.data.frame(geneALL)
names(geneALLdf)[1] <- "Gene"
rownames(geneALLdf) <- geneALLdf$Gene
resALL_data <- merge(as.data.frame(geneALLdf), as.data.frame(ResGeneID_ALL), by="row.names", sort=FALSE)
resALL_data <- resALL_data[-c(1:2)]
names(resALL_data)[1] <- "Gene"
rownames(resALL_data) <- resALL_data$Gene
#Convert ENSEMBL id's to gene symbol 
resALL_data$symbol <- mapIds(org.Hs.eg.db, keys = resALL_data$Gene, column = "SYMBOL", keytype = "ENSEMBL", multivals= "first")
#and drop NA's (no gene symbol in database)
#resALL_data <- drop_na(resALL_data)

#Use dplyr to make foldchange object
foldChange <- resALL_data %>% 
  dplyr::select(symbol, log2FoldChange) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(log2FoldChange=mean(log2FoldChange))
names(foldChange)[1] <- "Gene"
foldChange <- deframe(foldChange)

##---Gene Ontology---##
GOupBP <- enrichGO(gene = geneUP,
                universe = univ,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)
GOdownBP <- enrichGO(gene = geneDOWN,
                universe = univ,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)
GOupCC <- enrichGO(gene = geneUP,
                   universe = univ,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)
GOdownCC <- enrichGO(gene = geneDOWN,
                     universe = univ,
                     keyType = "ENSEMBL",
                     OrgDb = org.Hs.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05,
                     readable = TRUE)
GOupMF <- enrichGO(gene = geneUP,
                   universe = univ,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)
GOdownMF <- enrichGO(gene = geneDOWN,
                     universe = univ,
                     keyType = "ENSEMBL",
                     OrgDb = org.Hs.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05,
                     readable = TRUE)
allGOBP <- enrichGO(gene = geneALL,
                    universe = univ,
                    keyType = "ENSEMBL",
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)
allGOCC <- enrichGO(gene = geneALL,
                    universe = univ,
                    keyType = "ENSEMBL",
                    OrgDb = org.Hs.eg.db,
                    ont = "CC",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)
allGOMF <- enrichGO(gene = geneALL,
                    universe = univ,
                    keyType = "ENSEMBL",
                    OrgDb = org.Hs.eg.db,
                    ont = "MF",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)

#Write results into table
summary_upBP <- data.frame(GOupBP)
summary_downBP <- data.frame(GOdownBP)
summary_upCC <- data.frame(GOupCC)
summary_downCC <- data.frame(GOdownCC)
summary_upMF <- data.frame(GOupMF)
summary_downMF <- data.frame(GOdownMF)
summary_allBP <- data.frame(allGOBP)
summary_allCC <- data.frame(allGOCC)
summary_allMF <- data.frame(allGOMF)

#Save summary tables
write.table(as.data.frame(summary_upBP), sep = "\t", file = "PWSexp_upBP.tsv", quote = FALSE)
#write.table(as.data.frame(summary_downBP), sep = "\t", file = "PWSexp_downBP.tsv", quote = FALSE)
write.table(as.data.frame(summary_upCC), sep = "\t", file = "PWSexp_upCC.tsv", quote = FALSE)
write.table(as.data.frame(summary_downCC), sep = "\t", file = "PWSexp_downCC.tsv", quote = FALSE)
write.table(as.data.frame(summary_upMF), sep = "\t", file = "PWSexp_upMF.tsv", quote = FALSE)
write.table(as.data.frame(summary_downMF), sep = "\t", file = "PWSexp_downMF.tsv", quote = FALSE)
write.table(as.data.frame(summary_allBP), sep = "\t", file = "PWSexp_allBP.tsv", quote = FALSE)
write.table(as.data.frame(summary_allCC), sep = "\t", file = "PWSexp_allCC.tsv", quote = FALSE)
write.table(as.data.frame(summary_allMF), sep = "\t", file = "PWSexp_allMF.tsv", quote = FALSE)

#Simplify combining similar terms
SimpUpBP <- simplify(GOupBP, cutoff = 0.8, by = "p.adjust", select_fun = min)
summary_SimpUpBP <- data.frame(SimpUpBP)
#SimpDownBP <- simplify(GOdownBP, cutoff = 0.8, by = "p.adjust", select_fun = min)
#summary_SimpDownBP <- data.frame(SimpDownBP)
SimpUpCC <- simplify(GOupCC, cutoff = 0.8, by = "p.adjust", select_fun = min)
summary_SimpUpCC <- data.frame(SimpUpCC)
SimpDownCC <- simplify(GOdownCC, cutoff = 0.8, by = "p.adjust", select_fun = min)
summary_SimpDownCC <- data.frame(SimpDownCC)
SimpUpMF <- simplify(GOupMF, cutoff = 0.8, by = "p.adjust", select_fun = min)
summary_SimpUpMF <- data.frame(SimpUpMF)
SimpDownMF <- simplify(GOdownMF, cutoff = 0.8, by = "p.adjust", select_fun = min)
summary_SimpDownMF <- data.frame(SimpDownMF)
SimpAllBP <- simplify(allGOBP, cutoff = 0.7, by = "p.adjust", select_fun = min)
summary_SimpAllBP <- data.frame(SimpAllBP)
SimpAllCC <- simplify(allGOCC, cutoff = 0.7, by = "p.adjust", select_fun = min)
summary_SimpAllCC <- data.frame(SimpAllCC)
SimpAllMF <- simplify(allGOMF, cutoff = 0.7, by = "p.adjust", select_fun = min)
summary_SimpAllMF <- data.frame(SimpAllMF)
#Save summary tables
write.table(as.data.frame(summary_SimpAllBP), sep = "\t", file = "PWSexp_allBP_simplified.tsv", quote = FALSE)
write.table(as.data.frame(summary_SimpAllCC), sep = "\t", file = "PWSexp_allCC_simplified.tsv", quote = FALSE)
write.table(as.data.frame(summary_SimpAllMF), sep = "\t", file = "PWSexp_allMF_simplified.tsv", quote = FALSE)
write.table(as.data.frame(summary_SimpUpBP), sep = "\t", file = "PWSexp_upBP_simplified.tsv", quote = FALSE)
write.table(as.data.frame(summary_SimpUpCC), sep = "\t", file = "PWSexp_upCC_simplified.tsv", quote = FALSE)
write.table(as.data.frame(summary_SimpUpMF), sep = "\t", file = "PWSexp_upMF_simplified.tsv", quote = FALSE)
write.table(as.data.frame(summary_SimpDownCC), sep = "\t", file = "PWSexp_downCC_simplified.tsv", quote = FALSE)
write.table(as.data.frame(summary_SimpDownMF), sep = "\t", file = "PWSexp_downMF_simplified.tsv", quote = FALSE)

##---Dotplots---##
#For BP all
  #Calculate fold change to use on x-axis of dotplot
BPall_FC <- summary_SimpAllBP
BPall_FC$GeneRatio <- sapply(BPall_FC$GeneRatio, function(x) eval(parse(text=x)))
BPall_FC$BgRatio <- sapply(BPall_FC$BgRatio, function(x) eval(parse(text=x)))
BPall_FC$FoldChange <- BPall_FC$GeneRatio/BPall_FC$BgRatio
BPall_FC$logFoldChange <- log2(BPall_FC$FoldChange)
BPall_FC$neglogPadj <- -log(BPall_FC$p.adjust)
  #Order by p.adj
ordBPall_FC <- order(BPall_FC$p.adjust)
sortBPall_FC <- BPall_FC[ordBPall_FC, ]
  #Only use top 25 terms
sortBPall_FC <- sortBPall_FC[c(0:25),]
  #Order by logFoldChange
sortBPall_FC <- sortBPall_FC[order(sortBPall_FC$logFoldChange),]
sortBPall_FC$Description <- factor(sortBPall_FC$Description, levels = sortBPall_FC$Description)
  #Make dotplots
dotAllBPgg <- ggplot(data = sortBPall_FC, aes(x = logFoldChange, y = Description, 
                                            color = `neglogPadj`, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  labs(title="Biological Processes in delSNORD116 (Shared Genes)", y="", x="logFoldChange")
dotAllBPgg
#For CC all
  #Calculate fold change to use on x-axis of dotplot
CCall_FC <- summary_SimpAllCC
CCall_FC$GeneRatio <- sapply(CCall_FC$GeneRatio, function(x) eval(parse(text=x)))
CCall_FC$BgRatio <- sapply(CCall_FC$BgRatio, function(x) eval(parse(text=x)))
CCall_FC$FoldChange <- CCall_FC$GeneRatio/CCall_FC$BgRatio
CCall_FC$logFoldChange <- log2(CCall_FC$FoldChange)
CCall_FC$neglogPadj <- -log(CCall_FC$p.adjust)
  #Order by p.adj
ordCCall_FC <- order(CCall_FC$p.adjust)
sortCCall_FC <- CCall_FC[ordCCall_FC, ]
  #Order by logFoldChange
sortCCall_FC <- sortCCall_FC[order(sortCCall_FC$logFoldChange),]
sortCCall_FC$Description <- factor(sortCCall_FC$Description, levels = sortCCall_FC$Description)
  #Make dotplots
dotAllCCgg <- ggplot(data = sortCCall_FC, aes(x = logFoldChange, y = Description, 
                                              color = `neglogPadj`, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  labs(title="Cellular Components in delSNORD116 (Shared Genes)", y="", x="logFoldChange")
dotAllCCgg
#For MF all
  #Calculate fold change to use on x-axis of dotplot
MFall_FC <- summary_SimpAllMF
MFall_FC$GeneRatio <- sapply(MFall_FC$GeneRatio, function(x) eval(parse(text=x)))
MFall_FC$BgRatio <- sapply(MFall_FC$BgRatio, function(x) eval(parse(text=x)))
MFall_FC$FoldChange <- MFall_FC$GeneRatio/MFall_FC$BgRatio
MFall_FC$logFoldChange <- log2(MFall_FC$FoldChange)
MFall_FC$neglogPadj <- -log(MFall_FC$p.adjust)
  #Order by p.adj
ordMFall_FC <- order(MFall_FC$p.adjust)
sortMFall_FC <- MFall_FC[ordMFall_FC, ]
  #Order by logFoldChange
sortMFall_FC <- sortMFall_FC[order(sortMFall_FC$logFoldChange),]
sortMFall_FC$Description <- factor(sortMFall_FC$Description, levels = sortMFall_FC$Description)
  #Make dotplots
dotAllMFgg <- ggplot(data = sortMFall_FC, aes(x = logFoldChange, y = Description, 
                                              color = `neglogPadj`, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  labs(title="Molecular Functions in delSNORD116 (Shared Genes)", y="", x="logFoldChange")
dotAllMFgg
#For BP Up
  #Calculate fold change to use on x-axis of dotplot
BPup_FC <- summary_SimpUpBP
BPup_FC$GeneRatio <- sapply(BPup_FC$GeneRatio, function(x) eval(parse(text=x)))
BPup_FC$BgRatio <- sapply(BPup_FC$BgRatio, function(x) eval(parse(text=x)))
BPup_FC$FoldChange <- BPup_FC$GeneRatio/BPup_FC$BgRatio
BPup_FC$logFoldChange <- log2(BPup_FC$FoldChange)
BPup_FC$neglogPadj <- -log(BPup_FC$p.adjust)
  #Order by p.adj
ordBPup_FC <- order(BPup_FC$p.adjust)
sortBPup_FC <- BPup_FC[ordBPup_FC, ]
  #Only use top 25 terms
sortBPup_FC <- sortBPup_FC[c(0:25),]
  #Order by logFoldChange
sortBPup_FC <- sortBPup_FC[order(sortBPup_FC$logFoldChange),]
sortBPup_FC$Description <- factor(sortBPup_FC$Description, levels = sortBPup_FC$Description)
  #Make dotplots
dotUpBPgg <- ggplot(data = sortBPup_FC, aes(x = logFoldChange, y = Description, 
                                            color = `neglogPadj`, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  labs(title="Upregulated Biological Processes in delSNORD116", y="", x="logFoldChange")
dotUpBPgg
#For CC Up
  #Calculate fold change to use on x-axis of dotplot
CCup_FC <- summary_SimpUpCC
CCup_FC$GeneRatio <- sapply(CCup_FC$GeneRatio, function(x) eval(parse(text=x)))
CCup_FC$BgRatio <- sapply(CCup_FC$BgRatio, function(x) eval(parse(text=x)))
CCup_FC$FoldChange <- CCup_FC$GeneRatio/CCup_FC$BgRatio
CCup_FC$logFoldChange <- log2(CCup_FC$FoldChange)
CCup_FC$neglogPadj <- -log(CCup_FC$p.adjust)
  #Order by p.adj & then logFoldChange
ordCCup_FC <- order(CCup_FC$p.adjust)
sortCCup_FC <- CCup_FC[ordCCup_FC, ]
sortCCup_FC <- sortCCup_FC[order(sortCCup_FC$logFoldChange),]
sortCCup_FC$Description <- factor(sortCCup_FC$Description, levels = sortCCup_FC$Description)
  #Make dotplot
dotUpCCgg <- ggplot(data = sortCCup_FC, aes(x = logFoldChange, y = Description, 
                                           color = `neglogPadj`, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  labs(title="Upregulated Cellular Components in delSNORD116", y="", x="logFoldChange")
dotUpCCgg
#For CC Down
  #Calculate fold change to use on x-axis of dotplot
CCdown_FC <- summary_SimpDownCC
CCdown_FC$GeneRatio <- sapply(CCdown_FC$GeneRatio, function(x) eval(parse(text=x)))
CCdown_FC$BgRatio <- sapply(CCdown_FC$BgRatio, function(x) eval(parse(text=x)))
CCdown_FC$FoldChange <- CCdown_FC$GeneRatio/CCdown_FC$BgRatio
CCdown_FC$logFoldChange <- log2(CCdown_FC$FoldChange)
CCdown_FC$neglogPadj <- -log(CCdown_FC$p.adjust)
  #Order by p.adj & then logFoldChange
ordCCdown_FC <- order(CCdown_FC$p.adjust)
sortCCdown_FC <- CCdown_FC[ordCCdown_FC, ]
sortCCdown_FC <- sortCCdown_FC[order(sortCCdown_FC$logFoldChange),]
sortCCdown_FC$Description <- factor(sortCCdown_FC$Description, levels = sortCCdown_FC$Description)
  #Make dotplot
dotDownCCgg <- ggplot(data = sortCCdown_FC, aes(x = logFoldChange, y = Description, 
                                            color = `neglogPadj`, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  labs(title="Downregulated Cellular Components in delSNORD116", y="", x="logFoldChange")
dotDownCCgg
#For MF Up
  #Calculate fold change to use on x-axis of dotplot
MFup_FC <- summary_SimpUpMF
MFup_FC$GeneRatio <- sapply(MFup_FC$GeneRatio, function(x) eval(parse(text=x)))
MFup_FC$BgRatio <- sapply(MFup_FC$BgRatio, function(x) eval(parse(text=x)))
MFup_FC$FoldChange <- MFup_FC$GeneRatio/MFup_FC$BgRatio
MFup_FC$logFoldChange <- log2(MFup_FC$FoldChange)
MFup_FC$neglogPadj <- -log(MFup_FC$p.adjust)
  #Order by logFoldChange & then padj
ordMFup_FC <- order(MFup_FC$p.adjust)
sortMFup_FC <- MFup_FC[ordMFup_FC, ]
sortMFup_FC <- sortMFup_FC[order(sortMFup_FC$logFoldChange),]
sortMFup_FC$Description <- factor(sortMFup_FC$Description, levels = sortMFup_FC$Description)
  #Make dotplot
dotUpMFgg <- ggplot(data = sortMFup_FC, aes(x = logFoldChange, y = Description, 
                                            color = `neglogPadj`, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  labs(title="Upregulated Molecular Functions in delSNORD116", y="", x="logFoldChange")
dotUpMFgg
#For MF Down
  #Calculate fold change to use on x-axis of dotplot
MFdown_FC <- summary_SimpDownMF
MFdown_FC$GeneRatio <- sapply(MFdown_FC$GeneRatio, function(x) eval(parse(text=x)))
MFdown_FC$BgRatio <- sapply(MFdown_FC$BgRatio, function(x) eval(parse(text=x)))
MFdown_FC$FoldChange <- MFdown_FC$GeneRatio/MFdown_FC$BgRatio
MFdown_FC$logFoldChange <- log2(MFdown_FC$FoldChange)
MFdown_FC$neglogPadj <- -log(MFdown_FC$p.adjust)
  #Order by p.adj & then logFoldChange
ordMFdown_FC <- order(MFdown_FC$p.adjust)
sortMFdown_FC <- MFdown_FC[ordMFdown_FC, ]
sortMFdown_FC <- sortMFdown_FC[order(sortMFdown_FC$logFoldChange),]
sortMFdown_FC$Description <- factor(sortMFdown_FC$Description, levels = sortMFdown_FC$Description)
  #Make dotplot
dotDownMFgg <- ggplot(data = sortMFdown_FC, aes(x = logFoldChange, y = Description, 
                                                color = `neglogPadj`, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  labs(title="Downregulated Molecular Functions in delSNORD116", y="", x="logFoldChange")
dotDownMFgg

#Save plots as pdf
pdf("PWSexp_GO_simp.pdf", width = 10, height = 8)
dotUpBPgg
dotUpCCgg
dotDownCCgg
dotUpMFgg
dotDownMFgg
dev.off()
#Save plots as pdf
pdf("PWSexp_GO_simp_sharedgenes.pdf", width = 10, height = 8)
dotAllBPgg
dotAllCCgg
dotAllMFgg
dev.off()

##---Disease Ontology---##
#Make list of genes for universe
mart <- useMart("ensembl", host = "https://apr2018.archive.ensembl.org/", dataset = "hsapiens_gene_ensembl")
R <- data.table(univ)
names(R)[1] <- "ensembl_gene_id"
genes.table <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene"), values= R$ensembl_gene_id, mart= mart)
DOuniv <- as.character(genes.table$entrezgene)
#Break up into up & down regulated significant gene lists
DOgeneUP <- genes.table$ensembl_gene_id %in% geneUP
DOgeneUP <- subset(genes.table, DOgeneUP)
DOgeneUP <- drop_na(DOgeneUP) #Gets rid of 7 genes
DOgeneUP <- as.character(DOgeneUP$entrezgene)
DOgeneDOWN <- genes.table$ensembl_gene_id %in% geneDOWN
DOgeneDOWN <- subset(genes.table, DOgeneDOWN)
DOgeneDOWN <- drop_na(DOgeneDOWN) #Gets rid of 10 genes
DOgeneDOWN <- as.character(DOgeneDOWN$entrezgene)
#Combine up + down list
DOgeneALL <- append(DOgeneUP, DOgeneDOWN)
#Run disease ontology
DOup <- enrichDGN(gene = DOgeneUP,
                   universe = DOuniv,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)
DOdown <- enrichDGN(gene = DOgeneDOWN,
                    universe = DOuniv,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)
DOall <- enrichDGN(gene = DOgeneALL,
                    universe = DOuniv,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)

#Write results into table
summary_upDO <- data.frame(DOup)
summary_downDO <- data.frame(DOdown)
summary_allDO <- data.frame(DOall)
#Save summary tables
write.table(as.data.frame(summary_upDO), sep = "\t", file = "PWSexp_upDO.tsv", quote = FALSE)
#write.table(as.data.frame(summary_downDO), sep = "\t", file = "PWSexp_downDO.tsv", quote = FALSE)
write.table(as.data.frame(summary_allDO), sep = "\t", file = "PWSexp_allDO.tsv", quote = FALSE)

#Creating dotplots
#For DO up
  #Calculate fold change to use on x-axis of dotplot
DOup_FC <- summary_upDO
DOup_FC$GeneRatio <- sapply(DOup_FC$GeneRatio, function(x) eval(parse(text=x)))
DOup_FC$BgRatio <- sapply(DOup_FC$BgRatio, function(x) eval(parse(text=x)))
DOup_FC$FoldChange <- DOup_FC$GeneRatio/DOup_FC$BgRatio
DOup_FC$logFoldChange <- log2(DOup_FC$FoldChange)
DOup_FC$neglogPadj <- -log(DOup_FC$p.adjust)
  #Order by p.adj
ordDOup_FC <- order(DOup_FC$p.adjust)
sortDOup_FC <- DOup_FC[ordDOup_FC, ]
  #Only use top 25 terms
sortDOup_FC <- sortDOup_FC[c(0:25),]
  #Order by logFoldChange
sortDOup_FC <- sortDOup_FC[order(sortDOup_FC$logFoldChange),]
sortDOup_FC$Description <- factor(sortDOup_FC$Description, levels = sortDOup_FC$Description)
  #Make dotplots
dotDOup <- ggplot(data = sortDOup_FC, aes(x = logFoldChange, y = Description, 
                                            color = `neglogPadj`, size = Count)) + 
  geom_point() +
  scale_colour_gradient2(low = "#E2E2E2", mid = "#F99287", high = "#DC3220", midpoint = 10.2) +
  theme_bw() + 
  labs(title="Upregulated Disease Ontology in smDEL", y="", x="log2(FoldEnrichment)")
dotDOup
#For DO all
  #Calculate fold change to use on x-axis of dotplot
DOall_FC <- summary_allDO
DOall_FC$GeneRatio <- sapply(DOall_FC$GeneRatio, function(x) eval(parse(text=x)))
DOall_FC$BgRatio <- sapply(DOall_FC$BgRatio, function(x) eval(parse(text=x)))
DOall_FC$FoldChange <- DOall_FC$GeneRatio/DOall_FC$BgRatio
DOall_FC$logFoldChange <- log2(DOall_FC$FoldChange)
DOall_FC$neglogPadj <- -log(DOall_FC$p.adjust)
  #Order by p.adj
ordDOall_FC <- order(DOall_FC$p.adjust)
sortDOall_FC <- DOall_FC[ordDOall_FC, ]
  #Only use top 25 terms
sortDOall_FC <- sortDOall_FC[c(0:25),]
  #Order by logFoldChange
sortDOall_FC <- sortDOall_FC[order(sortDOall_FC$logFoldChange),]
sortDOall_FC$Description <- factor(sortDOall_FC$Description, levels = sortDOall_FC$Description)
  #Make dotplots
dotDOall <- ggplot(data = sortDOall_FC, aes(x = logFoldChange, y = Description, 
                                          color = `neglogPadj`, size = Count)) + 
  geom_point() +
  scale_colour_gradient2(low = "#E2E2E2", mid = "#37DC50", high = "#007B13", midpoint = 6.5) +
  theme_bw() + 
  labs(title="Disease Ontology in smDEL (H9 & CT2 Shared Genes)", y="", x="log2(FoldEnrichment)")
dotDOall

#Save plots as pdf
pdf("PWSexp_DO.pdf", width = 10, height = 8)
dotDOup
dev.off()
pdf("PWSexp_DOall.pdf", width = 10, height = 8)
dotDOall
dev.off()

##---Gene-Concept Network plots---##
cnetBP1 <- cnetplot(
  allGOBP,
  showCategory = 5,
  foldChange = foldChange,
  layout = "star",
  colorEdge = TRUE,
  circular = FALSE,
  node_label = "all",
  cex_category = 1,
  cex_gene = 1,
  node_label_size = NULL,
  cex_label_category = 1,
  cex_label_gene = 1,
  shadowtext = "all")
cnetBP1
circBP1 <- cnetplot(
  allGOBP,
  showCategory = 10,
  foldChange = foldChange,
  layout = "circle",
  colorEdge = TRUE,
  circular = TRUE,
  node_label = "all",
  cex_category = 1,
  cex_gene = 1,
  node_label_size = NULL,
  cex_label_category = 1,
  cex_label_gene = 1)
circBP1
cnetBP1simp <- cnetplot(
  SimpAllBP,
  showCategory = 5,
  foldChange = foldChange,
  layout = "star",
  colorEdge = TRUE,
  circular = FALSE,
  node_label = "all",
  cex_category = 1,
  cex_gene = 1,
  node_label_size = NULL,
  cex_label_category = 1,
  cex_label_gene = 1,
  shadowtext = "all")
cnetBP1simp
#Convert gene ID to symbol for DO
edox <- setReadable(DOall, 'org.Hs.eg.db', 'ENTREZID')
#Plot
categories <- c("Mental Retardation", "Short palm", "Short toe", "Hypospadias", "Intervertebral Disc Degeneration")
cnetDO <- cnetplot(
  edox,
  showCategory = categories,
  foldChange = foldChange,
  layout = "dh",
  colorEdge = TRUE,
  circular = FALSE,
  node_label = "all",
  cex_category = 1,
  cex_gene = 1,
  node_label_size = NULL,
  cex_label_category = 1,
  cex_label_gene = 1,
  shadowtext = "all")
cnetDO
#Save plots as pdf
pdf("PWSexp_cnetplots.pdf", width = 10, height = 8)

dev.off()

##---Enrichment Map---##
allGOBPe <- pairwise_termsim(allGOBP)
emap1 <- emapplot(
  allGOBPe,
  split = 'ONTOLOGY',
  showCategory = 50,
  layout = NULL,
  color = "p.adjust",
  node_label = "category",
  cex_label_category = 0.5,
  repel = TRUE)
emap1
#Save plots as pdf
pdf("PWSexp_emapplots.pdf", width = 10, height = 8)
emap1
dev.off()
