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
library("gplots")

#Set working directory
directory <- "../Cotney_Lab/PWS_RNASeq/STAR/lgDEL"
setwd(directory)

#Read in files of shared genes
geneUP <- read.delim("VennDiagramlgDEL_IntersectUpResultsENSEMBL.txt", sep="\t", head=F)
geneDOWN <- read.delim("VennDiagramlgDEL_IntersectDownResultsENSEMBL.txt", sep="\t", head=F)

#Load dds object
load("PWSlgDEL_DESeqSTAR.RData")
#Make list of genes for universe (all genes from alignment)
rownames(dds) <- gsub('\\..+$', '', rownames(dds))
univ <- row.names(dds)

#Break up into up & down regulated gene lists
geneUP <- geneUP$V1
geneDOWN <- geneDOWN$V1
#Create combined list for both up & down regulated lists
geneALL <- append(geneUP, geneDOWN)

#Read in combined results object
ResGeneID_ALL <- read.csv(file = 'CombinedH9CT2_Results_DESeqSTAR.csv', header = TRUE, row.names = "X")
#Calculate average log2FoldChange
ResGeneID_ALL$log2FoldChange <- rowMeans(ResGeneID_ALL[,c(3,51)])

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

#Use dplyr to make foldchange object
foldChange <- resALL_data %>% 
  dplyr::select(symbol, log2FoldChange) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(log2FoldChange=mean(log2FoldChange))
names(foldChange)[1] <- "Gene"
foldChange <- deframe(foldChange)
#Make binerized foldChange object
foldChange2 <- resALL_data %>% 
  dplyr::select(symbol, log2FoldChange) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(log2FoldChange=mean(log2FoldChange))
foldChange2 <- foldChange2 %>% mutate(log2FoldChange = ifelse(log2FoldChange < 0, -1, log2FoldChange))
foldChange2 <- foldChange2 %>% mutate(log2FoldChange = ifelse(log2FoldChange > 0, 1, log2FoldChange))
names(foldChange2)[1] <- "Gene"
foldChange2 <- deframe(foldChange2)

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
write.table(as.data.frame(summary_upBP), sep = "\t", file = "PWSexplgDEL_upBP.tsv", quote = FALSE)
write.table(as.data.frame(summary_downBP), sep = "\t", file = "PWSexplgDEL_downBP.tsv", quote = FALSE)
write.table(as.data.frame(summary_upCC), sep = "\t", file = "PWSexplgDEL_upCC.tsv", quote = FALSE)
write.table(as.data.frame(summary_downCC), sep = "\t", file = "PWSexplgDEL_downCC.tsv", quote = FALSE)
write.table(as.data.frame(summary_upMF), sep = "\t", file = "PWSexplgDEL_upMF.tsv", quote = FALSE)
write.table(as.data.frame(summary_downMF), sep = "\t", file = "PWSexplgDEL_downMF.tsv", quote = FALSE)
write.table(as.data.frame(summary_allBP), sep = "\t", file = "PWSexplgDEL_allBP.tsv", quote = FALSE)
write.table(as.data.frame(summary_allCC), sep = "\t", file = "PWSexplgDEL_allCC.tsv", quote = FALSE)
write.table(as.data.frame(summary_allMF), sep = "\t", file = "PWSexplgDEL_allMF.tsv", quote = FALSE)

#Simplify combining similar terms
SimpUpBP <- simplify(GOupBP, cutoff = 0.8, by = "p.adjust", select_fun = min)
summary_SimpUpBP <- data.frame(SimpUpBP)
SimpDownBP <- simplify(GOdownBP, cutoff = 0.8, by = "p.adjust", select_fun = min)
summary_SimpDownBP <- data.frame(SimpDownBP)
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
write.table(as.data.frame(summary_SimpAllBP), sep = "\t", file = "PWSexplgDEL_allBP_simplified.tsv", quote = FALSE)
write.table(as.data.frame(summary_SimpAllCC), sep = "\t", file = "PWSexplgDEL_allCC_simplified.tsv", quote = FALSE)
write.table(as.data.frame(summary_SimpAllMF), sep = "\t", file = "PWSexplgDEL_allMF_simplified.tsv", quote = FALSE)
write.table(as.data.frame(summary_SimpUpBP), sep = "\t", file = "PWSexplgDEL_upBP_simplified.tsv", quote = FALSE)
write.table(as.data.frame(summary_SimpUpCC), sep = "\t", file = "PWSexplgDEL_upCC_simplified.tsv", quote = FALSE)
write.table(as.data.frame(summary_SimpUpMF), sep = "\t", file = "PWSexplgDEL_upMF_simplified.tsv", quote = FALSE)
write.table(as.data.frame(summary_SimpDownBP), sep = "\t", file = "PWSexplgDEL_downBP_simplified.tsv", quote = FALSE)
write.table(as.data.frame(summary_SimpDownCC), sep = "\t", file = "PWSexplgDEL_downCC_simplified.tsv", quote = FALSE)
write.table(as.data.frame(summary_SimpDownMF), sep = "\t", file = "PWSexplgDEL_downMF_simplified.tsv", quote = FALSE)

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
  #Order by logFoldChange
sortBPall_FC <- sortBPall_FC[order(sortBPall_FC$logFoldChange),]
sortBPall_FC$Description <- factor(sortBPall_FC$Description, levels = sortBPall_FC$Description)
  #Make dotplots
dotAllBPgg <- ggplot(data = sortBPall_FC, aes(x = logFoldChange, y = Description, 
                                            color = `neglogPadj`, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  labs(title="Biological Processes in lgDEL (Shared Genes)", y="", x="logFoldChange")
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
  labs(title="Cellular Components in lgDEL (Shared Genes)", y="", x="logFoldChange")
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
  labs(title="Molecular Functions in lgDEL (Shared Genes)", y="", x="logFoldChange")
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
  labs(title="Upregulated Biological Processes in lgDEL", y="", x="logFoldChange")
dotUpBPgg
#For BP Down
  #Calculate fold change to use on x-axis of dotplot
BPdown_FC <- summary_SimpDownBP
BPdown_FC$GeneRatio <- sapply(BPdown_FC$GeneRatio, function(x) eval(parse(text=x)))
BPdown_FC$BgRatio <- sapply(BPdown_FC$BgRatio, function(x) eval(parse(text=x)))
BPdown_FC$FoldChange <- BPdown_FC$GeneRatio/BPdown_FC$BgRatio
BPdown_FC$logFoldChange <- log2(BPdown_FC$FoldChange)
BPdown_FC$neglogPadj <- -log(BPdown_FC$p.adjust)
  #Order by p.adj
ordBPdown_FC <- order(BPdown_FC$p.adjust)
sortBPdown_FC <- BPdown_FC[ordBPdown_FC, ]
  #Order by logFoldChange
sortBPdown_FC <- sortBPdown_FC[order(sortBPdown_FC$logFoldChange),]
sortBPdown_FC$Description <- factor(sortBPdown_FC$Description, levels = sortBPdown_FC$Description)
  #Make dotplots
dotDownBPgg <- ggplot(data = sortBPdown_FC, aes(x = logFoldChange, y = Description, 
                                            color = `neglogPadj`, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  labs(title="Downregulated Biological Processes in lgDEL", y="", x="logFoldChange")
dotDownBPgg
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
  labs(title="Upregulated Cellular Components in lgDEL", y="", x="logFoldChange")
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
  labs(title="Downregulated Cellular Components in lgDEL", y="", x="logFoldChange")
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
  labs(title="Upregulated Molecular Functions in lgDEL", y="", x="logFoldChange")
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
  labs(title="Downregulated Molecular Functions in lgDEL", y="", x="logFoldChange")
dotDownMFgg

#Save plots as pdf
pdf("PWSexplgDEL_GO_simp.pdf", width = 10, height = 8)
dotUpBPgg
dotDownBPgg
dotUpCCgg
dotDownCCgg
dotUpMFgg
dotDownMFgg
dev.off()
#Save plots as pdf
pdf("PWSexplgDEL_GO_simp_sharedgenes.pdf", width = 10, height = 8)
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
DOgeneUP <- drop_na(DOgeneUP) #Gets rid of 31 genes
DOgeneUP <- as.character(DOgeneUP$entrezgene)
DOgeneDOWN <- genes.table$ensembl_gene_id %in% geneDOWN
DOgeneDOWN <- subset(genes.table, DOgeneDOWN)
DOgeneDOWN <- drop_na(DOgeneDOWN) #Gets rid of 42 genes
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
write.table(as.data.frame(summary_upDO), sep = "\t", file = "PWSexplgDEL_upDO.tsv", quote = FALSE)
write.table(as.data.frame(summary_downDO), sep = "\t", file = "PWSexplgDEL_downDO.tsv", quote = FALSE)
write.table(as.data.frame(summary_allDO), sep = "\t", file = "PWSexplgDEL_allDO.tsv", quote = FALSE)

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
  #Order by logFoldChange
sortDOup_FC <- sortDOup_FC[order(sortDOup_FC$logFoldChange),]
sortDOup_FC$Description <- factor(sortDOup_FC$Description, levels = sortDOup_FC$Description)
  #Make dotplots
dotDOup <- ggplot(data = sortDOup_FC, aes(x = logFoldChange, y = Description, 
                                            color = `neglogPadj`, size = Count)) + 
  geom_point() +
  scale_colour_gradient2(low = "#E2E2E2", mid = "#F99287", high = "#DC3220", midpoint = 3.77) +
  theme_bw() + 
  labs(title="Upregulated Disease Ontology in lgDEL", y="", x="log2(FoldEnrichment)")
dotDOup
#For DO down
  #Calculate fold change to use on x-axis of dotplot
DOdown_FC <- summary_downDO
DOdown_FC$GeneRatio <- sapply(DOdown_FC$GeneRatio, function(x) eval(parse(text=x)))
DOdown_FC$BgRatio <- sapply(DOdown_FC$BgRatio, function(x) eval(parse(text=x)))
DOdown_FC$FoldChange <- DOdown_FC$GeneRatio/DOdown_FC$BgRatio
DOdown_FC$logFoldChange <- log2(DOdown_FC$FoldChange)
DOdown_FC$neglogPadj <- -log(DOdown_FC$p.adjust)
  #Order by p.adj
ordDOdown_FC <- order(DOdown_FC$p.adjust)
sortDOdown_FC <- DOdown_FC[ordDOdown_FC, ]
  #Order by logFoldChange
sortDOdown_FC <- sortDOdown_FC[order(sortDOdown_FC$logFoldChange),]
sortDOdown_FC$Description <- factor(sortDOdown_FC$Description, levels = sortDOdown_FC$Description)
  #Make dotplots
dotDOdown <- ggplot(data = sortDOdown_FC, aes(x = logFoldChange, y = Description, 
                                          color = `neglogPadj`, size = Count)) + 
  geom_point() +
  scale_colour_gradient2(low = "#E2E2E2", mid = "#6EACFF", high = "#004CB3", midpoint = 10) +
  theme_bw() + 
  labs(title="Downregulated Disease Ontology in lgDEL", y="", x="log2(FoldEnrichment)")
dotDOdown
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
  #Order by logFoldChange
sortDOall_FC <- sortDOall_FC[order(sortDOall_FC$logFoldChange),]
sortDOall_FC$Description <- factor(sortDOall_FC$Description, levels = sortDOall_FC$Description)
  #Make dotplots
dotDOall <- ggplot(data = sortDOall_FC, aes(x = logFoldChange, y = Description, 
                                          color = `neglogPadj`, size = Count)) + 
  geom_point() +
  scale_colour_gradient2(low = "#E2E2E2", mid = "#37DC50", high = "#007B13", midpoint = 5.2) +
  theme_bw() + 
  labs(title="Disease Ontology in lgDEL (Shared Genes)", y="", x="log2(FoldEnrichment)")
dotDOall
#Save plots as pdf
pdf("PWSexplgDEL_DO.pdf", width = 10, height = 8)
dotDOup
dotDOdown
dev.off()
pdf("PWSexplgDEL_DOall.pdf", width = 10, height = 8)
dotDOall
dev.off()
pdf("skinnylgDEL_DOdown.pdf", width = 5, height = 8)
dotDOdown
dev.off()

##---Gene-Concept Network plots---##
cnetMF <- cnetplot(
  SimpAllMF,
  showCategory = 5,
  foldChange = foldChange2,
  colorEdge = TRUE,
  circular = FALSE,
  node_label = "all",
  cex_category = 1,
  cex_gene = 1,
  node_label_size = NULL,
  cex_label_category = 1,
  cex_label_gene = 0.6,
  shadowtext = "all") + scale_color_gradient2(low = "#005AB5", mid = "#EFEFEF",  high = "#DC3220")
cnetMF
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
edox <- setReadable(DOdown, 'org.Hs.eg.db', 'ENTREZID')
categories <- c("Obesity, Abdominal", "Acromicric Dysplasia", "Narrow palm", "Hand deformities", "Depressed nasal ridge", "Delayed Puberty", "Abnormality of the genital system")
#Plot
cnetDO <- cnetplot(
  edox,
  showCategory = categories,
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
cnetDO
#Save plot as pdf
pdf("DOdown_cnetplot.pdf", width = 10, height = 8)
cnetDO
dev.off()
pdf("MFall_cnetplot.pdf", width = 15, height = 8)
cnetMF
dev.off()
pdf("MFall_cnetplot_binerized.pdf", width = 15, height = 8)
cnetMF
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

##Determine expression of ribosomal genes from MF/disease ontologies
  #To find list of all genes under specific GO term
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'go_id', 'name_1006', 'definition_1006'), filters = 'go', values = 'GO:0003735', mart = ensembl)
gene_ont_MFribo <- gene.data[which(gene.data$go_id == "GO:0003735"),]
gene_ont_MFribo$hgnc_symbol[gene_ont_MFribo$hgnc_symbol == ""] <- NA
gene_ont_MFribo <- gene_ont_MFribo %>% drop_na(hgnc_symbol)
#Drop top two MT genes
MFribo_list <- unique(gene_ont_MFribo$hgnc_symbol)
  write.table(MFribo_list, "allGO0003735_symbols.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
  #To find list of genes present from list of shared genes
myRibo_list <- summary_allMF[1,8]
  write.table(myRibo_list, "myGenes_GO0003735_symbols.txt", quote=FALSE, row.names = FALSE, col.names = FALSE) #Put in Notepad++ and replace "/" with "\n"
myRibo_list <- scan("myGenes_GO0003735_symbols.txt", what="", sep="\n")
#Make overall list without my genes
MFribo_minusMine <- MFribo_list %>% filter(!MFribo_list$V1 %in% myRibo_list$V1)
  write.table(MFribo_minusMine, "GO0003735_symbols_minusMyGenes.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
#Read in & format GTEx table (converted .gct file downloaded from GTEx to .txt file for reading in)
GTEx <- read.delim("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.txt", sep="\t", head=T)
colnames(GTEx) <- c("GeneID", "GeneName", "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)",	"Adrenal Gland",	"Artery - Aorta",	"Artery - Coronary",	"Artery - Tibial",	"Bladder",	"Brain - Amygdala",	"Brain - Anterior cingulate cortex (BA24)",	"Brain - Caudate (basal ganglia)",	"Brain - Cerebellar Hemisphere",	"Brain - Cerebellum",	"Brain - Cortex",	"Brain - Frontal Cortex (BA9)",	"Brain - Hippocampus",	"Brain - Hypothalamus",	"Brain - Nucleus accumbens (basal ganglia)",	"Brain - Putamen (basal ganglia)",	"Brain - Spinal cord (cervical c-1)",	"Brain - Substantia nigra",	"Breast - Mammary Tissue",	"Cells - Cultured fibroblasts",	"Cells - EBV-transformed lymphocytes",	"Cervix - Ectocervix",	"Cervix - Endocervix",	"Colon - Sigmoid",	"Colon - Transverse",	"Esophagus - Gastroesophageal Junction",	"Esophagus - Mucosa",	"Esophagus - Muscularis",	"Fallopian Tube",	"Heart - Atrial Appendage",	"Heart - Left Ventricle",	"Kidney - Cortex",	"Kidney - Medulla",	"Liver",	"Lung",	"Minor Salivary Gland",	"Muscle - Skeletal",	"Nerve - Tibial",	"Ovary",	"Pancreas",	"Pituitary",	"Prostate",	"Skin - Not Sun Exposed (Suprapubic)",	"Skin - Sun Exposed (Lower leg)",	"Small Intestine - Terminal Ileum",	"Spleen",	"Stomach",	"Testis",	"Thyroid",	"Uterus",	"Vagina",	"Whole Blood")
#Parse genes from GTEx table
myRibo_GTEx <- GTEx[(GTEx$GeneName) %in% myRibo_list, ]
row.names(myRibo_GTEx) <- myRibo_GTEx$GeneName
myRibo_GTEx <- myRibo_GTEx[c(3:56)]
#Plot
heatmapMat_myRibo <- data.matrix(myRibo_GTEx)
hr_myRibo <- hclust(as.dist(1-cor(t(heatmapMat_myRibo), method="spearman")), method="complete")
hc_myRibo <- hclust(as.dist(1-cor(heatmapMat_myRibo, method="spearman")), method="complete")
heatmap.2( heatmapMat_myRibo, key=T, scale="row", 
           trace="none", Rowv=as.dendrogram(hr_myRibo), Colv=as.dendrogram(hc_myRibo), dendrogram = "row", cexRow = .4, cexCol = .7,
           col = colorRampPalette(c("blue", "white", "red"))(n = 1000))
pdf("myRiboHeatmap_skinny.pdf", width = 9, height = 4)
heatmap.2( heatmapMat_myRibo, key=T, scale="row", 
           trace="none", Rowv=as.dendrogram(hr_myRibo), Colv=as.dendrogram(hc_myRibo), dendrogram = "row", cexRow = .4, cexCol = .7,
           col = colorRampPalette(c("blue", "white", "red"))(n = 1000))
dev.off()
pdf("myRiboHeatmap_reg.pdf", width = 10, height = 8)
heatmap.2( heatmapMat_myRibo, key=T, scale="row", 
           trace="none", Rowv=as.dendrogram(hr_myRibo), Colv=as.dendrogram(hc_myRibo), dendrogram = "row", cexRow = .4, cexCol = .7,
           col = colorRampPalette(c("blue", "white", "red"))(n = 1000))
dev.off()

