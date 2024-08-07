##---DESeq---##
#Load relevant libraries
library("DESeq2")
library("data.table")
library("EnhancedVolcano")
library("dplyr")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("biomaRt")
library("pheatmap")
library("RColorBrewer")
library("dendextend")
library("gplots")
library("ggplot2")
library("ggVennDiagram")
library("ggpubr")

#Set working directory
directory <- "../Cotney_Lab/PWS_RNASeq/HISAT2/lgDEL"
setwd(directory)

#Import file
counts <- read.delim("PWSRNAseq_HISAT2_featureCounts.txt", sep="\t", head=T, skip = 1, row.names = "Geneid")

#Simplify column names
colnames(counts) <- c("chr", "start", "end", "strand", "length", "CT2_1", "CT2_2", "CT2_3", "CT2_4", "CT2_5", "CT2_6", "CT2_lgDEL_1", "CT2_lgDEL_2", "CT2_lgDEL_3", "CT2_lgDEL_4", "CT2_lgDEL_5", "CT2_lgDEL_6", "CT2_smDEL_1", "CT2_smDEL_2", "CT2_smDEL_3", "CT2_smDEL_4", "CT2_smDEL_5", "CT2_smDEL_6", "H9_1", "H9_2", "H9_3", "H9_4", "H9_5", "H9_6", "H9_lgDEL_1", "H9_lgDEL_2", "H9_lgDEL_3", "H9_lgDEL_4", "H9_lgDEL_5", "H9_lgDEL_6", "H9_smDEL_1", "H9_smDEL_2", "H9_smDEL_3", "H9_smDEL_4", "H9_smDEL_5")
#Get rid of first 5 columns
counts <- counts[ ,6:ncol(counts)]
#Turn into matrix
counts <- as.matrix(counts)

#Load in metadata
colData <- as.data.frame(read.table("Metadata.txt", sep="\t", head=T, row.names = "ID"))
colData$Condition <- as.factor(colData$Condition)
colData$Background_Condition <- as.factor(paste(colData$Genetic.Background, colData$Condition, sep="_"))

#Create DESeq objects
ddstotal <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~Background_Condition)
ddstotal

#Trim matrix, remove low counts
summary(rowSums(counts(ddstotal)))
#Filter using cutoff (can modify if necessary)
dds <- ddstotal [rowSums(counts(ddstotal)) > 1 ,]

#Run DESeq on data
dds <- DESeq(dds)
  #To load
  load("PWSlgDEL_DESeqHISAT2.RData")
#Get log values
rld <- rlog(dds, blind=FALSE)
head(assay(rld))
hist(assay(rld))

##---PCA Plot---##
plotPCA(rld, intgroup = "Condition")
plotPCA(rld, intgroup = "Genetic.Background")
plotPCA(rld, intgroup = "Sequencing.Batch")
plotPCA(rld, intgroup = "Differentiation.Start.Date")
plotPCA(rld, intgroup = "Neuron.Collection.Date")
plotPCA(rld, intgroup = "Plate.Type")
plotPCA(rld, intgroup = "RNA.Extraction.Date")
plotPCA(rld, intgroup = "RIN.value")
plotPCA(rld, intgroup = "RNA.Concentration..ng.ul.")

#Save as pdf
pdf("PWSlgDEL_Exp_PCAs.pdf", width = 10, height = 8)
plotPCA(rld, intgroup = "Condition")
plotPCA(rld, intgroup = "Genetic.Background")
plotPCA(rld, intgroup = "Sequencing.Batch")
plotPCA(rld, intgroup = "Differentiation.Start.Date")
plotPCA(rld, intgroup = "Neuron.Collection.Date")
plotPCA(rld, intgroup = "Plate.Type")
plotPCA(rld, intgroup = "RNA.Extraction.Date")
plotPCA(rld, intgroup = "RIN.value")
plotPCA(rld, intgroup = "RNA.Concentration..ng.ul.")
dev.off()

##---Data "Housekeeping"---##
#Removing decimal in gene name
rownames(rld) <- gsub('\\..+$', '', rownames(rld))
rownames(dds) <- gsub('\\..+$', '', rownames(dds))
#Check that it worked
counts(dds)
#Run results
resH9lgDEL <- results(dds, contrast = c("Background_Condition","H9_lgDEL","H9_WT"))
resH9lgDEL
resCT2lgDEL <- results(dds, contrast = c("Background_Condition","CT2_lgDEL","CT2_WT"))
resCT2lgDEL
#Merge results table & counts table
resDataH9lgDEL <- merge(as.data.frame(resH9lgDEL), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resDataCT2lgDEL <- merge(as.data.frame(resCT2lgDEL), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resDataH9lgDEL)[1] <- "Gene"
names(resDataCT2lgDEL)[1] <- "Gene"
rownames(resDataH9lgDEL) <- resDataH9lgDEL$Gene
rownames(resDataCT2lgDEL) <- resDataCT2lgDEL$Gene
#Change ESMBL ID's to gene names
mart <- useMart("ensembl", host = "https://dec2017.archive.ensembl.org/", dataset = "hsapiens_gene_ensembl")
  #For H9lgDEL
R1 <- data.table(resDataH9lgDEL)
rownames(R1) <- R1$Gene
genes.table_H9lgDEL <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","description","gene_biotype","chromosome_name","start_position","end_position","strand"), values= R1$Gene, mart= mart)
  #For CT2lgDEL
R2 <- data.table(resDataCT2lgDEL)
rownames(R2) <- R2$Gene
genes.table_CT2lgDEL <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","description","gene_biotype","chromosome_name","start_position","end_position","strand"), values= R2$Gene, mart= mart)
#Order lists by ENSEMBL ID so they merge correctly
resOrderedgene_H9lgDEL <- R1[order(R1$Gene),]
names(genes.table_H9lgDEL)[1] <- "Gene"
resOrderedgene_CT2lgDEL <- R2[order(R2$Gene),]
names(genes.table_CT2lgDEL)[1] <- "Gene"
#Merge all things
ResGeneID_H9lgDEL <- merge(as.data.frame(resOrderedgene_H9lgDEL), as.data.frame(genes.table_H9lgDEL), by="Gene", sort=TRUE, all = FALSE)
rownames(ResGeneID_H9lgDEL) <- ResGeneID_H9lgDEL$Gene
ResGeneID_CT2lgDEL <- merge(as.data.frame(resOrderedgene_CT2lgDEL), as.data.frame(genes.table_CT2lgDEL), by="Gene", sort=TRUE, all = FALSE)
rownames(ResGeneID_CT2lgDEL) <- ResGeneID_CT2lgDEL$Gene
#Fill in blank external gene names
  #For H9lgDEL
ResGeneID_H9lgDEL$external_gene_name <- as.character(ResGeneID_H9lgDEL$external_gene_name)
ResGeneID_H9lgDEL$external_gene_name[ResGeneID_H9lgDEL$external_gene_name==""] <- NA
ResGeneID_H9lgDEL$external_gene_name <- as.factor(ResGeneID_H9lgDEL$external_gene_name)
GeneSymbolDE_H9lgDEL <- mapIds(org.Hs.eg.db, keys = ResGeneID_H9lgDEL$Gene, column = "SYMBOL", keytype = "ENSEMBL", multivals= "first")
for (i in c(1:length(GeneSymbolDE_H9lgDEL))){
  
  if (is.na(GeneSymbolDE_H9lgDEL[i])){
    
    GeneSymbolDE_H9lgDEL[i]<-gsub("\\..*","",ResGeneID_H9lgDEL$Gene[i])
    
  }
  
}
ResGeneID_H9lgDEL$external_gene_name <- GeneSymbolDE_H9lgDEL
  #For CT2lgDEL
ResGeneID_CT2lgDEL$external_gene_name <- as.character(ResGeneID_CT2lgDEL$external_gene_name)
ResGeneID_CT2lgDEL$external_gene_name[ResGeneID_CT2lgDEL$external_gene_name==""] <- NA
ResGeneID_CT2lgDEL$external_gene_name <- as.factor(ResGeneID_CT2lgDEL$external_gene_name)
GeneSymbolDE_CT2lgDEL <- mapIds(org.Hs.eg.db, keys = ResGeneID_CT2lgDEL$Gene, column = "SYMBOL", keytype = "ENSEMBL", multivals= "first")
for (i in c(1:length(GeneSymbolDE_CT2lgDEL))){
  
  if (is.na(GeneSymbolDE_CT2lgDEL[i])){
    
    GeneSymbolDE_CT2lgDEL[i]<-gsub("\\..*","",ResGeneID_CT2lgDEL$Gene[i])
    
  }
  
}
ResGeneID_CT2lgDEL$external_gene_name <- GeneSymbolDE_CT2lgDEL
#Merge both results objects
ResGeneID_ALL <- merge(as.data.frame(ResGeneID_H9lgDEL), as.data.frame(ResGeneID_CT2lgDEL), by="Gene", sort=TRUE, all = FALSE)
rownames(ResGeneID_ALL) <- ResGeneID_ALL$Gene
#Drop unnecessary columns
ResGeneID_ALL <- ResGeneID_ALL[-c(56:97)]
#Simplify column names
colnames(ResGeneID_ALL) <- c("Gene", "baseMean_H9", "log2FoldChange_H9", "lfcSE_H9", "stat_H9", "pvalue_H9", "padj_H9", "CT2_1", "CT2_2", "CT2_3", "CT2_4", "CT2_5", "CT2_6", "CT2_lgDEL_1", "CT2_lgDEL_2", "CT2_lgDEL_3", "CT2_lgDEL_4", "CT2_lgDEL_5", "CT2_lgDEL_6", "CT2_smDEL_1", "CT2_smDEL_2", "CT2_smDEL_3", "CT2_smDEL_4", "CT2_smDEL_5", "CT2_smDEL_6", "H9_1", "H9_2", "H9_3", "H9_4", "H9_5", "H9_6", "H9_lgDEL_1", "H9_lgDEL_2", "H9_lgDEL_3", "H9_lgDEL_4", "H9_lgDEL_5", "H9_lgDEL_6", "H9_smDEL_1", "H9_smDEL_2", "H9_smDEL_3", "H9_smDEL_4", "H9_smDEL_5", "external_gene_name", "description", "gene_biotype", "chromosome_name", "start_position", "end_position", "strand", "baseMean_CT2", "log2FoldChange_CT2", "lfcSE_CT2", "stat_CT2", "pvalue_CT2", "padj_CT2")
#Save summary table
write.csv(as.data.frame(ResGeneID_ALL), file = "CombinedH9CT2lgDEL_Results_DESeqHISAT2.csv")
  #To read in
  ResGeneID_ALL <- read.csv(file = 'CombinedH9CT2lgDEL_Results_DESeqHISAT2.csv', header = TRUE, row.names = "X")
  
##---Boxplots---##
#Create boxplots of gene expression for 15q11-q13 region
#Import csv made via UCSC genome browser
chr15q <- read.csv("chr15q11-q13_gencodev25.csv")
#Set column names
colnames(chr15q) <- c("transcript_id", "name2", "ensembl_gene_id", "gene_name")
#Remove decimal in gene_id
chr15q$ensembl_gene_id <- gsub('\\..+$', '', chr15q$ensembl_gene_id)
#Create list of unique ENSEMBL genes
chr15q_genelist <- unique(chr15q$ensembl_gene_id)
#Parse DE results file by desired gene list
chr15qdat <- ResGeneID_ALL[(ResGeneID_ALL$Gene) %in% chr15q_genelist, ]
#Put genes in correct order
chr15qOrder <- chr15qdat[order(chr15qdat$start_position),]
#Add 1 to all counts --> Change columns based on samples!
plusone <- function(x){
  return (x + 1)
}
chr15qOrder_plus1 <- chr15qOrder
chr15qOrder_plus1[,c(8:42)]<- data.frame(lapply(chr15qOrder_plus1[,c(8:42)],plusone))
#Calculate mean of WT values
chr15qOrder_plus1$wtmed_H9 <- rowMeans(chr15qOrder_plus1[,c(26:31)])
chr15qOrder_plus1$wtmed_CT2 <- rowMeans(chr15qOrder_plus1[,c(8:13)])
#Divide large deletion samples by WT mean
chr15qOrder_plus1[,c(32:37)] <- chr15qOrder_plus1[,c(32:37)]/chr15qOrder_plus1$wtmed_H9
chr15qOrder_plus1[,c(14:19)] <- chr15qOrder_plus1[,c(14:19)]/chr15qOrder_plus1$wtmed_CT2
#Log transform (same columns as division step above)
chr15qOrder_plus1_log2 <- chr15qOrder_plus1
chr15qOrder_plus1_log2[,c(32:37)] <- log2(chr15qOrder_plus1_log2[,c(32:37)])
chr15qOrder_plus1_log2[,c(14:19)] <- log2(chr15qOrder_plus1_log2[,c(14:19)])
#Transpose
chr15qOrderH9lgDEL_plus1_log2_t <- t(chr15qOrder_plus1_log2[,c(32:37)])
chr15qOrderCT2lgDEL_plus1_log2_t <- t(chr15qOrder_plus1_log2[,c(14:19)])
#Make the gene name the title of the columns
colnames(chr15qOrderH9lgDEL_plus1_log2_t) <- chr15qOrder_plus1_log2$external_gene_name
colnames(chr15qOrderCT2lgDEL_plus1_log2_t) <- chr15qOrder_plus1_log2$external_gene_name
  #Save final 15q gene expression table
  write.csv(as.data.frame(chr15qOrder_plus1_log2), file = "CombinedlgDEL_DESeqHISAT2_15qgenes_gencodev25.csv")
#Filter DE vs nonDE genes for chr15
chr15qde_H9lgDEL <- chr15qdat[which(chr15qdat$padj_H9 < 0.05),]
chr15qde_CT2lgDEL <- chr15qdat[which(chr15qdat$padj_CT2 < 0.05),]
chr15qNONde_H9lgDEL <- chr15qdat[!(row.names(chr15qdat) %in% row.names(chr15qde_H9lgDEL)),]
chr15qNONde_CT2lgDEL <- chr15qdat[!(row.names(chr15qdat) %in% row.names(chr15qde_CT2lgDEL)),]
#Plot with ggplot
  #Reshape data
chr_meltH9lgDEL <- reshape2::melt(chr15qOrderH9lgDEL_plus1_log2_t)
head(chr_meltH9lgDEL)
chr_meltCT2lgDEL <- reshape2::melt(chr15qOrderCT2lgDEL_plus1_log2_t)
head(chr_meltCT2lgDEL)
  #Set factors
jitter <- position_jitter(width = 0.18, height = 0.1)
chr15qOrder_plus1_log2$external_gene_name <- factor(chr15qOrder_plus1_log2$external_gene_name, levels = chr15qOrder_plus1_log2$external_gene_name)
  #Plot
H9_plot <- ggplot(chr_meltH9lgDEL, aes(x = Var2, y = value)) + 
  geom_boxplot(aes(color = Var2), filter(chr_meltH9lgDEL, Var2 %in% chr15qde_H9lgDEL$external_gene_name), color = "#AB5000", outlier.shape = NA) + geom_boxplot(aes(color = Var2), filter(chr_meltH9lgDEL, Var2 %in% chr15qNONde_H9lgDEL$external_gene_name), color = "black", outlier.shape = NA) +
  geom_point(aes(color = Var2), filter(chr_meltH9lgDEL, Var2 %in% chr15qde_H9lgDEL$external_gene_name), size = 1, color = "#AB5000", alpha = 0.5, position = jitter) + geom_point(aes(color = Var2), filter(chr_meltH9lgDEL, Var2 %in% chr15qNONde_H9lgDEL$external_gene_name), size = 1, color = "black", alpha = 0.5, position = jitter) +
  scale_x_discrete(limits = levels(chr15qOrder_plus1_log2$external_gene_name)) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6), panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "lightgray"), panel.background = element_rect(fill = "white"), axis.line = element_line(color="black")) +
  labs(title="Chr15q11-q13 Gene Expression: H9lgDEL", x ="Gene", y = "log2(foldChange)")
H9_plot
CT2_plot <- ggplot(chr_meltCT2lgDEL, aes(x = Var2, y = value)) + 
  geom_boxplot(aes(color = Var2), filter(chr_meltCT2lgDEL, Var2 %in% chr15qde_CT2lgDEL$external_gene_name), color = "#5D3A9B", outlier.shape = NA) + geom_boxplot(aes(color = Var2), filter(chr_meltCT2lgDEL, Var2 %in% chr15qNONde_CT2lgDEL$external_gene_name), color = "black", outlier.shape = NA) +
  geom_point(aes(color = Var2), filter(chr_meltCT2lgDEL, Var2 %in% chr15qde_CT2lgDEL$external_gene_name), size = 1, color = "#5D3A9B", alpha = 0.5, position = jitter) + geom_point(aes(color = Var2), filter(chr_meltCT2lgDEL, Var2 %in% chr15qNONde_CT2lgDEL$external_gene_name), size = 1, color = "black", alpha = 0.5, position = jitter) +
  scale_x_discrete(limits = levels(chr15qOrder_plus1_log2$external_gene_name)) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6), panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "lightgray"), panel.background = element_rect(fill = "white"), axis.line = element_line(color="black")) +
  labs(title="Chr15q11-q13 Gene Expression: CT2lgDEL", x ="Gene", y = "log2(foldChange)")
CT2_plot
#Save as pdf
pdf("Chr15q11-q13GeneExp_lgDEL_gencodev25.pdf", width = 10, height = 8)
H9_plot
CT2_plot
dev.off()
#Create simplified plots for figures
  #Change gene name for SNHG14 & SNORD115-45
chr_meltH9lgDEL$Var2 <- gsub('ENSG00000224078', 'SNHG14', chr_meltH9lgDEL$Var2)
chr_meltCT2lgDEL$Var2 <- gsub('ENSG00000224078', 'SNHG14', chr_meltCT2lgDEL$Var2)
chr_meltH9lgDEL$Var2 <- gsub('ENSG00000212380', 'SNORD115-45', chr_meltH9lgDEL$Var2)
chr_meltCT2lgDEL$Var2 <- gsub('ENSG00000212380', 'SNORD115-45', chr_meltCT2lgDEL$Var2)
  #Change names of above genes in appropriate filtered file so they end up color coded on graph correctly
chr15qde_H9lgDEL$external_gene_name <- gsub('ENSG00000224078', 'SNHG14', chr15qde_H9lgDEL$external_gene_name)
chr15qde_CT2lgDEL$external_gene_name <- gsub('ENSG00000224078', 'SNHG14', chr15qde_CT2lgDEL$external_gene_name)
chr15qde_H9lgDEL$external_gene_name <- gsub('ENSG00000212380', 'SNORD115-45', chr15qde_H9lgDEL$external_gene_name)
chr15qde_CT2lgDEL$external_gene_name <- gsub('ENSG00000212380', 'SNORD115-45', chr15qde_CT2lgDEL$external_gene_name)
  #Remove genes that don't have gene symbols & LINC's
chr_meltH9lgDEL <- chr_meltH9lgDEL[!grepl("ENSG", chr_meltH9lgDEL$Var2), ]
chr_meltCT2lgDEL <- chr_meltCT2lgDEL[!grepl("ENSG", chr_meltCT2lgDEL$Var2), ]
chr_meltH9lgDEL <- chr_meltH9lgDEL[!grepl("LINC", chr_meltH9lgDEL$Var2), ]
chr_meltCT2lgDEL <- chr_meltCT2lgDEL[!grepl("LINC", chr_meltCT2lgDEL$Var2), ]
  #Use more narrow region
chr_meltH9lgDEL <- chr_meltH9lgDEL[(7:366),] 
chr_meltCT2lgDEL <- chr_meltCT2lgDEL[(7:366),]
  #Make unique list of gene names
GnsList <- unique(chr_meltH9lgDEL$Var2)
  #Set factor level to order x-axis by gene start site
GnsList <- factor(GnsList, levels = GnsList)
  #Plot
H9_plot2 <- ggplot(chr_meltH9lgDEL, aes(x = Var2, y = value)) + 
  geom_boxplot(aes(color = Var2), filter(chr_meltH9lgDEL, Var2 %in% chr15qde_H9lgDEL$external_gene_name), color = "#AB5000", outlier.shape = NA) + geom_boxplot(aes(color = Var2), filter(chr_meltH9lgDEL, Var2 %in% chr15qNONde_H9lgDEL$external_gene_name), color = "black", outlier.shape = NA) +
  geom_point(aes(color = Var2), filter(chr_meltH9lgDEL, Var2 %in% chr15qde_H9lgDEL$external_gene_name), size = 2, color = "#AB5000", alpha = 0.5, position = jitter) + geom_point(aes(color = Var2), filter(chr_meltH9lgDEL, Var2 %in% chr15qNONde_H9lgDEL$external_gene_name), size = 2, color = "black", alpha = 0.5, position = jitter) +
  scale_x_discrete(limits = levels(GnsList)) + ylim(-12.5,5) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6), panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "lightgray"), panel.background = element_rect(fill = "white"), axis.line = element_line(color="black")) +
  labs(title="Chr15q11-q13 Gene Expression: H9lgDEL", x ="Gene", y = "log2(foldChange)")
H9_plot2
CT2_plot2 <- ggplot(chr_meltCT2lgDEL, aes(x = Var2, y = value)) + 
  geom_boxplot(aes(color = Var2), filter(chr_meltCT2lgDEL, Var2 %in% chr15qde_CT2lgDEL$external_gene_name), color = "#5D3A9B", outlier.shape = NA) + geom_boxplot(aes(color = Var2), filter(chr_meltCT2lgDEL, Var2 %in% chr15qNONde_CT2lgDEL$external_gene_name), color = "black", outlier.shape = NA) +
  geom_point(aes(color = Var2), filter(chr_meltCT2lgDEL, Var2 %in% chr15qde_CT2lgDEL$external_gene_name), size = 2, color = "#5D3A9B", alpha = 0.5, position = jitter) + geom_point(aes(color = Var2), filter(chr_meltCT2lgDEL, Var2 %in% chr15qNONde_CT2lgDEL$external_gene_name), size = 2, color = "black", alpha = 0.5, position = jitter) +
  scale_x_discrete(limits = levels(GnsList)) + ylim(-12.5,5) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6), panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "lightgray"), panel.background = element_rect(fill = "white"), axis.line = element_line(color="black")) +
  labs(title="Chr15q11-q13 Gene Expression: CT2lgDEL", x ="Gene", y = "log2(foldChange)")
CT2_plot2
  #Save as pdf
pdf("Chr15q11-q13GeneExp_lgDEL_gencodev25_narrow.pdf", width = 10, height = 8)
H9_plot2
CT2_plot2
dev.off()
#Arrange plots
both_plot <- ggarrange(H9_plot2, CT2_plot2, ncol = 1, nrow = 2)
both_plot
  #Save as pdf
pdf("Chr15q11-q13GeneExp_lgDEL_gencodev25_facet.pdf", width = 10, height = 8)
both_plot
dev.off()

##---Venn Diagram---##
H9lgDELsig <- ResGeneID_H9lgDEL[which(ResGeneID_H9lgDEL$padj < 0.05),]
upH9lgDEL <- H9lgDELsig[which(H9lgDELsig$log2FoldChange > 0),]
downH9lgDEL <- H9lgDELsig[which(H9lgDELsig$log2FoldChange < 0),]
CT2lgDELsig <- ResGeneID_CT2lgDEL[which(ResGeneID_CT2lgDEL$padj < 0.05),]
upCT2lgDEL <- CT2lgDELsig[which(CT2lgDELsig$log2FoldChange > 0),]
downCT2lgDEL <- CT2lgDELsig[which(CT2lgDELsig$log2FoldChange < 0),]
gene_list <- list(A = sample(downH9lgDEL$external_gene_name),
                  B = sample(downCT2lgDEL$external_gene_name),
                  C = sample(upCT2lgDEL$external_gene_name),
                  D = sample(upH9lgDEL$external_gene_name))
p1 <- ggVennDiagram(gene_list, 
                    category.names = c("DownH9lgDEL","DownCT2lgDEL","UpCT2lgDEL","UpH9lgDEL"),
                    label = "count")
p1
#Double check numbers correlate when using ENSEMBL id's versus external_gene_names
gene_list_ENSEMBL <- list(A = sample(as.character(downH9lgDEL$Gene)),
                          B = sample(as.character(downCT2lgDEL$Gene)),
                          C = sample(as.character(upCT2lgDEL$Gene)),
                          D = sample(as.character(upH9lgDEL$Gene)))
p2 <- ggVennDiagram(gene_list_ENSEMBL, 
                   category.names = c("DownH9lgDEL","DownCT2lgDEL","UpCT2lgDEL","UpH9lgDEL"),
                   label = "count")
p2
#Expand axis to show long labels
p1 + scale_x_continuous(expand = expansion(mult = .2))
p2 + scale_x_continuous(expand = expansion(mult = .2))
#Save as pdf
pdf("VennDiagramlgDEL.pdf", width = 10, height = 8)
p1
p2
dev.off()

#Pull out intersection values
intersect_Venn <- process_region_data(Venn(gene_list))
intersect_Venn_ENSEMBL <- process_region_data(Venn(gene_list_ENSEMBL))
#Save as csv
intersect_Venn <- as.data.frame(intersect_Venn)
intersect_Venn <- apply(intersect_Venn,2,as.character)
write.csv(intersect_Venn, file = "VennDiagramlgDEL_IntersectResults.csv")
intersect_Venn_ENSEMBL <- as.data.frame(intersect_Venn_ENSEMBL)
intersect_Venn_ENSEMBL <- apply(intersect_Venn_ENSEMBL,2,as.character)
write.csv(intersect_Venn_ENSEMBL, file = "VennDiagramlgDEL_IntersectResultsENSEMBL.csv")
#Subset & save lists
down_lgDEL <- list(intersect_Venn[5,3])
lapply(down_lgDEL, write, "VennDiagramlgDEL_IntersectDownResults.txt", append=TRUE, ncolumns=1000)
up_lgDEL <- list(intersect_Venn[10,3])
lapply(up_lgDEL, write, "VennDiagramlgDEL_IntersectUpResults.txt", append=TRUE, ncolumns=1000)
down_lgDEL_ENSEMBL <- list(intersect_Venn_ENSEMBL[5,3])
lapply(down_lgDEL_ENSEMBL, write, "VennDiagramlgDEL_IntersectDownResultsENSEMBL.txt", append=TRUE, ncolumns=1000)
up_lgDEL_ENSEMBL <- list(intersect_Venn_ENSEMBL[10,3])
lapply(up_lgDEL_ENSEMBL, write, "VennDiagramlgDEL_IntersectUpResultsENSEMBL.txt", append=TRUE, ncolumns=1000)

##Test for significance of overlaps
#Create gene lists
all_H9_genes <- list(ResGeneID_H9lgDEL$external_gene_name)
write.table(all_H9_genes, "AllH9Genes.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
all_CT2_genes <- list(ResGeneID_CT2lgDEL$external_gene_name)
write.table(all_CT2_genes, "AllCT2Genes.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
#Copy text files onto the cluster and run venn_intersect_test_rbg.sbatch (Done with ENSEMBL numbers)
#Copy files from permutations to computer and read in
A_B <- read.table(file='A_B.txt', header=FALSE)
C_D <- read.table(file='C_D.txt', header=FALSE)
#Calculate median value
median(A_B$V1) #equal to 124
median(C_D$V1) #equal to 94
#Plot histograms
hist_AB <- ggplot(A_B, aes(x=V1)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = 124), col='red', size=.6) +
  geom_vline(aes(xintercept = 373), col='purple', size=1, linetype="dashed") +
  scale_x_continuous(breaks = seq(0, 400, 20)) +
  labs(title="Histogram of A_B Permutation",x="Number of Overlapping Genes", y = "Frequency")
hist_AB
hist_CD <- ggplot(C_D, aes(x=V1)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = 94), col='red', size=.6) +
  geom_vline(aes(xintercept = 465), col='purple', size=1, linetype="dashed") +
  scale_x_continuous(breaks = seq(0, 480, 20)) +
  labs(title="Histogram of C_D Permutation",x="Number of Overlapping Genes", y = "Frequency")
hist_CD
#Save as pdf
pdf("VennOverlapPermutationsHistogramlgDEL.pdf", width = 10, height = 8)
hist_AB
hist_CD
dev.off()

##---Data "Housekeeping" part 2---##
#Take intersect gene lists and use Notepad++ to make them into a list with one gene per row (Find & Replace ", " with \n)
#Read in edited text files as the objects
down_lgDEL <- as.data.frame(read.table("VennDiagramlgDEL_IntersectDownResults.txt", sep="\t", head=F))
names(down_lgDEL)[1] <- "external_gene_name"
up_lgDEL <- as.data.frame(read.table("VennDiagramlgDEL_IntersectUpResults.txt", sep="\t", head=F))
names(up_lgDEL)[1] <- "external_gene_name"
down_lgDEL_ENSEMBL <- as.data.frame(read.table("VennDiagramlgDEL_IntersectDownResultsENSEMBL.txt", sep="\t", head=F))
names(down_lgDEL_ENSEMBL)[1] <- "Gene"
up_lgDEL_ENSEMBL <- as.data.frame(read.table("VennDiagramlgDEL_IntersectUpResultsENSEMBL.txt", sep="\t", head=F))
names(up_lgDEL_ENSEMBL)[1] <- "Gene"
#Combine lists of genes at intersects to use downstream
GeneName_Intersect <- rbind(up_lgDEL, down_lgDEL)
GeneID_Intersect <- rbind(up_lgDEL_ENSEMBL, down_lgDEL_ENSEMBL)
#Order list
GeneID_Intersect <- GeneID_Intersect[order(GeneID_Intersect$Gene),]
#Parse results file by intersect list
ResGeneID_ALL <- ResGeneID_ALL[order(ResGeneID_ALL$Gene),]
ResGene_ALL_sharedID <- ResGeneID_ALL[GeneID_Intersect,]

##---Volcano Plot---##
#Calculate average log2FoldChange
ResGene_ALL_sharedID$log2FoldChange_avg<-rowMeans(ResGene_ALL_sharedID[,c(3,51)])
#Make volcano plots
EnVol_H9 <- EnhancedVolcano(ResGene_ALL_sharedID,
                lab = ResGene_ALL_sharedID$external_gene_name,
                x = 'log2FoldChange_avg',
                FCcutoff = 0,
                y = 'padj_H9',
                pCutoff = 5e-2,
                title = 'H9lgDEL vs H9 WT',
                subtitle = 'DESeq2 Results')
EnVol_H9
EnVol_CT2 <- EnhancedVolcano(ResGene_ALL_sharedID,
                            lab = ResGene_ALL_sharedID$external_gene_name,
                            x = 'log2FoldChange_avg',
                            FCcutoff = 0,
                            y = 'padj_CT2',
                            pCutoff = 5e-2,
                            title = 'CT2lgDEL vs CT2 WT',
                            subtitle = 'DESeq2 Results')
EnVol_CT2
#Save as pdf
pdf("Volc_SharedGeneslgDEL.pdf", width = 10, height = 8)
EnVol_H9
EnVol_CT2
dev.off()

##---Heatmap: All shared DE genes---##
allDEgenes_shared <- ResGene_ALL_sharedID$Gene
rld4heatmapde_shared <- assay(rld)[allDEgenes_shared,]
  #row clustering order
hr_shared <- hclust(as.dist(1-cor(t(rld4heatmapde_shared), method="spearman")), method="complete")
    #column clustering order
hc_shared <- hclust(as.dist(1-cor(rld4heatmapde_shared, method="spearman")), method="complete")
heatmap.2( rld4heatmapde_shared, key=T, scale="row", 
           trace="none", Rowv=as.dendrogram(hr_shared), Colv=as.dendrogram(hc_shared), cexRow = .2, cexCol = .5,
           col = colorRampPalette(c("blue", "black", "yellow"))(n = 1000))
#Make pdf
pdf("AllsharedDEGenesheatmap.pdf", w=11.5, h=8, pointsize=8)
heatmap.2( rld4heatmapde_shared, key=T, scale="row", 
           trace="none", Rowv=as.dendrogram(hr_shared), Colv=as.dendrogram(hc_shared), cexRow = .2, cexCol = .5,
           col = colorRampPalette(c("blue", "black", "yellow"))(n = 1000))
dev.off()
#Break up heatmap based on number of clusters
dend2 <- as.dendrogram(hr_shared)
plot(dend2)
dend2 <- color_branches(hr_shared, k = 5)
dend2 <- color_labels(hr_shared, k = 5)
plot(dend2)
#Make into pdf
pdf("SharedGenesdendro.pdf", w=11.5, h=8)
plot(dend2)
dev.off()
#Take genes for gene ontologies based on clusters
clusterids <- cutree(hr_shared, k = 5)
clusterids[hr_shared$order]
rld4heatmapdecluster <-cbind(rld4heatmapde_shared,clusterID=clusterids)
rld4heatmapdecluster_data <- rld4heatmapdecluster[hr_shared$order,]
write.csv(rld4heatmapdecluster_data, file="Shared_DEGenes_clusterids.csv")

##---Heatmap: Top DE genes---##
resOrdered_shared <- ResGene_ALL_sharedID[order(ResGene_ALL_sharedID$log2FoldChange_avg),]
topGenes_neg <- head(rownames(resOrdered_shared),50)
topGenes_pos <- tail(rownames(resOrdered_shared),50)
topGenes_shared <- rbind(topGenes_pos, topGenes_neg)
mat1 <- assay(rld)[topGenes_shared,]
mat1 <- mat1 - rowMeans(mat1)
gns1 <- as.matrix(ResGeneID_ALL[topGenes_shared,]$external_gene_name)
row.names(gns1) <- rownames(ResGeneID_ALL[topGenes_shared,])
row.names(mat1)[match(row.names(gns1), row.names(mat1))] <- gns1
df <- as.data.frame(subset(colData(dds), select = c("Condition", "Genetic.Background", "sizeFactor")))
ann_colors <- list(
  sizeFactor = c("white", "firebrick"),
  Genetic.Background = c(CT2 = "#7570B3", H9 = "#E7298A"),
  Condition = c(WT = "#419644", smDEL = "#64B4F5", lgDEL = "#2A6697"))
pheatmap(mat1, annotation_col = df,
         main = "Shared WT vs lgDEL Top100Genes",
         angle_col = 45, fontsize_row = 3, fontsize_col = 5,
         annotation_colors = ann_colors)
#Save plot using export function

#Save dds object
save(dds, file = "PWSlgDEL_DESeqHISAT2.RData")

#Save summary tables
write.csv(as.data.frame(ResGeneID_H9lgDEL), file = "H9lgDEL_Results_DESeqHISAT2.csv")
write.csv(as.data.frame(ResGeneID_CT2lgDEL), file = "CT2lgDEL_Results_DESeqHISAT2.csv")

#Combine saved plots using Adobe
