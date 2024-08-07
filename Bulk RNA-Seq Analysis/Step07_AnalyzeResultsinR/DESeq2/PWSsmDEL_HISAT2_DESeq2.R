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
library("reshape2")
library("ggpubr")

#Set working directory
directory <- "../Cotney_Lab/PWS_RNASeq/HISAT2/smDEL"
setwd(directory)

#Import file
counts <- read.delim("PWSexp_featureCounts.txt", sep="\t", head=T, skip = 1, row.names = "Geneid")

#Simplify column names
colnames(counts) <- c("chr", "start", "end", "strand", "length", "CT2_1", "CT2_2", "CT2_3", "CT2_4", "CT2_5", "CT2_6", "CT2_lgDEL_1", "CT2_lgDEL_2", "CT2_lgDEL_3", "CT2_lgDEL_4", "CT2_lgDEL_5", "CT2_lgDEL_6", "CT2_smDEL_1", "CT2_smDEL_2", "CT2_smDEL_3", "CT2_smDEL_4", "CT2_smDEL_5", "CT2_smDEL_6", "H9_1", "H9_2", "H9_3", "H9_4", "H9_5", "H9_6", "H9_lgDEL_1", "H9_lgDEL_2", "H9_lgDEL_5", "H9_smDEL_1", "H9_smDEL_2", "H9_smDEL_3", "H9_smDEL_4", "H9_smDEL_5")
#Get rid of first 5 columns
counts <- counts[ ,6:ncol(counts)]
#Turn into matrix
counts <- as.matrix(counts)

#Load in metadata
colData <- as.data.frame(read.table("Metadata.txt", sep="\t", head=T, row.names = "ID"))
colData$Condition <- as.factor(colData$Condition)
colData$Background_Condition <- as.factor(paste(colData$Genetic.Background, colData$Condition, sep="_"))

#Create DESeq objects
ddstotal <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ Background_Condition)
ddstotal

#Trim matrix, remove low counts
summary(rowSums(counts(ddstotal)))
#Filter using cutoff (can modify if necessary)
dds <- ddstotal [rowSums(counts(ddstotal)) > 1 ,]

#Run DESeq on data
dds <- DESeq(dds)
  #To load
  load("PWS_DESeqHISAT2.RData")
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
pdf("PWS_Exp_PCAs.pdf", width = 10, height = 8)
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
resH9smDEL <- results(dds, contrast = c("Background_Condition","H9_smDEL","H9_WT"))
resH9smDEL
resCT2smDEL <- results(dds, contrast = c("Background_Condition","CT2_smDEL","CT2_WT"))
resCT2smDEL
#Merge results table & counts table
resDataH9smDEL <- merge(as.data.frame(resH9smDEL), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resDataCT2smDEL <- merge(as.data.frame(resCT2smDEL), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resDataH9smDEL)[1] <- "Gene"
names(resDataCT2smDEL)[1] <- "Gene"
rownames(resDataH9smDEL) <- resDataH9smDEL$Gene
rownames(resDataCT2smDEL) <- resDataCT2smDEL$Gene
#Change ESMBL ID's to gene names
mart <- useMart("ensembl", host = "https://useast.ensembl.org", dataset = "hsapiens_gene_ensembl")
  #For H9smDEL
R1 <- data.table(resDataH9smDEL)
rownames(R1) <- R1$Gene
genes.table_H9smDEL <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","description","gene_biotype","chromosome_name","start_position","end_position","strand"), values= R1$Gene, mart= mart)
  #For CT2smDEL
R2 <- data.table(resDataCT2smDEL)
rownames(R2) <- R2$Gene
genes.table_CT2smDEL <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","description","gene_biotype","chromosome_name","start_position","end_position","strand"), values= R2$Gene, mart= mart)
#Order lists by ENSEMBL ID so they merge correctly
resOrderedgene_H9smDEL <- R1[order(R1$Gene),]
names(genes.table_H9smDEL)[1] <- "Gene"
resOrderedgene_CT2smDEL <- R2[order(R2$Gene),]
names(genes.table_CT2smDEL)[1] <- "Gene"
#Merge all things
ResGeneID_H9smDEL <- merge(as.data.frame(resOrderedgene_H9smDEL), as.data.frame(genes.table_H9smDEL), by="Gene", sort=TRUE, all = FALSE)
rownames(ResGeneID_H9smDEL) <- ResGeneID_H9smDEL$Gene
ResGeneID_CT2smDEL <- merge(as.data.frame(resOrderedgene_CT2smDEL), as.data.frame(genes.table_CT2smDEL), by="Gene", sort=TRUE, all = FALSE)
rownames(ResGeneID_CT2smDEL) <- ResGeneID_CT2smDEL$Gene
#Fill in blank external gene names
  #For H9smDEL
ResGeneID_H9smDEL$external_gene_name <- as.character(ResGeneID_H9smDEL$external_gene_name)
ResGeneID_H9smDEL$external_gene_name[ResGeneID_H9smDEL$external_gene_name==""] <- NA
ResGeneID_H9smDEL$external_gene_name <- as.factor(ResGeneID_H9smDEL$external_gene_name)
GeneSymbolDE_H9smDEL <- mapIds(org.Hs.eg.db, keys = ResGeneID_H9smDEL$Gene, column = "SYMBOL", keytype = "ENSEMBL", multivals= "first")
for (i in c(1:length(GeneSymbolDE_H9smDEL))){
  
  if (is.na(GeneSymbolDE_H9smDEL[i])){
    
    GeneSymbolDE_H9smDEL[i]<-gsub("\\..*","",ResGeneID_H9smDEL$Gene[i])
    
  }
  
}
ResGeneID_H9smDEL$external_gene_name <- GeneSymbolDE_H9smDEL
  #For CT2smDEL
ResGeneID_CT2smDEL$external_gene_name <- as.character(ResGeneID_CT2smDEL$external_gene_name)
ResGeneID_CT2smDEL$external_gene_name[ResGeneID_CT2smDEL$external_gene_name==""] <- NA
ResGeneID_CT2smDEL$external_gene_name <- as.factor(ResGeneID_CT2smDEL$external_gene_name)
GeneSymbolDE_CT2smDEL <- mapIds(org.Hs.eg.db, keys = ResGeneID_CT2smDEL$Gene, column = "SYMBOL", keytype = "ENSEMBL", multivals= "first")
for (i in c(1:length(GeneSymbolDE_CT2smDEL))){
  
  if (is.na(GeneSymbolDE_CT2smDEL[i])){
    
    GeneSymbolDE_CT2smDEL[i]<-gsub("\\..*","",ResGeneID_CT2smDEL$Gene[i])
    
  }
  
}
ResGeneID_CT2smDEL$external_gene_name <- GeneSymbolDE_CT2smDEL
#Merge both results objects
ResGeneID_ALL <- merge(as.data.frame(ResGeneID_H9smDEL), as.data.frame(ResGeneID_CT2smDEL), by="Gene", sort=TRUE, all = FALSE)
rownames(ResGeneID_ALL) <- ResGeneID_ALL$Gene
#Drop unnecessary columns
ResGeneID_ALL <- ResGeneID_ALL[-c(53:91)]
#Simplify column names
colnames(ResGeneID_ALL) <- c("Gene", "baseMean_H9", "log2FoldChange_H9", "lfcSE_H9", "stat_H9", "pvalue_H9", "padj_H9", "CT2_1", "CT2_2", "CT2_3", "CT2_4", "CT2_5", "CT2_6", "CT2_lgDEL_1", "CT2_lgDEL_2", "CT2_lgDEL_3", "CT2_lgDEL_4", "CT2_lgDEL_5", "CT2_lgDEL_6", "CT2_smDEL_1", "CT2_smDEL_2", "CT2_smDEL_3", "CT2_smDEL_4", "CT2_smDEL_5", "CT2_smDEL_6", "H9_1", "H9_2", "H9_3", "H9_4", "H9_5", "H9_6", "H9_lgDEL_1", "H9_lgDEL_2", "H9_lgDEL_5", "H9_smDEL_1", "H9_smDEL_2", "H9_smDEL_3", "H9_smDEL_4", "H9_smDEL_5", "external_gene_name", "description", "gene_biotype", "chromosome_name", "start_position", "end_position", "strand", "baseMean_CT2", "log2FoldChange_CT2", "lfcSE_CT2", "stat_CT2", "pvalue_CT2", "padj_CT2")
#Save summary table
write.csv(as.data.frame(ResGeneID_ALL), file = "CombinedH9CT2_Results_DESeqHISAT2.csv")
  #To read in
  ResGeneID_ALL <- read.csv(file = 'CombinedH9CT2_Results_DESeqHISAT2.csv', header = TRUE, row.names = "X")

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
chr15qOrder_plus1[,c(8:39)]<- data.frame(lapply(chr15qOrder_plus1[,c(8:39)],plusone))
#Calculate mean of WT values
chr15qOrder_plus1$wtmed_H9 <- rowMeans(chr15qOrder_plus1[,c(26:31)])
chr15qOrder_plus1$wtmed_CT2 <- rowMeans(chr15qOrder_plus1[,c(8:13)])
#Divide large deletion samples by WT mean
chr15qOrder_plus1[,c(35:39)] <- chr15qOrder_plus1[,c(35:39)]/chr15qOrder_plus1$wtmed_H9
chr15qOrder_plus1[,c(20:25)] <- chr15qOrder_plus1[,c(20:25)]/chr15qOrder_plus1$wtmed_CT2
#Log transform (same columns as division step above)
chr15qOrder_plus1_log2 <- chr15qOrder_plus1
chr15qOrder_plus1_log2[,c(35:39)] <- log2(chr15qOrder_plus1_log2[,c(35:39)])
chr15qOrder_plus1_log2[,c(20:25)] <- log2(chr15qOrder_plus1_log2[,c(20:25)])
#Transpose
chr15qOrderH9smDEL_plus1_log2_t <- t(chr15qOrder_plus1_log2[,c(35:39)])
chr15qOrderCT2smDEL_plus1_log2_t <- t(chr15qOrder_plus1_log2[,c(20:25)])
#Make the gene name the title of the columns
colnames(chr15qOrderH9smDEL_plus1_log2_t) <- chr15qOrder_plus1_log2$external_gene_name
colnames(chr15qOrderCT2smDEL_plus1_log2_t) <- chr15qOrder_plus1_log2$external_gene_name
  #Save final 15q gene expression table
  write.csv(as.data.frame(chr15qOrder_plus1_log2), file = "CombinedsmDEL_DESeqHISAT2_15qgenes_gencodev25.csv")
#Filter DE vs nonDE genes for chr15
chr15qde_H9smDEL <- chr15qdat[which(chr15qdat$padj_H9 < 0.05),]
chr15qde_CT2smDEL <- chr15qdat[which(chr15qdat$padj_CT2 < 0.05),]
chr15qNONde_H9smDEL <- chr15qdat[!(row.names(chr15qdat) %in% row.names(chr15qde_H9smDEL)),]
chr15qNONde_CT2smDEL <- chr15qdat[!(row.names(chr15qdat) %in% row.names(chr15qde_CT2smDEL)),]
#Plot with ggplot
  #Reshape data
chr_meltH9smDEL <- reshape2::melt(chr15qOrderH9smDEL_plus1_log2_t)
head(chr_meltH9smDEL)
chr_meltCT2smDEL <- reshape2::melt(chr15qOrderCT2smDEL_plus1_log2_t)
head(chr_meltCT2smDEL)
  #Set factors
jitter <- position_jitter(width = 0.18, height = 0.1)
chr15qOrder_plus1_log2$external_gene_name <- factor(chr15qOrder_plus1_log2$external_gene_name, levels = chr15qOrder_plus1_log2$external_gene_name)
  #Plot
H9_plot <- ggplot(chr_meltH9smDEL, aes(x = Var2, y = value)) + 
  geom_boxplot(aes(color = Var2), filter(chr_meltH9smDEL, Var2 %in% chr15qde_H9smDEL$external_gene_name), color = "#E66100", outlier.shape = NA) + geom_boxplot(aes(color = Var2), filter(chr_meltH9smDEL, Var2 %in% chr15qNONde_H9smDEL$external_gene_name), color = "black", outlier.shape = NA) +
  geom_point(aes(color = Var2), filter(chr_meltH9smDEL, Var2 %in% chr15qde_H9smDEL$external_gene_name), size = 1, color = "#E66100", alpha = 0.5, position = jitter) + geom_point(aes(color = Var2), filter(chr_meltH9smDEL, Var2 %in% chr15qNONde_H9smDEL$external_gene_name), size = 1, color = "black", alpha = 0.5, position = jitter) +
  scale_x_discrete(limits = levels(chr15qOrder_plus1_log2$external_gene_name)) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6), panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "lightgray"), panel.background = element_rect(fill = "white"), axis.line = element_line(color="black")) +
  labs(title="Chr15q11-q13 Gene Expression: H9smDEL", x ="Gene", y = "log2(foldChange)")
H9_plot
CT2_plot <- ggplot(chr_meltCT2smDEL, aes(x = Var2, y = value)) + 
  geom_boxplot(aes(color = Var2), filter(chr_meltCT2smDEL, Var2 %in% chr15qde_CT2smDEL$external_gene_name), color = "#9356FF", outlier.shape = NA) + geom_boxplot(aes(color = Var2), filter(chr_meltCT2smDEL, Var2 %in% chr15qNONde_CT2smDEL$external_gene_name), color = "black", outlier.shape = NA) +
  geom_point(aes(color = Var2), filter(chr_meltCT2smDEL, Var2 %in% chr15qde_CT2smDEL$external_gene_name), size = 1, color = "#9356FF", alpha = 0.5, position = jitter) + geom_point(aes(color = Var2), filter(chr_meltCT2smDEL, Var2 %in% chr15qNONde_CT2smDEL$external_gene_name), size = 1, color = "black", alpha = 0.5, position = jitter) +
  scale_x_discrete(limits = levels(chr15qOrder_plus1_log2$external_gene_name)) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6), panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "lightgray"), panel.background = element_rect(fill = "white"), axis.line = element_line(color="black")) +
  labs(title="Chr15q11-q13 Gene Expression: CT2smDEL", x ="Gene", y = "log2(foldChange)")
CT2_plot
#Save as pdf
pdf("Chr15q11-q13GeneExp_smDEL_gencodev25.pdf", width = 10, height = 8)
H9_plot
CT2_plot
dev.off()
#Create simplified plots for figures
  #Change gene name for SNHG14 & SNORD115-45
chr_meltH9smDEL$Var2 <- gsub('ENSG00000224078', 'SNHG14', chr_meltH9smDEL$Var2)
chr_meltCT2smDEL$Var2 <- gsub('ENSG00000224078', 'SNHG14', chr_meltCT2smDEL$Var2)
chr_meltH9smDEL$Var2 <- gsub('ENSG00000212380', 'SNORD115-45', chr_meltH9smDEL$Var2)
chr_meltCT2smDEL$Var2 <- gsub('ENSG00000212380', 'SNORD115-45', chr_meltCT2smDEL$Var2)
  #Change names of above genes in appropriate filtered file so they end up color coded on graph correctly
chr15qde_H9smDEL$external_gene_name <- gsub('ENSG00000224078', 'SNHG14', chr15qde_H9smDEL$external_gene_name)
chr15qde_CT2smDEL$external_gene_name <- gsub('ENSG00000224078', 'SNHG14', chr15qde_CT2smDEL$external_gene_name)
chr15qNONde_H9smDEL$external_gene_name <- gsub('ENSG00000212380', 'SNORD115-45', chr15qNONde_H9smDEL$external_gene_name)
chr15qNONde_CT2smDEL$external_gene_name <- gsub('ENSG00000212380', 'SNORD115-45', chr15qNONde_CT2smDEL$external_gene_name)
  #Remove genes that don't have gene symbols & LINC's
chr_meltH9smDEL <- chr_meltH9smDEL[!grepl("ENSG", chr_meltH9smDEL$Var2), ]
chr_meltCT2smDEL <- chr_meltCT2smDEL[!grepl("ENSG", chr_meltCT2smDEL$Var2), ]
chr_meltH9smDEL <- chr_meltH9smDEL[!grepl("LINC", chr_meltH9smDEL$Var2), ]
chr_meltCT2smDEL <- chr_meltCT2smDEL[!grepl("LINC", chr_meltCT2smDEL$Var2), ]
  #Use more narrow region
chr_meltH9smDEL <- chr_meltH9smDEL[(7:305),] 
chr_meltCT2smDEL <- chr_meltCT2smDEL[(7:366),]
  #Make unique list of gene names
GnsList <- unique(chr_meltH9smDEL$Var2)
  #Set factor level to order x-axis by gene start site
GnsList <- factor(GnsList, levels = GnsList)
  #Plot
H9_plot2 <- ggplot(chr_meltH9smDEL, aes(x = Var2, y = value)) + 
  geom_boxplot(aes(color = Var2), filter(chr_meltH9smDEL, Var2 %in% chr15qde_H9smDEL$external_gene_name), color = "#E66100", outlier.shape = NA) + geom_boxplot(aes(color = Var2), filter(chr_meltH9smDEL, Var2 %in% chr15qNONde_H9smDEL$external_gene_name), color = "black", outlier.shape = NA) +
  geom_point(aes(color = Var2), filter(chr_meltH9smDEL, Var2 %in% chr15qde_H9smDEL$external_gene_name), size = 2, color = "#E66100", alpha = 0.5, position = jitter) + geom_point(aes(color = Var2), filter(chr_meltH9smDEL, Var2 %in% chr15qNONde_H9smDEL$external_gene_name), size = 2, color = "black", alpha = 0.5, position = jitter) +
  scale_x_discrete(limits = levels(GnsList)) + scale_y_continuous(limits = c(-4.1,4.1), breaks = c(-4,-2,0,2,4)) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6), panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "lightgray"), panel.background = element_rect(fill = "white"), axis.line = element_line(color="black")) +
  labs(title="Chr15q11-q13 Gene Expression: H9smDEL", x ="Gene", y = "log2(foldChange)")
H9_plot2
CT2_plot2 <- ggplot(chr_meltCT2smDEL, aes(x = Var2, y = value)) + 
  geom_boxplot(aes(color = Var2), filter(chr_meltCT2smDEL, Var2 %in% chr15qde_CT2smDEL$external_gene_name), color = "#9356FF", outlier.shape = NA) + geom_boxplot(aes(color = Var2), filter(chr_meltCT2smDEL, Var2 %in% chr15qNONde_CT2smDEL$external_gene_name), color = "black", outlier.shape = NA) +
  geom_point(aes(color = Var2), filter(chr_meltCT2smDEL, Var2 %in% chr15qde_CT2smDEL$external_gene_name), size = 2, color = "#9356FF", alpha = 0.5, position = jitter) + geom_point(aes(color = Var2), filter(chr_meltCT2smDEL, Var2 %in% chr15qNONde_CT2smDEL$external_gene_name), size = 2, color = "black", alpha = 0.5, position = jitter) +
  scale_x_discrete(limits = levels(GnsList)) + scale_y_continuous(limits = c(-4.1,4.1), breaks = c(-4,-2,0,2,4)) +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=6), panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "lightgray"), panel.background = element_rect(fill = "white"), axis.line = element_line(color="black")) +
  labs(title="Chr15q11-q13 Gene Expression: CT2smDEL", x ="Gene", y = "log2(foldChange)")
CT2_plot2
  #Save as pdf
pdf("Chr15q11-q13GeneExp_smDEL_gencodev25_narrow.pdf", width = 10, height = 8)
H9_plot2
CT2_plot2
dev.off()
#Arrange plots
both_plot <- ggarrange(H9_plot2, CT2_plot2, ncol = 1, nrow = 2)
both_plot
  #Save as pdf
pdf("Chr15q11-q13GeneExp_smDEL_gencodev25_facet.pdf", width = 10, height = 8)
both_plot
dev.off()

##---Venn Diagram---##
H9smDELsig <- ResGeneID_H9smDEL[which(ResGeneID_H9smDEL$padj < 0.05),]
upH9smDEL <- H9smDELsig[which(H9smDELsig$log2FoldChange > 0),]
downH9smDEL <- H9smDELsig[which(H9smDELsig$log2FoldChange < 0),]
CT2smDELsig <- ResGeneID_CT2smDEL[which(ResGeneID_CT2smDEL$padj < 0.05),]
upCT2smDEL <- CT2smDELsig[which(CT2smDELsig$log2FoldChange > 0),]
downCT2smDEL <- CT2smDELsig[which(CT2smDELsig$log2FoldChange < 0),]
gene_list <- list(A = sample(downH9smDEL$external_gene_name),
                  B = sample(downCT2smDEL$external_gene_name),
                  C = sample(upCT2smDEL$external_gene_name),
                  D = sample(upH9smDEL$external_gene_name))
p1 <- ggVennDiagram(gene_list, 
                    category.names = c("DownH9smDEL","DownCT2smDEL","UpCT2smDEL","UpH9smDEL"),
                    label = "count")
p1
#Double check numbers correlate when using ENSEMBL id's versus external_gene_names
gene_list_ENSEMBL <- list(A = sample(as.character(downH9smDEL$Gene)),
                          B = sample(as.character(downCT2smDEL$Gene)),
                          C = sample(as.character(upCT2smDEL$Gene)),
                          D = sample(as.character(upH9smDEL$Gene)))
p2 <- ggVennDiagram(gene_list_ENSEMBL, 
                   category.names = c("DownH9smDEL","DownCT2smDEL","UpCT2smDEL","UpH9smDEL"),
                   label = "count")
p2
#Expand axis to show long labels
p1 + scale_x_continuous(expand = expansion(mult = .2))
#Save as pdf
pdf("VennDiagram.pdf", width = 10, height = 8)
p1
dev.off()
#Pull out intersection values
intersect_Venn <- process_region_data(Venn(gene_list))
intersect_Venn_ENSEMBL <- process_region_data(Venn(gene_list_ENSEMBL))
#Save as csv
intersect_Venn <- as.data.frame(intersect_Venn)
intersect_Venn <- apply(intersect_Venn,2,as.character)
write.csv(intersect_Venn, file = "VennDiagram_IntersectResults.csv")
intersect_Venn_ENSEMBL <- as.data.frame(intersect_Venn_ENSEMBL)
intersect_Venn_ENSEMBL <- apply(intersect_Venn_ENSEMBL,2,as.character)
write.csv(intersect_Venn_ENSEMBL, file = "VennDiagram_IntersectResultsENSEMBL.csv")
#Subset & save lists
down_smDEL <- list(intersect_Venn[5,3])
lapply(down_smDEL, write, "VennDiagram_IntersectDownResults.txt", append=TRUE, ncolumns=1000)
up_smDEL <- list(intersect_Venn[10,3])
lapply(up_smDEL, write, "VennDiagram_IntersectUpResults.txt", append=TRUE, ncolumns=1000)
down_smDEL_ENSEMBL <- list(intersect_Venn_ENSEMBL[5,3])
lapply(down_smDEL_ENSEMBL, write, "VennDiagram_IntersectDownResultsENSEMBL.txt", append=TRUE, ncolumns=1000)
up_smDEL_ENSEMBL <- list(intersect_Venn_ENSEMBL[10,3])
lapply(up_smDEL_ENSEMBL, write, "VennDiagram_IntersectUpResultsENSEMBL.txt", append=TRUE, ncolumns=1000)
#Create gene lists to run test for significance of overlaps
all_H9_genes <- list(ResGeneID_H9smDEL$external_gene_name)
write.table(all_H9_genes, "AllH9Genes.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
all_CT2_genes <- list(ResGeneID_CT2smDEL$external_gene_name)
write.table(all_CT2_genes, "AllCT2Genes.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
#Read in files from permutations
A_B <- read.table(file='A_B.txt', header=FALSE)
C_D <- read.table(file='C_D.txt', header=FALSE)
#Calculate median value
median(A_B$V1) #equal to 33
median(C_D$V1) #equal to 50
#Plot histograms
hist_AB <- ggplot(A_B, aes(x=V1)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = 33), col='red', size=.6) +
  geom_vline(aes(xintercept = 135), col='purple', size=1, linetype="dashed") +
  scale_x_continuous(breaks = seq(0, 150, 10)) +
  labs(title="Histogram of A_B Permutation",x="Number of Overlapping Genes", y = "Frequency")
hist_AB
hist_CD <- ggplot(C_D, aes(x=V1)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = 50), col='red', size=.6) +
  geom_vline(aes(xintercept = 183), col='purple', size=1, linetype="dashed") +
  scale_x_continuous(breaks = seq(0, 185, 10)) +
  labs(title="Histogram of C_D Permutation",x="Number of Overlapping Genes", y = "Frequency")
hist_CD
#Save as pdf
pdf("VennOverlapPermutationsHistogram.pdf", width = 10, height = 8)
hist_AB
hist_CD
dev.off()

##---Data "Housekeeping" part 2---##
#Take intersect gene lists and use Notepad++ to make them into a list with one gene per row (Find & Replace ", " with \n)
#Read in edited text files as the objects
down_smDEL <- as.data.frame(read.table("VennDiagram_IntersectDownResults.txt", sep="\t", head=F))
names(down_smDEL)[1] <- "external_gene_name"
up_smDEL <- as.data.frame(read.table("VennDiagram_IntersectUpResults.txt", sep="\t", head=F))
names(up_smDEL)[1] <- "external_gene_name"
down_smDEL_ENSEMBL <- as.data.frame(read.table("VennDiagram_IntersectDownResultsENSEMBL.txt", sep="\t", head=F))
names(down_smDEL_ENSEMBL)[1] <- "Gene"
up_smDEL_ENSEMBL <- as.data.frame(read.table("VennDiagram_IntersectUpResultsENSEMBL.txt", sep="\t", head=F))
names(up_smDEL_ENSEMBL)[1] <- "Gene"
#Combine lists of genes at intersects to use downstream
GeneName_Intersect <- rbind(up_smDEL, down_smDEL)
GeneID_Intersect <- rbind(up_smDEL_ENSEMBL, down_smDEL_ENSEMBL)
#Order list
GeneID_Intersect <- GeneID_Intersect[order(GeneID_Intersect$Gene),]
#Parse results file by intersect list
ResGeneID_ALL <- ResGeneID_ALL[order(ResGeneID_ALL$Gene),]
ResGene_ALL_sharedID <- ResGeneID_ALL[GeneID_Intersect,]

##---Volcano Plot---##
#Calculate average log2FoldChange
ResGene_ALL_sharedID$log2FoldChange_avg<-rowMeans(ResGene_ALL_sharedID[,c(3,48)])
#Make volcano plots
EnVol_H9 <- EnhancedVolcano(ResGene_ALL_sharedID,
                lab = ResGene_ALL_sharedID$external_gene_name,
                x = 'log2FoldChange_avg',
                FCcutoff = 0,
                y = 'padj_H9',
                pCutoff = 5e-2,
                title = 'H9smDEL vs H9 WT',
                subtitle = 'DESeq2 Results')
EnVol_H9
EnVol_CT2 <- EnhancedVolcano(ResGene_ALL_sharedID,
                            lab = ResGene_ALL_sharedID$external_gene_name,
                            x = 'log2FoldChange_avg',
                            FCcutoff = 0,
                            y = 'padj_CT2',
                            pCutoff = 5e-2,
                            title = 'CT2smDEL vs CT2 WT',
                            subtitle = 'DESeq2 Results')
EnVol_CT2
#Save as pdf
pdf("Volc_SharedGenes.pdf", width = 10, height = 8)
EnVol_H9
EnVol_CT2
dev.off()

##---Heatmaps---##
#All shared DE genes
allDEgenes_shared <- ResGene_ALL_sharedID$Gene
rld4heatmapde_shared <- assay(rld)[allDEgenes_shared,]
    #Drop lgDEL (because it's wonky)
rld4heatmapde_shared <- as.data.frame(rld4heatmapde_shared)
rld4heatmapde_shared <- rld4heatmapde_shared[-c(7:12)]
rld4heatmapde_shared <- rld4heatmapde_shared[-c(19:21)]
rld4heatmapde_shared <- as.matrix(rld4heatmapde_shared)  
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
dend2 <- color_branches(hr_shared, k = 4)
dend2 <- color_labels(hr_shared, k = 4)
plot(dend2)
#Make into pdf
pdf("SharedGenesdendro.pdf", w=11.5, h=8)
plot(dend2)
dev.off()
#Take genes for gene ontologies based on clusters
clusterids <- cutree(hr_shared, k = 4)
clusterids[hr_shared$order]
rld4heatmapdecluster <-cbind(rld4heatmapde_shared,clusterID=clusterids)
rld4heatmapdecluster_data <- rld4heatmapdecluster[hr_shared$order,]
write.csv(rld4heatmapdecluster_data, file="Shared_DEGenes_clusterids.csv")

#Top DE genes
resOrdered_shared <- ResGene_ALL_sharedID[order(ResGene_ALL_sharedID$log2FoldChange_avg),]
topGenes_neg <- head(rownames(resOrdered_shared),50)
topGenes_pos <- tail(rownames(resOrdered_shared),50)
topGenes_shared <- rbind(topGenes_pos, topGenes_neg)
mat1 <- assay(rld)[topGenes_shared,]
mat1 <- as.data.frame(mat1)
mat1 <- mat1[-c(7:12)]
mat1 <- mat1[-c(19:21)]
mat1 <- as.matrix(mat1)
mat1 <- mat1 - rowMeans(mat1)
gns1 <- as.matrix(ResGeneID_ALL[topGenes_shared,]$external_gene_name)
row.names(gns1) <- rownames(ResGeneID_ALL[topGenes_shared,])
row.names(mat1)[match(row.names(gns1), row.names(mat1))] <- gns1
df <- as.data.frame(subset(colData(dds), select = c("Condition", "Genetic.Background", "sizeFactor")))
ann_colors <- list(
  sizeFactor = c("white", "firebrick"),
  Genetic.Background = c(CT2 = "#1B9E77", H9 = "#D95F02"),
  Condition = c(WT = "#7570B3", smDEL = "#E7298A"))
pheatmap(mat1, annotation_col = df,
         main = "Shared WT vs smDEL Top100Genes",
         angle_col = 45, fontsize_row = 3, fontsize_col = 5,
         annotation_colors = ann_colors)
#Save plot using export function

#Save dds object
save(dds, file = "PWS_DESeqHISAT2.RData")

#Save summary tables
write.csv(as.data.frame(ResGeneID_H9smDEL), file = "H9smDEL_Results_DESeqHISAT2.csv")
write.csv(as.data.frame(ResGeneID_CT2smDEL), file = "CT2smDEL_Results_DESeqHISAT2.csv")

##---Boxplots---##
#Create boxplots of gene expression for 15q11-q13 region
#Import csv made via UCSC genome browser
chr15q <- read.csv("chr15q11-q13.csv")
#Set column names
names(chr15q)[1] <- "transcript_id"
names(chr15q)[2] <- "external_gene_name"
#Change ENSEMBL transcript id's to ENSEMBL gene id's
mart <- useMart("ensembl", host = "https://useast.ensembl.org", dataset = "hsapiens_gene_ensembl")
chr15qGns <- getBM(filters = "ensembl_transcript_id_version", 
                   attributes= c("ensembl_transcript_id_version","ensembl_gene_id", "hgnc_symbol", "description","gene_biotype","chromosome_name","start_position","end_position","strand"), 
                   values= chr15q$transcript_id,
                   mart= mart)
#Convert id's
GeneSymbol<-mapIds(org.Hs.eg.db, keys = chr15qGns$ensembl_gene_id, column = "SYMBOL", keytype = "ENSEMBL", multivals= "first") 
for (i in c(1:length(GeneSymbol))){
  
  if (is.na(GeneSymbol[i])){
    
    GeneSymbol[i]<-gsub("\\..*","",chr15qGns$ensembl_gene_id[i])
    
  }
  
}
chr15qGns$external_gene_name <- GeneSymbol
#Make list of unique ENSEMBL id's
GnsUniq <- unique(chr15qGns$ensembl_gene_id)
#Compare from DE analysis
  #For H9smDEL
    #Filter DE genes for chr15
chr15qde_H9smDEL <- subset(ResGeneID_H9smDEL[(ResGeneID_H9smDEL$Gene) %in% GnsUniq,], padj < 0.05)
    #Parse results file by desired gene list
chr15qdat_H9smDEL <- ResGeneID_H9smDEL[(ResGeneID_H9smDEL$Gene) %in% GnsUniq, ]
    #Put genes in correct order
chr15qOrder_H9smDEL <- chr15qdat_H9smDEL[order(chr15qdat_H9smDEL$start_position),]
    #Add 1 to all counts --> Change columns based on samples!
plusone <- function(x){
  return (x + 1)
}
chr15qOrder_H9smDEL_plus1 <- chr15qOrder_H9smDEL
chr15qOrder_H9smDEL_plus1[,c(8:39)]<- data.frame(lapply(chr15qOrder_H9smDEL[,c(8:39)],plusone))
    #Calculate mean of WT values
chr15qOrder_H9smDEL_plus1$wtmed<-rowMeans(chr15qOrder_H9smDEL_plus1[,c(26:31)])
    #Divide deletion samples by WT mean
chr15qOrder_H9smDEL_plus1[,c(35:39)]<-chr15qOrder_H9smDEL_plus1[,c(35:39)]/chr15qOrder_H9smDEL_plus1$wtmed
    #Log transform (same columns as division step above)
chr15qOrderH9smDEL_plus1_log2 <- chr15qOrder_H9smDEL_plus1
chr15qOrderH9smDEL_plus1_log2[,c(35:39)]<-log2(chr15qOrder_H9smDEL_plus1[,c(35:39)])
    #Transpose
chr15qOrderH9smDEL_plus1_log2_t<-t(chr15qOrderH9smDEL_plus1_log2[,c(35:39)])
    #Make the gene name the title of the columns
colnames(chr15qOrderH9smDEL_plus1_log2_t)<-chr15qOrderH9smDEL_plus1_log2$external_gene_name
    #Save final 15q gene expression table
write.csv(as.data.frame(chr15qOrderH9smDEL_plus1_log2), file = "H9smDEL_DESeqHISAT2_15qgenes.csv")
#Compare from DE analysis
  #For CT2smDEL
    #Filter DE genes for chr15
chr15qde_CT2smDEL <- subset(ResGeneID_CT2smDEL[(ResGeneID_CT2smDEL$Gene) %in% GnsUniq,], padj < 0.05)
    #Parse results file by desired gene list
chr15qdat_CT2smDEL <- ResGeneID_CT2smDEL[(ResGeneID_CT2smDEL$Gene) %in% GnsUniq, ]
    #Put genes in correct order
chr15qOrder_CT2smDEL <- chr15qdat_CT2smDEL[order(chr15qdat_CT2smDEL$start_position),]
    #Add 1 to all counts --> Change columns based on samples!
plusone <- function(x){
  return (x + 1)
}
chr15qOrder_CT2smDEL_plus1 <- chr15qOrder_CT2smDEL
chr15qOrder_CT2smDEL_plus1[,c(8:39)]<- data.frame(lapply(chr15qOrder_CT2smDEL[,c(8:39)],plusone))
    #Calculate mean of WT values
chr15qOrder_CT2smDEL_plus1$wtmed<-rowMeans(chr15qOrder_CT2smDEL_plus1[,c(8:13)])
    #Divide deletion samples by WT mean
chr15qOrder_CT2smDEL_plus1[,c(20:25)]<-chr15qOrder_CT2smDEL_plus1[,c(20:25)]/chr15qOrder_CT2smDEL_plus1$wtmed
    #Log transform (same columns as division step above)
chr15qOrderCT2smDEL_plus1_log2 <- chr15qOrder_CT2smDEL_plus1
chr15qOrderCT2smDEL_plus1_log2[,c(20:25)]<-log2(chr15qOrder_CT2smDEL_plus1[,c(20:25)])
    #Transpose
chr15qOrderCT2smDEL_plus1_log2_t<-t(chr15qOrderCT2smDEL_plus1_log2[,c(20:25)])
    #Make the gene name the title of the columns
colnames(chr15qOrderCT2smDEL_plus1_log2_t)<-chr15qOrderCT2smDEL_plus1_log2$external_gene_name
    #Save final 15q gene expression table
write.csv(as.data.frame(chr15qOrderCT2smDEL_plus1_log2), file = "CT2smDEL_DESeqHISAT2_15qgenes.csv")
#Plot with ggplot
  #For H9smDEL
    #Reshape data
chr_meltH9smDEL <- melt(chr15qOrderH9smDEL_plus1_log2_t)
head(chr_meltH9smDEL)
    #Plot
H9smDELplot <- ggplot(chr_meltH9smDEL, aes(x = Var2, y = value))+ geom_boxplot() + geom_jitter(
    shape=16, position=position_jitter(0.2)) + theme(
      axis.text.x = element_text(angle = 60, hjust=1, size=6)) + labs(
        title="Chr15q11-q13 Gene Expression: H9smDEL", x ="Gene", y = "log2foldchange")
H9smDELplot
  #For CT2smDEL
    #Reshape data
chr_meltCT2smDEL <- melt(chr15qOrderCT2smDEL_plus1_log2_t)
head(chr_meltCT2smDEL)
    #Plot
CT2smDELplot <- ggplot(chr_meltCT2smDEL, aes(x = Var2, y = value))+ geom_boxplot() + geom_jitter(
    shape=16, position=position_jitter(0.2)) + theme(
      axis.text.x = element_text(angle = 60, hjust=1, size=6)) + labs(
        title="Chr15q11-q13 Gene Expression: CT2smDEL", x ="Gene", y = "log2foldchange")
CT2smDELplot
#Save as pdf
pdf("Chr15q11-q13_GeneExp_smDEL.pdf", width = 10, height = 8)
H9smDELplot
CT2smDELplot
dev.off()

#Combine saved plots using Adobe
