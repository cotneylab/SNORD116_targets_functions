##---DESeq---##
#Load relevant libraries
library("DESeq2")
library("data.table")
library("EnhancedVolcano")
library("dplyr")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("biomaRt")
library("RColorBrewer")
library("dendextend")
library("gplots")
library("ggplot2")
library("ggVennDiagram")
library("reshape2")
library("ggupset")
library("tidyr")
library("ComplexUpset")
library("devtools")
library("ggpubr")

#Set working directory
directory <- "../Cotney_Lab/PWS_RNASeq/STAR/lgDEL"
setwd(directory)

#Import file
counts <- read.delim("PWSRNAseq_STAR_featureCounts.txt", sep="\t", head=T, skip = 1, row.names = "Geneid")

#Simplify column names
colnames(counts) <- c("chr", "start", "end", "strand", "length", "CT2_1", "CT2_2", "CT2_3", "CT2_4", "CT2_5", "CT2_6", "CT2_lgDEL_1", "CT2_lgDEL_2", "CT2_lgDEL_3", "CT2_lgDEL_4", "CT2_lgDEL_5", "CT2_lgDEL_6", "CT2_smDEL_1", "CT2_smDEL_2", "CT2_smDEL_3", "CT2_smDEL_4", "CT2_smDEL_5", "CT2_smDEL_6", "H9_1", "H9_2", "H9_3", "H9_4", "H9_5", "H9_6", "H9_lgDEL_1", "H9_lgDEL_2", "H9_lgDEL_3", "H9_lgDEL_4", "H9_lgDEL_5", "H9_lgDEL_6", "H9_smDEL_1", "H9_smDEL_2", "H9_smDEL_3", "H9_smDEL_4", "H9_smDEL_5")
#Get rid of first 5 columns
counts <- counts[ ,6:ncol(counts)]
#Turn into matrix
counts <- as.matrix(counts)

#Load in metadata
colData <- as.data.frame(read.table("Metadata.txt", sep="\t", head=T, row.names = "ID"))
colData$Condition <- as.factor(colData$Condition)
colData$Genetic.Background <- as.factor(colData$Genetic.Background)
colData$Background_Condition <- as.factor(paste(colData$Genetic.Background, colData$Condition, sep="_"))

#Make objects with no smDEL for comparisons
counts_lgDEL <- as.data.frame(counts)
counts_lgDEL <- counts_lgDEL[-c(13:18,31:35)]
counts_lgDEL <- as.matrix(counts_lgDEL)
colData_lgDEL <- colData[-(13:18),]
colData_lgDEL <- colData_lgDEL[-(25:29),]

#Create DESeq objects
ddstotal <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ Background_Condition)
ddstotal
ddsCond <- DESeqDataSetFromMatrix(countData = counts_lgDEL, colData = colData_lgDEL, design = ~ Condition)
ddsCond
ddsCplusB <- DESeqDataSetFromMatrix(countData = counts_lgDEL, colData = colData_lgDEL, design = ~ Genetic.Background + Condition + Genetic.Background:Condition)
ddsCplusB

#Trim matrix, remove low counts
summary(rowSums(counts(ddstotal)))
summary(rowSums(counts(ddsCond)))
summary(rowSums(counts(ddsCplusB)))
#Filter using cutoff (can modify if necessary)
dds <- ddstotal [rowSums(counts(ddstotal)) > 1 ,]
ddsC <- ddsCond [rowSums(counts(ddsCond)) > 1 ,]
ddsCB <- ddsCplusB [rowSums(counts(ddsCplusB)) > 1 ,]

#Run DESeq on data
dds <- DESeq(dds)
ddsC <- DESeq(ddsC)
ddsCB <- DESeq(ddsCB)
#Save dds objects
save(dds, file = "PWSlgDEL_DESeqSTAR.RData")
save(ddsC, file = "PWSlgDEL_DESeqSTAR_ConditionOnly.RData")
save(ddsCB, file = "PWSlgDEL_DESeqSTAR_CplusB.RData")
  #To load
  load("PWSlgDEL_DESeqSTAR.RData")
  load("PWSlgDEL_DESeqSTAR_ConditionOnly.RData")
  load("PWSlgDEL_DESeqSTAR_CplusB.RData")

#Get log values
rld <- rlog(dds, blind=FALSE)
head(assay(rld))
hist(assay(rld))
rldC <- rlog(ddsC, blind=FALSE)
head(assay(rldC))
hist(assay(rldC))
rldCB <- rlog(ddsCB, blind=FALSE)
head(assay(rldCB))
hist(assay(rldCB))

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

plotPCA(rldC, intgroup = "Condition")
plotPCA(rldC, intgroup = "Genetic.Background")
plotPCA(rldCB, intgroup = "Condition")
plotPCA(rldCB, intgroup = "Genetic.Background")

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

pdf("PWS_lgDELvWT.pdf", width = 10, height = 8)
plotPCA(rldC, intgroup = "Condition")
plotPCA(rldC, intgroup = "Genetic.Background")
plotPCA(rldCB, intgroup = "Condition")
plotPCA(rldCB, intgroup = "Genetic.Background")
dev.off()

##---Data "Housekeeping"---##
#Removing decimal in gene name
rownames(rld) <- gsub('\\..+$', '', rownames(rld))
rownames(dds) <- gsub('\\..+$', '', rownames(dds))
rownames(rldC) <- gsub('\\..+$', '', rownames(rldC))
rownames(ddsC) <- gsub('\\..+$', '', rownames(ddsC))
rownames(rldCB) <- gsub('\\..+$', '', rownames(rldCB))
rownames(ddsCB) <- gsub('\\..+$', '', rownames(ddsCB))
#Check that it worked
counts(dds)
#Run results
resH9lgDEL <- results(dds, contrast = c("Background_Condition","H9_lgDEL","H9_WT"))
resH9lgDEL
resCT2lgDEL <- results(dds, contrast = c("Background_Condition","CT2_lgDEL","CT2_WT"))
resCT2lgDEL
C_reslgDEL <- results(ddsC, contrast = c("Condition", "lgDEL", "WT"))
C_reslgDEL
CB_reslgDEL <- results(ddsCB, contrast = c("Condition", "lgDEL", "WT"))
CB_reslgDEL

#Merge results table & counts table
  #For original (Background_Condition)
resDataH9lgDEL <- merge(as.data.frame(resH9lgDEL), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resDataCT2lgDEL <- merge(as.data.frame(resCT2lgDEL), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resDataH9lgDEL)[1] <- "Gene"
names(resDataCT2lgDEL)[1] <- "Gene"
rownames(resDataH9lgDEL) <- resDataH9lgDEL$Gene
rownames(resDataCT2lgDEL) <- resDataCT2lgDEL$Gene
  #For condition only
C_resDatalgDEL <- merge(as.data.frame(C_reslgDEL), as.data.frame(counts(ddsC, normalized=TRUE)), by="row.names", sort=FALSE)
names(C_resDatalgDEL)[1] <- "Gene"
rownames(C_resDatalgDEL) <- C_resDatalgDEL$Gene
  #For Genetic.Background + Condition + Genetic.Background:Condition
CB_resDatalgDEL <- merge(as.data.frame(CB_reslgDEL), as.data.frame(counts(ddsCB, normalized=TRUE)), by="row.names", sort=FALSE)
names(CB_resDatalgDEL)[1] <- "Gene"
rownames(CB_resDatalgDEL) <- CB_resDatalgDEL$Gene
#Change ESMBL ID's to gene names (use appropriate version based on GENCODE annotation used for alignment)
mart <- useMart("ensembl", host = "https://dec2017.archive.ensembl.org/", dataset = "hsapiens_gene_ensembl")
  #For H9lgDEL
R1 <- data.table(resDataH9lgDEL)
rownames(R1) <- R1$Gene
genes.table_H9lgDEL <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","description","gene_biotype","chromosome_name","start_position","end_position","strand"), values= R1$Gene, mart= mart)
  #For CT2lgDEL
R2 <- data.table(resDataCT2lgDEL)
rownames(R2) <- R2$Gene
genes.table_CT2lgDEL <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","description","gene_biotype","chromosome_name","start_position","end_position","strand"), values= R2$Gene, mart= mart)
  #For condition only
RC <- data.table(C_resDatalgDEL)
rownames(RC) <- RC$Gene
genes.table_lgDEL_C <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","description","gene_biotype","chromosome_name","start_position","end_position","strand"), values= RC$Gene, mart= mart)
  #For Genetic.Background + Condition + Genetic.Background:Condition
RCB <- data.table(CB_resDatalgDEL)
rownames(RCB) <- RCB$Gene
genes.table_lgDEL_CB <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","description","gene_biotype","chromosome_name","start_position","end_position","strand"), values= RCB$Gene, mart= mart)
#Order lists by ENSEMBL ID so they merge correctly
  #For original (Background_Condition)
resOrderedgene_H9lgDEL <- R1[order(R1$Gene),]
names(genes.table_H9lgDEL)[1] <- "Gene"
resOrderedgene_CT2lgDEL <- R2[order(R2$Gene),]
names(genes.table_CT2lgDEL)[1] <- "Gene"
  #For condition
C_resOrderedgene_lgDEL <- RC[order(RC$Gene),]
names(genes.table_lgDEL_C)[1] <- "Gene"
  #For Genetic.Background + Condition + Genetic.Background:Condition
CB_resOrderedgene_lgDEL <- RCB[order(RCB$Gene),]
names(genes.table_lgDEL_CB)[1] <- "Gene"
#Merge all things
  #For original (Background_Condition)
ResGeneID_H9lgDEL <- merge(as.data.frame(resOrderedgene_H9lgDEL), as.data.frame(genes.table_H9lgDEL), by="Gene", sort=TRUE, all = FALSE)
rownames(ResGeneID_H9lgDEL) <- ResGeneID_H9lgDEL$Gene
ResGeneID_CT2lgDEL <- merge(as.data.frame(resOrderedgene_CT2lgDEL), as.data.frame(genes.table_CT2lgDEL), by="Gene", sort=TRUE, all = FALSE)
rownames(ResGeneID_CT2lgDEL) <- ResGeneID_CT2lgDEL$Gene
  #For condition only
ResGeneID_lgDEL_C <- merge(as.data.frame(C_resOrderedgene_lgDEL), as.data.frame(genes.table_lgDEL_C), by="Gene", sort=TRUE, all = FALSE)
rownames(ResGeneID_lgDEL_C) <- ResGeneID_lgDEL_C$Gene
  #For Genetic.Background + Condition + Genetic.Background:Condition
ResGeneID_lgDEL_CB <- merge(as.data.frame(CB_resOrderedgene_lgDEL), as.data.frame(genes.table_lgDEL_CB), by="Gene", sort=TRUE, all = FALSE)
rownames(ResGeneID_lgDEL_CB) <- ResGeneID_lgDEL_CB$Gene
#Fill in blank external gene names
  #Original (Background_Condition)
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
  #Condition only
ResGeneID_lgDEL_C$external_gene_name <- as.character(ResGeneID_lgDEL_C$external_gene_name)
ResGeneID_lgDEL_C$external_gene_name[ResGeneID_lgDEL_C$external_gene_name==""] <- NA
ResGeneID_lgDEL_C$external_gene_name <- as.factor(ResGeneID_lgDEL_C$external_gene_name)
GeneSymbolDE_lgDEL_C <- mapIds(org.Hs.eg.db, keys = ResGeneID_lgDEL_C$Gene, column = "SYMBOL", keytype = "ENSEMBL", multivals= "first")
for (i in c(1:length(GeneSymbolDE_lgDEL_C))){
  
  if (is.na(GeneSymbolDE_lgDEL_C[i])){
    
    GeneSymbolDE_lgDEL_C[i]<-gsub("\\..*","",ResGeneID_lgDEL_C$Gene[i])
    
  }
  
}
ResGeneID_lgDEL_C$external_gene_name <- GeneSymbolDE_lgDEL_C
  #For Genetic.Background + Condition + Genetic.Background:Condition
ResGeneID_lgDEL_CB$external_gene_name <- as.character(ResGeneID_lgDEL_CB$external_gene_name)
ResGeneID_lgDEL_CB$external_gene_name[ResGeneID_lgDEL_CB$external_gene_name==""] <- NA
ResGeneID_lgDEL_CB$external_gene_name <- as.factor(ResGeneID_lgDEL_CB$external_gene_name)
GeneSymbolDE_lgDEL_CB <- mapIds(org.Hs.eg.db, keys = ResGeneID_lgDEL_CB$Gene, column = "SYMBOL", keytype = "ENSEMBL", multivals= "first")
for (i in c(1:length(GeneSymbolDE_lgDEL_CB))){
  
  if (is.na(GeneSymbolDE_lgDEL_CB[i])){
    
    GeneSymbolDE_lgDEL_CB[i]<-gsub("\\..*","",ResGeneID_lgDEL_CB$Gene[i])
    
  }
  
}
ResGeneID_lgDEL_CB$external_gene_name <- GeneSymbolDE_lgDEL_CB
#For Original (Background_Condition)
  #Merge both results objects
ResGeneID_ALL <- merge(as.data.frame(ResGeneID_H9lgDEL), as.data.frame(ResGeneID_CT2lgDEL), by="Gene", sort=TRUE, all = FALSE)
rownames(ResGeneID_ALL) <- ResGeneID_ALL$Gene
  #Drop unnecessary columns
ResGeneID_ALL <- ResGeneID_ALL[-c(56:97)]
  #Simplify column names
colnames(ResGeneID_ALL) <- c("Gene", "baseMean_H9", "log2FoldChange_H9", "lfcSE_H9", "stat_H9", "pvalue_H9", "padj_H9", "CT2_1", "CT2_2", "CT2_3", "CT2_4", "CT2_5", "CT2_6", "CT2_lgDEL_1", "CT2_lgDEL_2", "CT2_lgDEL_3", "CT2_lgDEL_4", "CT2_lgDEL_5", "CT2_lgDEL_6", "CT2_smDEL_1", "CT2_smDEL_2", "CT2_smDEL_3", "CT2_smDEL_4", "CT2_smDEL_5", "CT2_smDEL_6", "H9_1", "H9_2", "H9_3", "H9_4", "H9_5", "H9_6", "H9_lgDEL_1", "H9_lgDEL_2", "H9_lgDEL_3", "H9_lgDEL_4", "H9_lgDEL_5", "H9_lgDEL_6", "H9_smDEL_1", "H9_smDEL_2", "H9_smDEL_3", "H9_smDEL_4", "H9_smDEL_5", "external_gene_name", "description", "gene_biotype", "chromosome_name", "start_position", "end_position", "strand", "baseMean_CT2", "log2FoldChange_CT2", "lfcSE_CT2", "stat_CT2", "pvalue_CT2", "padj_CT2")
#Save summary tables
write.csv(as.data.frame(ResGeneID_ALL), file = "CombinedH9CT2_Results_DESeqSTAR.csv")
write.csv(as.data.frame(ResGeneID_lgDEL_C), file = "Results_DESeqSTAR_conditiononly.csv")
write.csv(as.data.frame(ResGeneID_lgDEL_CB), file = "Results_DESeqSTAR_condplusbkgrnd.csv")
  #To read in tables
  ResGeneID_ALL <- read.csv(file = 'CombinedH9CT2_Results_DESeqSTAR.csv', header = TRUE, row.names = "X")
  ResGeneID_lgDEL_C <- read.csv(file = 'Results_DESeqSTAR_conditiononly.csv', header = TRUE, row.names = "X")
  ResGeneID_lgDEL_CB <- read.csv(file = 'Results_DESeqSTAR_condplusbkgrnd.csv', header = TRUE, row.names = "X")

##---Venn Diagram---##
  #For original (Background_Condition)
H9lgDELsig <- ResGeneID_ALL[which(ResGeneID_ALL$padj_H9 < 0.05),]
upH9lgDEL <- H9lgDELsig[which(H9lgDELsig$log2FoldChange_H9 > 0),]
downH9lgDEL <- H9lgDELsig[which(H9lgDELsig$log2FoldChange_H9 < 0),]
CT2lgDELsig <- ResGeneID_ALL[which(ResGeneID_ALL$padj_CT2 < 0.05),]
upCT2lgDEL <- CT2lgDELsig[which(CT2lgDELsig$log2FoldChange_CT2 > 0),]
downCT2lgDEL <- CT2lgDELsig[which(CT2lgDELsig$log2FoldChange_CT2 < 0),]
gene_list <- list(A = sample(downH9lgDEL$Gene),
                  B = sample(downCT2lgDEL$Gene),
                  C = sample(upCT2lgDEL$Gene),
                  D = sample(upH9lgDEL$Gene))
p1 <- ggVennDiagram(gene_list, 
                    category.names = c("DownH9lgDEL","DownCT2lgDEL","UpCT2lgDEL","UpH9lgDEL"),
                    label = "count")
p1
#Expand axis to show long labels
p1 + scale_x_continuous(expand = expansion(mult = .2))
#Save as pdf
pdf("VennDiagramlgDEL.pdf", width = 10, height = 8)
p1
dev.off()
#Pull out intersection values
intersect_Venn_ENSEMBL <- process_region_data(Venn(gene_list))
#Save as csv
intersect_Venn_ENSEMBL <- as.data.frame(intersect_Venn_ENSEMBL)
intersect_Venn_ENSEMBL <- apply(intersect_Venn_ENSEMBL,2,as.character)
write.csv(intersect_Venn_ENSEMBL, file = "VennDiagramlgDEL_IntersectResultsENSEMBL.csv")
#Subset & save lists
down_lgDEL_ENSEMBL <- list(intersect_Venn_ENSEMBL[5,3])
lapply(down_lgDEL_ENSEMBL, write, "VennDiagramlgDEL_IntersectDownResultsENSEMBL.txt", append=TRUE, ncolumns=1000)
up_lgDEL_ENSEMBL <- list(intersect_Venn_ENSEMBL[10,3])
lapply(up_lgDEL_ENSEMBL, write, "VennDiagramlgDEL_IntersectUpResultsENSEMBL.txt", append=TRUE, ncolumns=1000)

##Test for significance of overlaps
#Create gene lists
all_H9_genes <- list(resDataH9lgDEL$Gene)
write.table(all_H9_genes, "AllH9Genes.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
all_CT2_genes <- list(resDataCT2lgDEL$Gene)
write.table(all_CT2_genes, "AllCT2Genes.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
#Copy text files onto the cluster and run venn_intersect_test_rbg.sbatch (Done with ENSEMBL numbers)
#Copy files from permutations to computer and read in
A_B <- read.table(file='A_B.txt', header=FALSE)
C_D <- read.table(file='C_D.txt', header=FALSE)
#Calculate median value
median(A_B$V1)
median(C_D$V1)
#Plot histograms
hist_AB <- ggplot(A_B, aes(x=V1)) + 
  geom_histogram(binwidth=2, color="black", fill="white") +
  geom_vline(aes(xintercept = median(V1)), col='#224930', linewidth=1, linetype="dashed", alpha=0.8) +
  geom_vline(aes(xintercept = 381), col='#005AB5', linewidth=1) +
  scale_x_continuous(breaks = seq(0, 400, 20)) +
  theme_bw() +
  labs(title="Histogram of Downregulated Overlap Permutation Test",x="Number of Overlapping Genes", y = "Frequency")
hist_AB
hist_CD <- ggplot(C_D, aes(x=V1)) + 
  geom_histogram(binwidth=2, color="black", fill="white") +
  geom_vline(aes(xintercept = median(V1)), col='#224930', linewidth=1, linetype="dashed", alpha=0.8) +
  geom_vline(aes(xintercept = 483), col='#DC3220', linewidth=1) +
  scale_x_continuous(breaks = seq(0, 480, 20)) +
  theme_bw() +
  labs(title="Histogram of Upregulated Overlap Permutation",x="Number of Overlapping Genes", y = "Frequency")
hist_CD
#Save as pdf
pdf("VennOverlapPermutationsHistogramlgDEL.pdf", width = 10, height = 8)
hist_AB
hist_CD
dev.off()

##---Data "Housekeeping" part 2---##
#Take intersect gene lists and use Notepad++ to make them into a list with one gene per row (Find & Replace ", " with \n)
#Read in edited text files as the objects
down_lgDEL_ENSEMBL <- as.data.frame(read.table("VennDiagramlgDEL_IntersectDownResultsENSEMBL.txt", sep="\t", head=F))
names(down_lgDEL_ENSEMBL)[1] <- "Gene"
up_lgDEL_ENSEMBL <- as.data.frame(read.table("VennDiagramlgDEL_IntersectUpResultsENSEMBL.txt", sep="\t", head=F))
names(up_lgDEL_ENSEMBL)[1] <- "Gene"
#Find external_gene_names
down_lgDEL <- ResGeneID_ALL[down_lgDEL_ENSEMBL$Gene,]
down_lgDEL <- as.data.frame(down_lgDEL$external_gene_name)
names(down_lgDEL)[1] <- "external_gene_name"
down_lgDEL <- down_lgDEL[order(down_lgDEL$external_gene_name),]
up_lgDEL <- ResGeneID_ALL[up_lgDEL_ENSEMBL$Gene,]
up_lgDEL <- as.data.frame(up_lgDEL$external_gene_name)
names(up_lgDEL)[1] <- "external_gene_name"
up_lgDEL <- up_lgDEL[order(up_lgDEL$external_gene_name),]
  #Save lists
  write.table(down_lgDEL, "VennDiagram_IntersectDownResults.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(up_lgDEL, "VennDiagram_IntersectUpResults.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
#Combine lists of genes at intersects to use downstream
GeneName_Intersect <- as.data.frame(append(up_lgDEL, down_lgDEL))
names(GeneName_Intersect)[1] <- "external_gene_name"
GeneID_Intersect <- rbind(up_lgDEL_ENSEMBL, down_lgDEL_ENSEMBL)
#Order list
GeneID_Intersect <- GeneID_Intersect[order(GeneID_Intersect$Gene),]
#Parse results file by intersect list
ResGeneID_ALL <- ResGeneID_ALL[order(ResGeneID_ALL$Gene),]
ResGene_ALL_sharedID <- ResGeneID_ALL[GeneID_Intersect,]

##---Volcano Plot---##
#For Original (Background_Condition)
  #Calculate average log2FoldChange
ResGene_ALL_sharedID$log2FoldChange_avg<-rowMeans(ResGene_ALL_sharedID[,c(3,51)])
  #Make volcano plots
EnVol_H9 <- EnhancedVolcano(ResGene_ALL_sharedID,
                            lab = ResGene_ALL_sharedID$external_gene_name,
                            x = 'log2FoldChange_avg',
                            xlim = c(-10,10),
                            FCcutoff = 0,
                            y = 'padj_H9',
                            ylim = c(0,100),
                            pCutoff = 5e-2,
                            title = 'H9lgDEL vs H9 WT',
                            subtitle = 'DESeq2 Results of Shared Genes',
                            points(x = -10, y = 100, pch = 19, clip = FALSE))
EnVol_H9
EnVol_CT2 <- EnhancedVolcano(ResGene_ALL_sharedID,
                             lab = ResGene_ALL_sharedID$external_gene_name,
                             x = 'log2FoldChange_avg',
                             xlim = c(-15,10),
                             FCcutoff = 0,
                             y = 'padj_CT2',
                             ylim = c(0,100),
                             pCutoff = 5e-2,
                             title = 'CT2lgDEL vs CT2 WT',
                             subtitle = 'DESeq2 Results')
EnVol_CT2
  #Save as pdf
pdf("Volc_SharedGeneslgDEL.pdf", width = 10, height = 8)
EnVol_H9
EnVol_CT2
dev.off()
##Creating a "prettier" plot for figures
  #Copy results object to new object
res.volc <- ResGene_ALL_sharedID
  #Calculate p.adj for -log10(padj) = 100
1/(10^100) = 1e-100
  #Set cutoff(s)
x_cutoff <- 5
y_cutoff <- 1e-100
  #Identify genes that pass the threshold(s)
up_x <- rownames(subset(res.volc, log2FoldChange_avg >= x_cutoff))
up_x
down_x <- rownames(subset(res.volc, log2FoldChange_avg <= x_cutoff * -1))
down_x
y_H9 <- rownames(subset(res.volc, padj_H9 <= y_cutoff))
y_H9
y_CT2 <- rownames(subset(res.volc, padj_CT2 <= y_cutoff))
y_CT2
  #Impute log2FCs/padj in the object to be +/- cutoff  
res.volc$log2FoldChange_avg <- ifelse(res.volc$log2FoldChange_avg >= x_cutoff, x_cutoff,
                                      ifelse(res.volc$log2FoldChange_avg <= x_cutoff * -1, x_cutoff * -1,
                                             res.volc$log2FoldChange_avg))
res.volc$padj_H9 <- ifelse(res.volc$padj_H9 <= y_cutoff, y_cutoff, res.volc$padj_H9)
res.volc$padj_CT2 <- ifelse(res.volc$padj_CT2 <= y_cutoff, y_cutoff, res.volc$padj_CT2)
  #Create custom shapes for the points which fall out of plot window
customshape <- rep(19, nrow(res.volc))
names(customshape) <- rep('group1', nrow(res.volc))
customshape[which(rownames(res.volc) %in% up_x)] <- -9658
names(customshape)[which(rownames(res.volc) %in% up_x)] <- 'group2'
customshape[which(rownames(res.volc) %in% down_x)] <- -9668
names(customshape)[which(rownames(res.volc) %in% down_x)] <- 'group3'
customshape[which(rownames(res.volc) %in% y_H9)] <- 17
names(customshape)[which(rownames(res.volc) %in% y_H9)] <- 'group4'
  #Creat custom sizes for the points
customsize <- rep(2.0, nrow(res.volc))
customsize[which(rownames(res.volc) %in% up_x)] <- 5
customsize[which(rownames(res.volc) %in% down_x)] <- 5
customsize[which(rownames(res.volc) %in% y_H9)] <- 4
  #Plot
EnVol_H9_fig <- EnhancedVolcano(res.volc,
                                lab = res.volc$external_gene_name,
                                labSize = 3,
                                drawConnectors = TRUE,
                                arrowheads = FALSE,
                                x = 'log2FoldChange_avg',
                                y = 'padj_H9',
                                xlab = bquote("avg"~Log[2] ~ "(foldChange)"),
                                ylab = bquote(~-Log[10] ~ (P.adj)),
                                shapeCustom = customshape,
                                pointSize = customsize,
                                col = c("gray", "darkgray", "#FFB47E", "#AB5000"),
                                FCcutoff = 1,
                                pCutoff = 0.05,
                                xlim = c(-5,5),
                                ylim = c(0,100),
                                legendPosition = 'top',
                                legendIconSize = 3,
                                legendLabSize = 12,
                                title = 'H9-lgDEL vs H9-WT',
                                subtitle = 'DESeq2 Results: DEGs Shared Between H9 & CT2 lgDEL')
  #Create custom shapes for the points which fall out of plot window
customshape <- rep(19, nrow(res.volc))
names(customshape) <- rep('group1', nrow(res.volc))
customshape[which(rownames(res.volc) %in% up_x)] <- -9658
names(customshape)[which(rownames(res.volc) %in% up_x)] <- 'group2'
customshape[which(rownames(res.volc) %in% down_x)] <- -9668
names(customshape)[which(rownames(res.volc) %in% down_x)] <- 'group3'
customshape[which(rownames(res.volc) %in% y_CT2)] <- 17
names(customshape)[which(rownames(res.volc) %in% y_CT2)] <- 'group4'
  #Creat custom sizes for the points
customsize <- rep(2.0, nrow(res.volc))
customsize[which(rownames(res.volc) %in% up_x)] <- 5
customsize[which(rownames(res.volc) %in% down_x)] <- 5
customsize[which(rownames(res.volc) %in% y_CT2)] <- 4
  #Plot
EnVol_CT2_fig <- EnhancedVolcano(res.volc,
                                lab = res.volc$external_gene_name,
                                labSize = 3,
                                drawConnectors = TRUE,
                                arrowheads = FALSE,
                                x = 'log2FoldChange_avg',
                                y = 'padj_CT2',
                                xlab = bquote("avg"~Log[2] ~ "(foldChange)"),
                                ylab = bquote(~-Log[10] ~ (P.adj)),
                                shapeCustom = customshape,
                                pointSize = customsize,
                                col = c("gray", "darkgray", "#CFB4FF", "#5D3A9B"),
                                FCcutoff = 1,
                                pCutoff = 0.05,
                                xlim = c(-5,5),
                                ylim = c(0,100),
                                legendPosition = 'top',
                                legendIconSize = 3,
                                legendLabSize = 12,
                                title = 'CT2-lgDEL vs CT2-WT',
                                subtitle = 'DESeq2 Results: DEGs Shared Between H9 & CT2 lgDEL')
  #Save pdf
grDevices::cairo_pdf("Volc_SharedGeneslgDEL_figure.pdf", width = 10, height = 8, onefile = TRUE)
EnVol_H9_fig
EnVol_CT2_fig
dev.off()

#Create list of genes from our region on chromosome 15 to put in bold on volcano plots
  #Read in file of genes in our region
chr15list <- read.csv("../../HISAT2/lgDEL/chr15q11-q13_gencodev25.csv")
  #Rename columns
colnames(chr15list) <- c("transcript_id", "name2", "ensembl_gene_id", "gene_name")
  #Remove decimal in gene_id
chr15list$ensembl_gene_id <- gsub('\\..+$', '', chr15list$ensembl_gene_id)
  #Create list of unique ENSEMBL genes
chr15list <- unique(chr15list$ensembl_gene_id)
  #Parse sharedID list
chr15_overlap <- ResGene_ALL_sharedID[(ResGene_ALL_sharedID$Gene) %in% chr15list, ]
  #Pull out all names needing to be bolded
chr15_overlap_list <- chr15_overlap$external_gene_name
  #Write out list as text file
lapply(chr15_overlap_list, write, "chr15_genes_onVolc.txt", append=TRUE, ncolumns=1000)
##Back to regularly scheduled Volcano Plots :)
#For condition only
  #Make volcano plot
EnVol_C <- EnhancedVolcano(ResGeneID_lgDEL_C,
                            lab = ResGeneID_lgDEL_C$external_gene_name,
                            x = 'log2FoldChange',
                            FCcutoff = 0,
                            y = 'padj',
                            title = 'lgDEL vs WT',
                            subtitle = 'Condition only')
EnVol_C
#For Genetic.Background + Condition + Genetic.Background:Condition
  #Make volcano plot
EnVol_CB <- EnhancedVolcano(ResGeneID_lgDEL_CB,
                           lab = ResGeneID_lgDEL_CB$external_gene_name,
                           x = 'log2FoldChange',
                           FCcutoff = 0,
                           y = 'padj',
                           title = 'lgDEL vs WT',
                           subtitle = 'Genetic.Background + Condition + Genetic.Background:Condition')
EnVol_CB
#Save as pdf
pdf("Volc_lgDELanalysis.pdf", width = 10, height = 8)
EnVol_C
EnVol_CB
dev.off()

##---Heatmap: All DE genes---##
#All shared DE genes - Original (Background_Condition)
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
#Top DE genes - Original (Background_Condition)
resOrdered_shared <- ResGene_ALL_sharedID[order(ResGene_ALL_sharedID$log2FoldChange_avg),]
resOrdered_sharedPCG <- resOrdered_shared %>% filter(resOrdered_shared$gene_biotype %in% "protein_coding")
topGenes_neg <- head(rownames(resOrdered_sharedPCG),25)
topGenes_pos <- tail(rownames(resOrdered_sharedPCG),25)
topGenes_shared <- rbind(topGenes_pos, topGenes_neg)
rld4heatmapde_top <- assay(rld)[topGenes_shared,]
rld4heatmapde_top <- rld4heatmapde_top[, !colnames(rld4heatmapde_top) %in% c("CT2_smDEL_1", "CT2_smDEL_2", "CT2_smDEL_3", "CT2_smDEL_4", "CT2_smDEL_5", "CT2_smDEL_6", "H9_smDEL_1", "H9_smDEL_2", "H9_smDEL_3", "H9_smDEL_4", "H9_smDEL_5")]
hr_shared <- hclust(as.dist(1-cor(t(rld4heatmapde_top), method="pearson")), method="complete")
df <- as.data.frame(hr_shared[["labels"]])
colnames(df)[1] <- "Gene"
df <- ResGene_ALL_sharedID %>% arrange(factor(Gene, levels = df$Gene))
df <- df[1:50,]
rownames(rld4heatmapde_top) <- df$external_gene_name
rld4heatmapde_H9 <- rld4heatmapde_top[, colnames(rld4heatmapde_top) %in% c("H9_1", "H9_2", "H9_3", "H9_4", "H9_5", "H9_6", "H9_lgDEL_1", "H9_lgDEL_2", "H9_lgDEL_3", "H9_lgDEL_4", "H9_lgDEL_5", "H9_lgDEL_6")]
rld4heatmapde_CT2 <- rld4heatmapde_top[, !colnames(rld4heatmapde_top) %in% c("H9_1", "H9_2", "H9_3", "H9_4", "H9_5", "H9_6", "H9_lgDEL_1", "H9_lgDEL_2", "H9_lgDEL_3", "H9_lgDEL_4", "H9_lgDEL_5", "H9_lgDEL_6")]
hc_H9 <- hclust(as.dist(1-cor(rld4heatmapde_H9, method="spearman")), method="complete")
hc_CT2 <- hclust(as.dist(1-cor(rld4heatmapde_CT2, method="spearman")), method="complete")
heatmap.2( rld4heatmapde_H9, key=T, scale="row", 
           trace="none", Rowv=as.dendrogram(hr_shared), Colv=as.dendrogram(hc_H9), dendrogram = "row", cexRow = .4, cexCol = .5,
           col = colorRampPalette(c("blue", "white", "red"))(n = 1000))
heatmap.2( rld4heatmapde_CT2, key=T, scale="row", 
           trace="none", Rowv=as.dendrogram(hr_shared), Colv=as.dendrogram(hc_CT2), dendrogram = "row", cexRow = .4, cexCol = .5,
           col = colorRampPalette(c("blue", "white", "red"))(n = 1000))

pdf("Top50Heatmap_skinny.pdf", width = 4, height = 9)
heatmap.2( rld4heatmapde_H9, key=T, scale="row", 
           trace="none", Rowv=as.dendrogram(hr_shared), Colv=as.dendrogram(hc_H9), dendrogram = "row", cexRow = .4, cexCol = .5,
           col = colorRampPalette(c("blue", "white", "red"))(n = 1000))
heatmap.2( rld4heatmapde_CT2, key=T, scale="row", 
           trace="none", Rowv=as.dendrogram(hr_shared), Colv=as.dendrogram(hc_CT2), dendrogram = "row", cexRow = .4, cexCol = .5,
           col = colorRampPalette(c("blue", "white", "red"))(n = 1000))
dev.off()

pdf("Top50Heatmap_reg.pdf", width = 10, height = 8)
heatmap.2( rld4heatmapde_H9, key=T, scale="row", 
           trace="none", Rowv=as.dendrogram(hr_shared), Colv=as.dendrogram(hc_H9), dendrogram = "row", cexRow = .4, cexCol = .5,
           col = colorRampPalette(c("blue", "white", "red"))(n = 1000))
heatmap.2( rld4heatmapde_CT2, key=T, scale="row", 
           trace="none", Rowv=as.dendrogram(hr_shared), Colv=as.dendrogram(hc_CT2), dendrogram = "row", cexRow = .4, cexCol = .5,
           col = colorRampPalette(c("blue", "white", "red"))(n = 1000))
dev.off()

##---Upset Plots---##
##For determining overlap of Background_Condition (same analysis as VennDiagram above, just showing it differently)
#Make lists of genes
downH9_gene <- as.data.frame(downH9lgDEL$Gene)
downCT2_gene <- as.data.frame(downCT2lgDEL$Gene)
upH9_gene <- as.data.frame(upH9lgDEL$Gene)
upCT2_gene <- as.data.frame(upCT2lgDEL$Gene)
#Change first column name
colnames(downH9_gene)[1] <- "Gene"
colnames(downCT2_gene)[1] <- "Gene"
colnames(upH9_gene)[1] <- "Gene"
colnames(upCT2_gene)[1] <- "Gene"
#Make column for TRUE/FALSE
downH9_gene$down_H9 <- TRUE
downCT2_gene$down_CT2 <- TRUE
upH9_gene$up_H9 <- TRUE
upCT2_gene$up_CT2 <- TRUE
#Merge lists
all_upset_df <- merge(downH9_gene, downCT2_gene, by = 'Gene', all = TRUE)
all_upset_df <- merge(all_upset_df, upH9_gene, by = 'Gene', all = TRUE)
all_upset_df <- merge(all_upset_df, upCT2_gene, by = 'Gene', all = TRUE)
#Set NA values to FALSE
all_upset_df[is.na(all_upset_df)] <- FALSE
#Plot
all_plot <- upset(all_upset_df, colnames(all_upset_df)[2:5], name = "DE Genes (p.adj<0.05)", width_ratio = 0.1, keep_empty_groups=TRUE,
                  matrix=(intersection_matrix(geom=geom_point(shape='square', size=3.5), segment=geom_segment(linetype='dotted'), outline_color=list(active='#FE6100', inactive='#929292'))
                          + scale_color_manual(values=c('TRUE'='#FFB000', 'FALSE'='#B5B4B4'), labels=c('TRUE'='yes', 'FALSE'='no'), breaks=c('TRUE', 'FALSE'), name='Intersection')
                          + scale_y_discrete(position='right')
                          + annotate(geom='text', label='Look here →', x='Comedy-Drama', y='Drama', size=5, hjust=1)),
                  queries=list(upset_query(intersect=c('down_H9', 'down_CT2'), color='#005AB5', fill='#005AB5', only_components=c('intersections_matrix', 'Intersection size')), upset_query(intersect=c('up_H9', 'up_CT2'), color='#DC3220', fill='#DC3220', only_components=c('intersections_matrix', 'Intersection size')))) + 
  ggtitle('Shared DE Genes in H9 & CT2 lgDEL vs WT')
all_plot
#Save plots as pdf
pdf("lgDEL_sharedGenes_Upset.pdf", width = 10, height = 8)
all_plot
dev.off()
##For determining overlap of 3 different analyses (results contrasts/specifications)
#Parse results files
sig <- ResGeneID_ALL[which(ResGeneID_ALL$padj_H9 < 0.05 & ResGeneID_ALL$padj_CT2 < 0.05),]
up <- sig[which(sig$log2FoldChange_H9 > 0 & sig$log2FoldChange_CT2 > 0),]
down <- sig[which(sig$log2FoldChange_H9 < 0 & sig$log2FoldChange_CT2 < 0),]
CondOnly_sig <- ResGeneID_lgDEL_C[which(ResGeneID_lgDEL_C$padj < 0.05),]
up_CondOnly <- CondOnly_sig[which(CondOnly_sig$log2FoldChange > 0),]
down_CondOnly <- CondOnly_sig[which(CondOnly_sig$log2FoldChange < 0),]
CplusB_sig <- ResGeneID_lgDEL_CB[which(ResGeneID_lgDEL_CB$padj < 0.05),]
up_CplusB <- CplusB_sig[which(CplusB_sig$log2FoldChange > 0),]
down_CplusB <- CplusB_sig[which(CplusB_sig$log2FoldChange < 0),]
#Make lists of genes
allSig_gene <- as.data.frame(sig$Gene)
allSig_CondOnly_gene <- as.data.frame(CondOnly_sig$Gene)
allSig_CplusB_gene <- as.data.frame(CplusB_sig$Gene)
down_gene <- as.data.frame(down$Gene)
down_CondOnly_gene <- as.data.frame(down_CondOnly$Gene)
down_CplusB_gene <- as.data.frame(down_CplusB$Gene)
up_gene <- as.data.frame(up$Gene)
up_CondOnly_gene <- as.data.frame(up_CondOnly$Gene)
up_CplusB_gene <- as.data.frame(up_CplusB$Gene)
#Change first column name
colnames(allSig_gene)[1] <- "Gene"
colnames(allSig_CondOnly_gene)[1] <- "Gene"
colnames(allSig_CplusB_gene)[1] <- "Gene"
colnames(down_gene)[1] <- "Gene"
colnames(down_CondOnly_gene)[1] <- "Gene"
colnames(down_CplusB_gene)[1] <- "Gene"
colnames(up_gene)[1] <- "Gene"
colnames(up_CondOnly_gene)[1] <- "Gene"
colnames(up_CplusB_gene)[1] <- "Gene"
#Make column for TRUE/FALSE
allSig_gene$C_B <- TRUE
allSig_CondOnly_gene$ConditionOnly <- TRUE
allSig_CplusB_gene$CplusB <- TRUE
down_gene$C_B <- TRUE
down_CondOnly_gene$ConditionOnly <- TRUE
down_CplusB_gene$CplusB <- TRUE
up_gene$C_B <- TRUE
up_CondOnly_gene$ConditionOnly <- TRUE
up_CplusB_gene$CplusB <- TRUE
#Merge lists
all_upset_df <- merge(allSig_gene, allSig_CondOnly_gene, by = 'Gene', all = TRUE)
all_upset_df <- merge(all_upset_df, allSig_CplusB_gene, by = 'Gene', all = TRUE)
down_upset_df <- merge(down_gene, down_CondOnly_gene, by = 'Gene', all = TRUE)
down_upset_df <- merge(down_upset_df, down_CplusB_gene, by = 'Gene', all = TRUE)
up_upset_df <- merge(up_gene, up_CondOnly_gene, by = 'Gene', all = TRUE)
up_upset_df <- merge(up_upset_df, up_CplusB_gene, by = 'Gene', all = TRUE)
#Set NA values to FALSE
all_upset_df[is.na(all_upset_df)] <- FALSE
down_upset_df[is.na(down_upset_df)] <- FALSE
up_upset_df[is.na(up_upset_df)] <- FALSE
#Plot
all_plot <- upset(all_upset_df, colnames(all_upset_df)[2:4], name = "DE Genes (p<0.05)", width_ratio = 0.1, keep_empty_groups=TRUE,
                  matrix=(intersection_matrix(geom=geom_point(shape='square', size=3.5), segment=geom_segment(linetype='dotted'), outline_color=list(active='#FE6100', inactive='#929292'))
                          + scale_color_manual(values=c('TRUE'='#FFB000', 'FALSE'='#B5B4B4'), labels=c('TRUE'='yes', 'FALSE'='no'), breaks=c('TRUE', 'FALSE'), name='Intersection')
                          + scale_y_discrete(position='right')
                          + annotate(geom='text', label='Look here →', x='Comedy-Drama', y='Drama', size=5, hjust=1)),
                  queries=list(upset_query(intersect=c('C_B', 'ConditionOnly', 'CplusB'), color='#DC267F', fill='#DC267F', only_components=c('intersections_matrix', 'Intersection size')))) + 
  ggtitle('All Shared DE Genes (Disregarding Directionality)')
all_plot
up_plot <- upset(up_upset_df, colnames(up_upset_df)[2:4], name = "DE Genes (p<0.05)", width_ratio = 0.1, keep_empty_groups=TRUE,
                 matrix=(intersection_matrix(geom=geom_point(shape='square', size=3.5), segment=geom_segment(linetype='dotted'), outline_color=list(active='#FE6100', inactive='#929292'))
                         + scale_color_manual(values=c('TRUE'='#FFB000', 'FALSE'='#B5B4B4'), labels=c('TRUE'='yes', 'FALSE'='no'), breaks=c('TRUE', 'FALSE'), name='Intersection')
                         + scale_y_discrete(position='right')
                         + annotate(geom='text', label='Look here →', x='Comedy-Drama', y='Drama', size=5, hjust=1)),
                 queries=list(upset_query(intersect=c('C_B', 'ConditionOnly', 'CplusB'), color='#DC267F', fill='#DC267F', only_components=c('intersections_matrix', 'Intersection size')))) + 
  ggtitle('Upregulated Shared DE Genes')
up_plot
down_plot <- upset(down_upset_df, colnames(down_upset_df)[2:4], name = "DE Genes (p<0.05)", width_ratio = 0.1, keep_empty_groups=TRUE,
                   matrix=(intersection_matrix(geom=geom_point(shape='square', size=3.5), segment=geom_segment(linetype='dotted'), outline_color=list(active='#FE6100', inactive='#929292'))
                           + scale_color_manual(values=c('TRUE'='#FFB000', 'FALSE'='#B5B4B4'), labels=c('TRUE'='yes', 'FALSE'='no'), breaks=c('TRUE', 'FALSE'), name='Intersection')
                           + scale_y_discrete(position='right')
                           + annotate(geom='text', label='Look here →', x='Comedy-Drama', y='Drama', size=5, hjust=1)),
                   queries=list(upset_query(intersect=c('C_B', 'ConditionOnly', 'CplusB'), color='#DC267F', fill='#DC267F', only_components=c('intersections_matrix', 'Intersection size')))) + 
  ggtitle('Downregulated Shared DE Genes')
down_plot
#Save plots as pdf
pdf("3analyses_lgDEL_Upset.pdf", width = 10, height = 8)
all_plot
up_plot
down_plot
dev.off()
#Subset gene list shared between all 3 analyses
genes <- all_upset_df[which(all_upset_df$C_B == 'TRUE' & all_upset_df$ConditionOnly == 'TRUE' & all_upset_df$CplusB == 'TRUE'),]
genesUP <- up_upset_df[which(up_upset_df$C_B == 'TRUE' & up_upset_df$ConditionOnly == 'TRUE' & up_upset_df$CplusB == 'TRUE'),]
genesDOWN <- down_upset_df[which(down_upset_df$C_B == 'TRUE' & down_upset_df$ConditionOnly == 'TRUE' & down_upset_df$CplusB == 'TRUE'),]
#Order lists
genes <- genes[order(genes$Gene),]
genesUP <- genesUP[order(genesUP$Gene),]
genesDOWN <- genesDOWN[order(genesDOWN$Gene),]
#Export GeneID list
write.table(genes$Gene, "allSigGene_3analyses_smDEL.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(genesUP$Gene, "UpGene_3analyses_smDEL.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(genesDOWN$Gene, "DownGene_3analyses_smDEL.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
#Determine external gene name
genesName <- ResGeneID_ALL[genes$Gene,]
genesUPName <- ResGeneID_ALL[genesUP$Gene,]
genesDOWNName <- ResGeneID_ALL[genesDOWN$Gene,]
#Drop unnecessary columns
genesName <- genesName[c(1,43)]
genesUPName <- genesUPName[c(1,43)]
genesDOWNName <- genesDOWNName[c(1,43)]
#Order list
genesName <- genesName[order(genesName$external_gene_name),]
genesUPName <- genesUPName[order(genesUPName$external_gene_name),]
genesDOWNName <- genesDOWNName[order(genesDOWNName$external_gene_name),]
#Export gene name lists
write.table(genesName$external_gene_name, "allSigGeneName_3analyses_smDEL.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(genesUPName$external_gene_name, "UpGeneName_3analyses_smDEL.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(genesDOWNName$external_gene_name, "DownGeneName_3analyses_smDEL.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)

##---Significance Testing---##
#Fishers exact to test for enrichment of imprinted genes in DEG results list
dat <- data.frame(
  "Imprint_no" = c(302,22168),
  "Imprint_yes" = c(5,270),
  row.names = c("DE", "Non-DE"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("Non-imprinted", "Imprinted")
dat
chisq.test(dat)$expected
test <- fisher.test(dat)
test

##Permutation for allSig overlaps between smDEL & lgDEL
#Subset significant genes from results objects
H9sigGenes <- ResGeneID_ALL[which(ResGeneID_ALL$padj_H9 < 0.05),]
CT2sigGenes <- ResGeneID_ALL[which(ResGeneID_ALL$padj_CT2 < 0.05),]
CsigGenes <- ResGeneID_lgDEL_C[which(ResGeneID_lgDEL_C$padj < 0.05),]
CplusBsigGenes <- ResGeneID_lgDEL_CB[which(ResGeneID_lgDEL_CB$padj < 0.05),]
#Pull list of genes
H9sigGeneList <- as.data.frame(H9sigGenes$Gene)
names(H9sigGeneList)[1] <- "Gene"
CT2sigGeneList <- as.data.frame(CT2sigGenes$Gene)
names(CT2sigGeneList)[1] <- "Gene"
CsigGeneList <- as.data.frame(CsigGenes$Gene)
names(CsigGeneList)[1] <- "Gene"
CplusBsigGeneList <- as.data.frame(CplusBsigGenes$Gene)
names(CplusBsigGeneList)[1] <- "Gene"
#Combine lists
lgDELsigGene <- rbind(H9sigGeneList, CT2sigGeneList, CsigGeneList, CplusBsigGeneList)
#Remove replicates
lgDELsigGene <- unique(lgDELsigGene$Gene)
#Save list
write.table(lgDELsigGene, "lgDELsigGene_permtest.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
#Copy text files onto the cluster and run permutation_test_rbg.sbatch (Done with ENSEMBL)
#Read in files from permutation
smtolg <- read.table("../smtolg.txt", header=FALSE)
lgtosm <- read.table("../lgtosm.txt", header=FALSE)
#Calculate median value
median(smtolg$V1) #equal to 14
median(lgtosm$V1) #equal to 11
#Plot histograms
hist_smtolg <- ggplot(smtolg, aes(x=V1)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = 14), col='red', size=.6) +
  geom_vline(aes(xintercept = 42), col='purple', size=1, linetype="dashed") +
  scale_x_continuous(breaks = seq(0, 55, 5)) +
  labs(title="Histogram of smtolg Permutation",x="Number of Overlapping Genes", y = "Frequency")
hist_smtolg
hist_lgtosm <- ggplot(lgtosm, aes(x=V1)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = 11), col='red', size=.6) +
  geom_vline(aes(xintercept = 42), col='purple', size=1, linetype="dashed") +
  scale_x_continuous(breaks = seq(0, 55, 5)) +
  labs(title="Histogram of lgtosm Permutation",x="Number of Overlapping Genes", y = "Frequency")
hist_lgtosm
#Save as pdf
pdf("../PermutationsHistogram.pdf", width = 10, height = 8)
hist_smtolg
hist_lgtosm
dev.off()

##---Scatter Plots---##
#Subset stat for H9 & CT2
stat_df <- as.data.frame(ResGene_ALL_sharedID[,c(5,53,43)])
#Take absolute values of stat
stat_df$abs_stat_H9 <- abs(stat_df$stat_H9)
stat_df$abs_stat_CT2 <- abs(stat_df$stat_CT2)
#Plot
stat_scatter <- ggplot(stat_df, aes(x=abs_stat_H9, y=abs_stat_CT2)) +
  geom_point(color = "black") + geom_smooth(method = lm, colour="black", fill = "lightgray", alpha = 0.25) + stat_cor(aes(label = paste(..p.label..)), label.x = 25, label.y = 30) + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~")), label.x = 25, label.y = 25) +
  ggtitle("Stat H9 vs CT2 lgDEL") + theme_classic() + labs(x = "|H9 stat|", y = "|CT2 stat|")
stat_scatter

absolute value or not??

#Save summary tables
write.csv(as.data.frame(resDataH9lgDEL), file = "H9lgDEL_Results_DESeqSTAR.csv")
write.csv(as.data.frame(resDataCT2lgDEL), file = "CT2lgDEL_Results_DESeqSTAR.csv")

#Combine saved plots using Adobe
