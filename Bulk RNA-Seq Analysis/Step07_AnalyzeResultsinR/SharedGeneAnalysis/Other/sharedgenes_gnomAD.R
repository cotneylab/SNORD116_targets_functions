##---Parsing gnomAD database for set of 42 genes---##
library("biomaRt")
library("dplyr")
library("ggplot2")
library("ggpubr")

#Set working directory
directory <- "../Cotney_Lab/PWS_RNASeq/STAR/"
setwd(directory)

#Import files
gnomAD <- read.delim("gnomad.v2.1.1.lof_metrics.by_gene.txt", sep="\t", head=T)
shared_genes <- read.delim("ENSEMBLsharedsmlgDEL.txt", sep="\t", head=F)

#Get gene names for GENCODEv19/Ensembl74_12.2013 (what gnomADv2.1.1 uses) - closest available release used below
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://feb2014.archive.ensembl.org")
genes.table <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_id"), values = shared_genes$V1, mart = mart)

#Make list of gene names
shared_genes_list <- as.data.frame(genes.table$external_gene_id)

#Rename column to match name of gene column in gnomAD df
colnames(shared_genes_list)[1] <- "gene"

#Merge lists to parse gnomAD
sharedGenes_gnomAD <- merge(shared_genes_list, gnomAD, by = 'gene', all = FALSE)

#Pull out relevant columns
sharedGenes_gnomAD_info <- sharedGenes_gnomAD[c(64, 1:2, 30, 21, 20, 17)]

#Change column names
colnames(sharedGenes_gnomAD_info) <- c("ensembl_gene_id","external_gene_id", "ensembl_transcript_id", "LOEUF", "pLI", "expected_SNVs_pLoF", "observed_SNVs_pLoF")

  #Save file
  write.csv(sharedGenes_gnomAD_info, file = "sharedgenes_gnomADv2.1.1_infotable.csv")

#Filter out the shared genes from rest of gnomAD
background_data <- gnomAD %>% filter(!gnomAD$gene %in% shared_genes_list$gene)

#Add column to dataframes to specify set
sharedGenes_gnomAD$dataset <- "sharedGenes"
background_data$dataset <- "background"

#Rebind lists
all_gnomAD <- rbind(sharedGenes_gnomAD, background_data)

#Check distribution of data
normality_check_LOEUF <- ggplot(all_gnomAD, aes(x=oe_lof_upper)) + geom_histogram(binwidth=0.1, color="black", fill="white")
normality_check_LOEUF
normality_check_pLI <- ggplot(all_gnomAD, aes(x=pLI)) + geom_histogram(binwidth=0.1, color="black", fill="white")
normality_check_pLI
  #Neither measurement is normally distributed, run Wilcoxen test

#Plot shared genes vs background
vi_plot_LOEUF <- ggplot(all_gnomAD, aes(dataset, oe_lof_upper, fill=dataset)) + geom_violin(trim=FALSE, show.legend = FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Distribution of LOEUF Scores", x = "" , y = "LOEUF Score") + theme_classic() + stat_compare_means(comparisons = list(c("background", "sharedGenes")), method = "wilcox.test", label = "p.format")
vi_plot_LOEUF
vi_plot_pLI <- ggplot(all_gnomAD, aes(dataset, pLI, fill=dataset)) + geom_violin(trim=FALSE, show.legend = FALSE) + geom_boxplot(width=0.02, fill="white") + labs(title="Distribution of pLI Scores", x = "" , y = "pLI Score") + theme_classic() + stat_compare_means(comparisons = list(c("background", "sharedGenes")), method = "wilcox.test", label = "p.format")
vi_plot_pLI

#Save as pdfs
pdf("violin_plots.pdf", width = 10, height = 8)
vi_plot_LOEUF
vi_plot_pLI
dev.off()

#Double check statistics manually
u_res_LOEUF <- wilcox.test(oe_lof_upper ~ dataset, data = all_gnomAD)
u_res_LOEUF
effect_size_LOEUF <- wilcox_effsize(oe_lof_upper ~ dataset, data = all_gnomAD)
effect_size_LOEUF
u_res_pLI <- wilcox.test(pLI ~ dataset, data = all_gnomAD)
u_res_pLI
effect_size_pLI <- wilcox_effsize(pLI ~ dataset, data = all_gnomAD)
effect_size_pLI

