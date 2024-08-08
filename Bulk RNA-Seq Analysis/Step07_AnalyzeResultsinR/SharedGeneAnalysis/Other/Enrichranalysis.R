##---Enrichr---##
#Option 1: Using list of shared genes, query Enrichr database
library("reticulate")
  #Install
gget <- import("gget")
  #Read in file of ENSEMBL ID's
shared <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/ENSEMBLsharedsmlgDEL.txt", sep="\t", head=F)
  #Make it into a list
genes <- shared$V1
  #Run Enrichr
df_shared <- gget$enrichr(genes, database = "Enrichr_Submissions_TF-Gene_Coocurrence", ensembl = TRUE)
df_shared_disgenet <- gget$enrichr(genes, database = "DisGeNET", ensembl = TRUE)
df_shared_grants <- gget$enrichr(genes, database = "Genes_Associated_with_NIH_Grants", ensembl = TRUE)
  #Calculate -log10(adjusted p-value)
df_shared_calc <- df_shared
df_shared_calc$neglogpadj <- -log10(df_shared_calc$adj_p_val)
  #Subset my TFs
resTF <- df_shared_calc[df_shared_calc$path_name %in% c('IRX5', 'MYBL2', 'NR5A2', 'PAX5', 'PAX6', 'PLAGL1', 'SOX21', 'ZIC2'),]
  #Count how many of my TFs are significant
test.stat1 <- nrow(resTF[resTF$adj_p_val < .05, ])
  #Save full table
  df_shared_calc <- apply(df_shared_calc,2,as.character)
  write.csv(as.data.frame(df_shared_calc), file = "../Cotney_Lab/PWS_RNASeq/STAR/Enrichr_Query_sharedGenes.csv")

#Option 2: Using list of shared genes, query Enrichr database
library("enrichR")
  #Read in file of ENSEMBL ID's
shared <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/ENSEMBLsharedsmlgDEL.txt", sep="\t", head=F)
  #Change ESEMBL ID's to HGNC symbols
mart <- useMart("ensembl", host = "https://dec2017.archive.ensembl.org/", dataset = "hsapiens_gene_ensembl")
Rs <- data.table(shared)
colnames(Rs) <- c("Gene")
genes.table_shared <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values= Rs$Gene, mart= mart)
  #Make into list (only ones with HGNC symbol)
genes.table_shared[genes.table_shared==""] <- NA
gene_set <- na.omit(genes.table_shared)
###This gets rid of 2 genes
gene_set <- gene_set$hgnc_symbol
  #Export list
write.table(gene_set, "../Cotney_Lab/PWS_RNASeq/STAR/sharedsmlgDEL_enrichr.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
  #Run Enrichr
enriched_shared <- as.data.frame(enrichr(gene_set, databases = "Enrichr_Submissions_TF-Gene_Coocurrence"))
  #Subset my TFs
enriched_TFshared <- enriched_shared[enriched_shared$Enrichr_Submissions_TF.Gene_Coocurrence.Term %in% c('IRX5', 'MYBL2', 'NR5A2', 'PAX5', 'PAX6', 'PLAGL1', 'SOX21', 'ZIC2'),]
  #Count how many of my TFs are significant
test.stat1 <- nrow(enriched_TFshared[enriched_TFshared$Enrichr_Submissions_TF.Gene_Coocurrence.Adjusted.P.value < .05, ])
  #Subset for general TFs
genTF_shared <- enriched_shared[enriched_shared$Enrichr_Submissions_TF.Gene_Coocurrence.Term %in% c('MYC', 'JUN', 'JUNB', 'FOS', 'FOSB', 'BDP1', 'GTF2F1', 'GTF2F2'),]
  #Count how many general TFs are significant
test.stat2 <- nrow(genTF_shared[genTF_shared$Enrichr_Submissions_TF.Gene_Coocurrence.Adjusted.P.value < .05, ])

#Run 100 permutations to test significance
library("biomaRt")
library("data.table")
library("purrr")
  #Read in Non-DEG list made for HOMER analysis
nonDEG <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/nonDEGall.txt", sep="\t", head=F)
  #Change ESEMBL ID's to HGNC symbols
mart <- useMart("ensembl", host = "https://dec2017.archive.ensembl.org/", dataset = "hsapiens_gene_ensembl")
R <- data.table(nonDEG)
colnames(R) <- c("Gene")
genes.table_nonDEG <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values= R$Gene, mart= mart)
  #Make into list (only ones with HGNC symbol)
genes.table_nonDEG[genes.table_nonDEG==""] <- NA
nonDEG <- na.omit(genes.table_nonDEG)
###This gets rid of 5302 genes
nonDEG <- nonDEG$hgnc_symbol
  #Set seed for reproducibility of results
set.seed(40)
  #Set number of observations to sample (in this case = number of genes)
n <- 40
  #Set the number of permutations
P <- 1000
  #Make lists
random_gene_lists <- map(1:P, ~ sample(nonDEG, 40, replace = FALSE))
  #Create empty variable
results <- list()
  #Run Enrichr
for (i in 1:length(random_gene_lists))
{
  gene_random <- random_gene_lists[[i]]
  permutations <- as.data.frame(enrichr(gene_random, databases = "Enrichr_Submissions_TF-Gene_Coocurrence"))
  results[[i]] <- permutations
}
  #Subset my TFs & general TFs
resTFperm <- list()
genTFperm <- list()
for (i in 1:length(random_gene_lists))
{
  myTF_perm <- results[[i]][results[[i]]$Enrichr_Submissions_TF.Gene_Coocurrence.Term %in% c('IRX5', 'MYBL2', 'NR5A2', 'PAX5', 'PAX6', 'PLAGL1', 'SOX21', 'ZIC2'),]
  resTFperm[[i]] <- myTF_perm
  genTF_perm <- results[[i]][results[[i]]$Enrichr_Submissions_TF.Gene_Coocurrence.Term %in% c('MYC', 'JUN', 'JUNB', 'FOS', 'FOSB', 'BDP1', 'GTF2F1', 'GTF2F2'),]
  genTFperm[[i]] <- genTF_perm
}
  #Initialize vectors to store all of the test-stats
Perm.test.stat1 <- rep(0, P)
Perm.test.stat2 <- rep(0, P)
  #Count how many of TFs are significant
for (i in 1:length(random_gene_lists))
{
  Perm.test.stat1[[i]] <- nrow(resTFperm[[i]][resTFperm[[i]]$Enrichr_Submissions_TF.Gene_Coocurrence.Adjusted.P.value < .05, ])
  Perm.test.stat2[[i]] <- nrow(genTFperm[[i]][genTFperm[[i]]$Enrichr_Submissions_TF.Gene_Coocurrence.Adjusted.P.value < .05, ])
}

##Testing Hypothesis: TFs that are not expressed in neurons should not be enriched in our set of 42 genes.
  #Read in counts table
counts <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/lgDEL/PWSRNAseq_STAR_featureCounts.txt", sep="\t", head=T, skip = 1, row.names = "Geneid")
  #Simplify column names
colnames(counts) <- c("chr", "start", "end", "strand", "length", "CT2_1", "CT2_2", "CT2_3", "CT2_4", "CT2_5", "CT2_6", "CT2_lgDEL_1", "CT2_lgDEL_2", "CT2_lgDEL_3", "CT2_lgDEL_4", "CT2_lgDEL_5", "CT2_lgDEL_6", "CT2_smDEL_1", "CT2_smDEL_2", "CT2_smDEL_3", "CT2_smDEL_4", "CT2_smDEL_5", "CT2_smDEL_6", "H9_1", "H9_2", "H9_3", "H9_4", "H9_5", "H9_6", "H9_lgDEL_1", "H9_lgDEL_2", "H9_lgDEL_3", "H9_lgDEL_4", "H9_lgDEL_5", "H9_lgDEL_6", "H9_smDEL_1", "H9_smDEL_2", "H9_smDEL_3", "H9_smDEL_4", "H9_smDEL_5")
  #Get rid of first 5 columns
counts <- counts[ ,6:ncol(counts)]
  #Create a column for sum
counts$rowsum <- rowSums(counts)
  #Subset unexpressed genes
no_counts <- counts[which(counts$rowsum == "0"),]
no_counts <- as.data.frame(rownames(no_counts))
colnames(no_counts) <- c("Ensembl.ID")
  #Remove decimal in gene name
no_counts$Ensembl.ID <- gsub('\\..+$', '', no_counts$Ensembl.ID)
  #Import file of all TFs: official list of human TFs from Lambert et al (PMID:29425488)
TFall <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/DatabaseExtract_v_1.01.txt", sep="\t", head=T)
  #Delete first column
TFall <- TFall[-c(1)]
  #Take only proteins thought to be TFs
TFall <- TFall[which(TFall$Is.TF. == "Yes"),]
  #Order list of all TFs & parse by no count list
TFall <- TFall[order(TFall$Ensembl.ID),]
unexpressedTF <- merge(TFall, no_counts, by="Ensembl.ID", sort=TRUE, all = FALSE)
  #Make into list
unexpressedTF_list <- unexpressedTF$HGNC.symbol
  #Subset unexpressed TFs from results
res_unexpTF <- df_shared[df_shared$path_name %in% c(unexpressedTF_list),]
  #Export list of unexpressed TFs
write.table(unexpressedTF_list, "../Cotney_Lab/PWS_RNASeq/STAR/unexpressedTFs.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
  ##Hypothesis is not supported.


##UNFINISHED##
#Run 100 permutations to test significance
library("biomaRt")
library("data.table")
  #Read in Non-DEG list made for HOMER analysis
nonDEG <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/nonDEGall.txt", sep="\t", head=F)
  #Change ESEMBL ID's to gene names
mart <- useMart("ensembl", host = "https://dec2017.archive.ensembl.org/", dataset = "hsapiens_gene_ensembl")
R <- data.table(nonDEG)
colnames(R) <- c("Gene")
genes.table_nonDEG <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"), values= R$Gene, mart= mart)
  #Make into list
nonDEG <- genes.table_nonDEG$external_gene_name
  #Set seed for reproducability of results
set.seed(40)
  #Set number of observations to sample (in this case = number of genes)
n <- 40
  #Set the number of permutations
P <- 100
  #Set the variable to resample from
variable <- nonDEG
  #Initialize a matrix to store the permutation data
PermSamples <- matrix(0, nrow = n, ncol = P)
  #Use loop to get permutation samples
for(i in 1:P)
{
  PermSamples[, i] <- sample(variable,
                             size = n,
                             replace = FALSE)
}
  #Convert to data frame
PermSamples_df <- as.data.frame(PermSamples)
  #Initialize vectors to store all of the test-stats
Perm.test.stat1 <- rep(0, P)
  #Create empty objects
genes_random <- vector("list", n)
df_random <- df_shared[FALSE,]
resTF_random <- resTF[FALSE,]
  #Loop through and calculate the test-stats
for (i in 1:P)
{
  #Make into a list
  genes_random[[i]] <- PermSamples_df$V[[i]]
  #Run Enrichr
  df_random[[i]] <- gget$enrichr(genes_random[[i]], database = "Enrichr_Submissions_TF-Gene_Coocurrence")
  #Subset my TFs
  resTF_random[[i]] <- df_random[[i]][df_random[[i]]$path_name %in% c('IRX5', 'MYBL2', 'NR5A2', 'PAX5', 'PAX6', 'PLAGL1', 'SOX21', 'ZIC2'),]
  #Calculate the perm-test-stat1 and save it
  Perm.test.stat1[[i]] <- nrow(resTF_random[[i]][resTF_random1[[i]]$adj_p_val < .05, ])
}
