##---Checking DEGs with SNORD116 del mouse---##
library("tidyverse")
library("ggplot2")

#Set working directory
directory <- "../OneDrive - UConn Health Center/NGN_SNORD116_Paper/Submission_NAR/Resubmission/"
setwd(directory)

#Load and format files
hum_mouse_ortho <- read.table(file = "mart_export.txt", sep = '\t', header = TRUE)
colnames(hum_mouse_ortho) <- c("ENSEMBL_ID_human", "ENSEMBL_ID_mouse", "mouse_gene_name", "human_gene_name", "source_hum_name")
mouse_zt6 <- read.csv("PMID23771028_mouse_zt6.csv", skip = 3, header = T)
mouse_zt6 <- mouse_zt6[c(1,4:7)]
colnames(mouse_zt6) <- c("mouse_gene_name", "WT.FPKM.zt6", "Snord116del.FPKM.zt6", "log2.fold_change.zt6", "q_value.zt6")
smDEL_up <- read.delim("../../../Cotney_Lab/PWS_RNASeq/STAR/smDEL/VennDiagram_IntersectUpResultsENSEMBL.txt", sep="\t", head=F)
smDEL_up$direction_hum <- "up"
smDEL_down <- read.delim("../../../Cotney_Lab/PWS_RNASeq/STAR/smDEL/VennDiagram_IntersectDownResultsENSEMBL.txt", sep="\t", head=F)
smDEL_down$direction_hum <- "down"
allsmDEL <- rbind(smDEL_up, smDEL_down)
colnames(allsmDEL)[1] <- "ENSEMBL_ID_human"

#Filter ortholog file for 1:1 orthology
one2one_ortho <- hum_mouse_ortho[!(duplicated(hum_mouse_ortho$ENSEMBL_ID_human) | duplicated(hum_mouse_ortho$ENSEMBL_ID_human, fromLast = TRUE)), ]
one2one_ortho <- one2one_ortho[!(duplicated(one2one_ortho$ENSEMBL_ID_mouse) | duplicated(one2one_ortho$ENSEMBL_ID_mouse, fromLast = TRUE)), ]
  #These are the genes that are not 1:1
  other <- hum_mouse_ortho[(duplicated(hum_mouse_ortho$ENSEMBL_ID_human) | duplicated(hum_mouse_ortho$ENSEMBL_ID_human, fromLast = TRUE)), ]
  other2 <- hum_mouse_ortho[(duplicated(hum_mouse_ortho$ENSEMBL_ID_mouse) | duplicated(hum_mouse_ortho$ENSEMBL_ID_mouse, fromLast = TRUE)), ]

#Merge with ortho list
all_tested_mouse <- merge(one2one_ortho, mouse_zt6, by="mouse_gene_name", sort = TRUE, all = FALSE)
  #Convert to lists of ENSEMBL genes and transfer to cluster
  all_tested_mouse_list <- all_tested_mouse$ENSEMBL_ID_human
  write.table(all_tested_mouse_list, file = "all_tested_mouse_ortho.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#Check how many human DEGs have ortholog
smDEL_ortho <- merge(one2one_ortho, allsmDEL, by="ENSEMBL_ID_human", sort = TRUE, all = FALSE)
  #287 human genes of the 317 total smDEL DEGs have one mouse ortholog
#Convert to lists of ENSEMBL genes and transfer to cluster
smDEL_ortho_list <- smDEL_ortho$ENSEMBL_ID_human
  write.table(smDEL_ortho_list, file = "smDEL_ortho_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Filter for significant values
mouse_zt6_sig <- all_tested_mouse[which(all_tested_mouse$q_value.zt6 < 0.05),]
  #5,348 mouse genes out of the 6,467 total mouse DEGs have one mouse ortholog

#Check for overlaps with human genes
shared_mouse2human <- merge(mouse_zt6_sig, allsmDEL, by="ENSEMBL_ID_human", sort = TRUE, all = FALSE)
  #116 genes are shared DEGs between human and mouse (disregarding directionality)

#Filter for shared directionality
shared_up <- shared_mouse2human[which(shared_mouse2human$log2.fold_change.zt6 > 0 & shared_mouse2human$direction_hum == "up"),]
shared_down <- shared_mouse2human[which(shared_mouse2human$log2.fold_change.zt6 < 0 & shared_mouse2human$direction_hum == "down"),]

#Plot histogram from permutation test
permtest <- as.data.frame(read.table("permtest.txt", sep="\t", head=F))
hist <- ggplot(permtest, aes(x=V1)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = median(V1)), col='green', linewidth=1, linetype="dashed") +
  geom_vline(aes(xintercept = 116), col='purple', linewidth=1) +
  labs(x="Human-Mouse Orthology Overlap", y = "Frequency")
hist
#Figure out p-value
greater <- permtest[which(permtest$V1 > 116),]
equal <- permtest[which(permtest$V1 == "116"),]
  #22 values equal to or greater than experimental value out of 1000 permutations

#Save as pdf
pdf("permtest_hum2mouse.pdf", width = 10, height = 8)
hist
dev.off()

#Import results data from smDEL to add to orthology table
smDEL_res <- read.csv("../../../Cotney_Lab/PWS_RNASeq/STAR/smDEL/CombinedH9CT2_Results_DESeqSTAR.csv", row.names = "X", header = TRUE)
#Drop columns & reorder
smDEL_res_new <- smDEL_res[c(1,3,7,48,52,41)]
colnames(smDEL_res_new)[1] <- "ENSEMBL_ID_human"
#Merge with ortho table
ortho_table_all <- merge(shared_mouse2human, smDEL_res_new, by="ENSEMBL_ID_human", sort = TRUE, all = FALSE)
#Add mouse direction
ortho_table_all$direction_mouse <- "up"
ortho_table_all <- ortho_table_all %>% mutate(direction_mouse = ifelse(log2.fold_change.zt6 < 0, "down", direction_mouse))
#Reorder & save orthology list
ortho_table_all_FIN <- ortho_table_all[c(3,2,1,4,15,16,10,8,9,11:14)]
colnames(ortho_table_all_FIN) <- c("ENSEMBL_ID_mouse", "mouse_gene_name", "ENSEMBL_ID_human", "human_gene_name", "description_human", "direction_mouse", "direction_hum", "log2FC_Zt6Snord116del_mouse", "qvalue_Zt6Snord116del_mouse", "log2FC_H9smDEL_hum", "padj_H9smDEL_hum", "log2FC_CT2smDEL_hum", "padj_smDEL_hum")
  write.csv(ortho_table_all_FIN, file = "human2mouse_sharedGenes.csv")
  