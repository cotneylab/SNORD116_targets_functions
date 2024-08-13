##---snoGloBe analysis---##
#Import files
targets115 <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/snord115_degshared", sep="\t", head=T))
targets116 <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/snord116_degshared", sep="\t", head=T))
targets_other <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/chr15otherSNORDs_degshared", sep="\t", head=T))
targets113 <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/snord113_degshared", sep="\t", head=T))
targets114 <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/snord114_degshared", sep="\t", head=T))

##Creating bed files for upload to UCSC browser
#Multiply count column by 100 to generate a "score" for the bed file
bed115 <- targets115name
bed115$score <- bed115$count * 100
bed116 <- targets116name
bed116$score <- bed116$count * 100
bedchr15other <- targets_other2
bedchr15other$score <- bedchr15other$count * 100
#Drop unnecessary columns & reorder
bed115 <- bed115[c(1:3,8,14,6)]
bed116 <- bed116[c(1:3,8,14,6)]
bedchr15other <- bedchr15other[c(1:3,8,14,6)]
#Write out file
write.table(bed115, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/snord115targets_snoglobe_hg38.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(bed116, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/snord116targets_snoglobe_hg38.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(bedchr15other, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/otherchr15SNORDstargets_snoglobe_hg38.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  #In notepad, edited file to include trackline & uploaded into browser if desired (https://genome.ucsc.edu/s/rbgilmore/snoGloBe_targets_RBG)

#Analysis of results
#Change ESEMBL ID's to external gene name
library("biomaRt")
library("data.table")
mart <- useMart("ensembl", host = "https://apr2018.archive.ensembl.org/", dataset = "hsapiens_gene_ensembl")
  #SNORD115
R115 <- data.table(targets115$target_id)
colnames(R115) <- c("target_id")
genes.table_115 <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"), values= R115$target_id, mart= mart)
names(genes.table_115)[1] <- "target_id"
  #SNORD116
R116 <- data.table(targets116$target_id)
colnames(R116) <- c("target_id")
genes.table_116 <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"), values= R116$target_id, mart= mart)
names(genes.table_116)[1] <- "target_id"
  #other chr15
R_other <- data.table(targets_other$target_id)
colnames(R_other) <- c("target_id")
genes.table_other <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"), values= R_other$target_id, mart= mart)
names(genes.table_other)[1] <- "target_id"
  #SNORD113
R113 <- data.table(targets113$target_id)
colnames(R113) <- c("target_id")
genes.table_113 <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"), values= R113$target_id, mart= mart)
names(genes.table_113)[1] <- "target_id"
  #SNORD114
R114 <- data.table(targets114$target_id)
colnames(R114) <- c("target_id")
genes.table_114 <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"), values= R114$target_id, mart= mart)
names(genes.table_114)[1] <- "target_id"

##Parse by gene target & plot
library("ggplot2")
library("dplyr")
  #SNORD115
genecount115 <- as.data.frame(targets115 %>% count(target_id))
genecount115name <- merge(genes.table_115, genecount115, by="target_id", sort=TRUE, all = TRUE)
geneplot115 <- ggplot(data = genecount115name, aes(x=external_gene_name, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD115: Shared Gene Targeting Events", x ="Target Gene Name", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5)
geneplot115
  #SNORD116
genecount116 <- as.data.frame(targets116 %>% count(target_id))
genecount116name <- merge(genes.table_116, genecount116, by="target_id", sort=TRUE, all = TRUE)
geneplot116 <- ggplot(data = genecount116name, aes(x=external_gene_name, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD116: Shared Gene Targeting Events", x ="Target Gene Name", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5)
geneplot116
  #other chr15
genecountchr15 <- as.data.frame(targets_other %>% count(target_id))
genecountchr15name <- merge(genes.table_other, genecountchr15, by="target_id", sort=TRUE, all = TRUE)
geneplotchr15 <- ggplot(data = genecountchr15name, aes(x=external_gene_name, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORDchr15: Shared Gene Targeting Events", x ="Target Gene Name", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.6, color="white", size=3.5)
geneplotchr15
  #SNORD113
genecount113 <- as.data.frame(targets113 %>% count(target_id))
genecount113name <- merge(genes.table_113, genecount113, by="target_id", sort=TRUE, all = TRUE)
geneplot113 <- ggplot(data = genecount113name, aes(x=external_gene_name, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD113: Shared Gene Targeting Events", x ="Target Gene Name", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5)
geneplot113
  #SNORD114
genecount114 <- as.data.frame(targets114 %>% count(target_id))
genecount114name <- merge(genes.table_114, genecount114, by="target_id", sort=TRUE, all = TRUE)
geneplot114 <- ggplot(data = genecount114name, aes(x=external_gene_name, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD114: Shared Gene Targeting Events", x ="Target Gene Name", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5)
geneplot114
#Save controls as pdf
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/GenesTargeted_controls.pdf", width = 10, height = 8)
geneplot113
geneplot114
dev.off()

##Remove SNORD copies outside of chr15 region & re-plot
library("tidyr")
  #SNORD115
targets115name <- as.data.frame(targets115 %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
targets115name$snoRNA_number <- sapply(targets115name$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
chr15targets115 <- targets115name %>% filter(targets115name$snoRNA_number != 'SNORD115')
chr15gene115 <- as.data.frame(chr15targets115 %>% count(target_id))
chr15gene115name <- merge(genes.table_115, chr15gene115, by="target_id", sort=TRUE, all = TRUE)
chr15geneplot115 <- ggplot(data = chr15gene115name, aes(x=external_gene_name, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="chr15-SNORD115: Shared Gene Targeting Events", x ="Target Gene Name", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5)
chr15geneplot115
  #SNORD116
targets116name <- as.data.frame(targets116 %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
targets116name$snoRNA_number <- sapply(targets116name$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
chr15targets116 <- targets116name %>% filter(targets116name$snoRNA_number != 'SNORD116')
chr15gene116 <- as.data.frame(chr15targets116 %>% count(target_id))
chr15gene116name <- merge(genes.table_116, chr15gene116, by="target_id", sort=TRUE, all = FALSE)
chr15geneplot116 <- ggplot(data = chr15gene116name, aes(x=external_gene_name, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="chr15-SNORD116: Shared Gene Targeting Events", x ="Target Gene Name", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5)
chr15geneplot116
#Save as pdf
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/GenesTargetedchr15.pdf", width = 10, height = 8)
geneplot115
geneplot116
chr15geneplot115
chr15geneplot116
geneplotchr15
dev.off()

##Parse by sno copy & plot
  #SNORD115
snocount115 <- as.data.frame(targets115 %>% count(sno_id))
snocount115name <- as.data.frame(snocount115 %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
snocount115name$snoRNA_number <- sapply(snocount115name$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
snocount115name[43, "snoRNA_number"] <- "SNORD115_nonchr15_1"
snocount115name[44, "snoRNA_number"] <- "SNORD115_nonchr15_2"
snocount115name <- as.data.frame(snocount115name %>% separate(snoRNA_number, c("SNORD115", "SNORD115_number"), sep = "-"))
snocount115name <- snocount115name %>% mutate(SNORD115_number = as.numeric(gsub("[^0-9]", "", SNORD115_number))) %>% arrange(SNORD115_number)
snocount115name$snoRNA_number <- paste(snocount115name$SNORD115, snocount115name$SNORD115_number, sep = "-")
snocount115name[49, "snoRNA_number"] <- "SNORD115_nonchr15_1"
snocount115name[50, "snoRNA_number"] <- "SNORD115_nonchr15_2"
snocount115name$snoRNA_number <- factor(snocount115name$snoRNA_number, levels = snocount115name$snoRNA_number)
snoplot115 <- ggplot(data = snocount115name, aes(x=snoRNA_number, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD115: Targeting Events by snoRNA Copy", x ="snoRNA Copy", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5) + scale_x_discrete(limits = levels(snocount115name$snoRNA_number))
snoplot115
  #SNORD116
snocount116 <- as.data.frame(targets116 %>% count(sno_id))
snocount116name <- as.data.frame(snocount116 %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
snocount116name$snoRNA_number <- sapply(snocount116name$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
snocount116name[2, "snoRNA_number"] <- "SNORD116_nonchr15_1"
snocount116name[23, "snoRNA_number"] <- "SNORD116_nonchr15_2"
snocount116name[28, "snoRNA_number"] <- "SNORD116_nonchr15_3"
snocount116name <- as.data.frame(snocount116name %>% separate(snoRNA_number, c("SNORD116", "SNORD116_number"), sep = "-"))
snocount116name <- snocount116name %>% mutate(SNORD116_number = as.numeric(gsub("[^0-9]", "", SNORD116_number))) %>% arrange(SNORD116_number)
snocount116name$snoRNA_number <- paste(snocount116name$SNORD116, snocount116name$SNORD116_number, sep = "-")
snocount116name[31, "snoRNA_number"] <- "SNORD116_nonchr15_1"
snocount116name[32, "snoRNA_number"] <- "SNORD116_nonchr15_2"
snocount116name[33, "snoRNA_number"] <- "SNORD116_nonchr15_3"
snocount116name$snoRNA_number <- factor(snocount116name$snoRNA_number, levels = snocount116name$snoRNA_number)
snoplot116 <- ggplot(data = snocount116name, aes(x=snoRNA_number, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD116: Targeting Events by snoRNA Copy", x ="snoRNA Copy", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5) + scale_x_discrete(limits = levels(snocount116name$snoRNA_number))
snoplot116
  #chr15 other
snocountchr15 <- as.data.frame(targets_other %>% count(sno_id))
snocountchr15name <- as.data.frame(snocountchr15 %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
snocountchr15name$snoRNA_number <- sapply(snocountchr15name$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
chr15_order <- c("SNORD107", "SNORD64", "SNORD108", "SNORD109A", "SNORD109B")
snocountchr15name <- snocountchr15name[order(match(snocountchr15name$snoRNA_number, chr15_order)), ]
snocountchr15name$snoRNA_number <- factor(snocountchr15name$snoRNA_number, levels = snocountchr15name$snoRNA_number)
snoplotchr15 <- ggplot(data = snocountchr15name, aes(x=snoRNA_number, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORDchr15: Targeting Events by snoRNA Copy", x ="snoRNA Copy", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5) + scale_x_discrete(limits = levels(snocountchr15name$snoRNA_number))
snoplotchr15
  #SNORD113
snocount113 <- as.data.frame(targets113 %>% count(sno_id))
snocount113name <- as.data.frame(snocount113 %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
snocount113name$snoRNA_number <- sapply(snocount113name$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
snocount113name <- as.data.frame(snocount113name %>% separate(snoRNA_number, c("SNORD113", "SNORD113_number"), sep = "-"))
snocount113name <- snocount113name %>% mutate(SNORD113_number = as.numeric(gsub("[^0-9]", "", SNORD113_number))) %>% arrange(SNORD113_number)
snocount113name$SNORD113_number[is.na(snocount113name$SNORD113_number)] <- snocount113name$snoENSEMBL_ID[is.na(snocount113name$SNORD113_number)]
snocount113name$snoRNA_number <- paste(snocount113name$SNORD113, snocount113name$SNORD113_number, sep = "-")
snocount113name$snoRNA_number <- factor(snocount113name$snoRNA_number, levels = snocount113name$snoRNA_number)
snoplot113 <- ggplot(data = snocount113name, aes(x=snoRNA_number, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD113: Targeting Events by snoRNA Copy", x ="snoRNA Copy", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5) + scale_x_discrete(limits = levels(snocount113name$snoRNA_number))
snoplot113
  #SNORD114
snocount114 <- as.data.frame(targets114 %>% count(sno_id))
snocount114name <- as.data.frame(snocount114 %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
snocount114name$snoRNA_number <- sapply(snocount114name$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
snocount114name <- as.data.frame(snocount114name %>% separate(snoRNA_number, c("SNORD114", "SNORD114_number"), sep = "-"))
snocount114name <- snocount114name %>% mutate(SNORD114_number = as.numeric(gsub("[^0-9]", "", SNORD114_number))) %>% arrange(SNORD114_number)
snocount114name$snoRNA_number <- paste(snocount114name$SNORD114, snocount114name$SNORD114_number, sep = "-")
snocount114name$snoRNA_number <- factor(snocount114name$snoRNA_number, levels = snocount114name$snoRNA_number)
snoplot114 <- ggplot(data = snocount114name, aes(x=snoRNA_number, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD114: Targeting Events by snoRNA Copy", x ="snoRNA Copy", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5) + scale_x_discrete(limits = levels(snocount114name$snoRNA_number))
snoplot114
#Save controls as pdf
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/TargetsbysnoRNAcopyControls.pdf", width = 10, height = 8)
snoplot113
snoplot114
dev.off()

##Remove SNORD copies outside of chr15 region & plot
  #SNORD115
chr15sno115 <- snocount115name %>% filter(!snocount115name$snoRNA_number %in% c("SNORD115_nonchr15_1", "SNORD115_nonchr15_2"))
chr15sno115$snoRNA_number <- factor(chr15sno115$snoRNA_number, levels = chr15sno115$snoRNA_number)
chr15snoplot115 <- ggplot(data = chr15sno115, aes(x=snoRNA_number, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="chr15-SNORD115: Targeting Events by snoRNA Copy", x ="snoRNA Copy", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5) + scale_x_discrete(limits = levels(chr15sno115$snoRNA_number))
chr15snoplot115
  #SNORD116
chr15sno116 <- snocount116name %>% filter(!snocount116name$snoRNA_number %in% c("SNORD116_nonchr15_1", "SNORD116_nonchr15_2", "SNORD116_nonchr15_3"))
chr15sno116$snoRNA_number <- factor(chr15sno116$snoRNA_number, levels = chr15sno116$snoRNA_number)
chr15snoplot116 <- ggplot(data = chr15sno116, aes(x=snoRNA_number, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="chr15-SNORD116: Targeting Events by snoRNA Copy", x ="snoRNA Copy", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5) + scale_x_discrete(limits = levels(chr15sno116$snoRNA_number))
chr15snoplot116
#Save as pdf
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/TargetsbysnoRNAcopychr15.pdf", width = 10, height = 8)
snoplot115
snoplot116
chr15snoplot115
chr15snoplot116
snoplotchr15
dev.off()
  #For SNORD116 figure
fig116 <- ggplot(data = chr15sno116, aes(x=snoRNA_number, y=n)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=7), panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "lightgray"), panel.background = element_rect(fill = "white"), axis.line = element_line(color="black")) + 
  labs(title="chr15-SNORD116: Targeting Events by snoRNA Copy", x ="snoRNA Copy", y = "Number of Targeting Events") + 
  geom_text(aes(label=n), vjust=1.1, color="white", size=3.5) + scale_x_discrete(limits = levels(chr15sno116$snoRNA_number))
fig116
  #Save as pdf
  pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116_targetsbycopy_fig.pdf", width = 10, height = 8)
  fig116
  dev.off()

##Subset by sno copy & gene
  #SNORD115
grouped_targets115 <- chr15targets115 %>% group_by(target_id, snoRNA_number)
summary_targets115 <- grouped_targets115 %>% summarize(count = n())
summary_targets115name <- merge(genes.table_115, summary_targets115, by="target_id", sort=TRUE, all = FALSE)
summary_targets115name <- as.data.frame(summary_targets115name %>% separate(snoRNA_number, c("SNORD115", "SNORD115_number"), sep = "-"))
summary_targets115name <- summary_targets115name %>% mutate(SNORD115_number = as.numeric(gsub("[^0-9]", "", SNORD115_number))) %>% arrange(SNORD115_number)
summary_targets115name$snoRNA_number <- paste(summary_targets115name$SNORD115, summary_targets115name$SNORD115_number, sep = "-")
  #SNORD116
grouped_targets116 <- chr15targets116 %>% group_by(target_id, snoRNA_number)
summary_targets116 <- grouped_targets116 %>% summarize(count = n())
summary_targets116name <- merge(genes.table_116, summary_targets116, by="target_id", sort=TRUE, all = FALSE)
summary_targets116name <- as.data.frame(summary_targets116name %>% separate(snoRNA_number, c("SNORD116", "SNORD116_number"), sep = "-"))
summary_targets116name <- summary_targets116name %>% mutate(SNORD116_number = as.numeric(gsub("[^0-9]", "", SNORD116_number))) %>% arrange(SNORD116_number)
summary_targets116name$snoRNA_number <- paste(summary_targets116name$SNORD116, summary_targets116name$SNORD116_number, sep = "-")
  #chr15 other
targets_other2 <- as.data.frame(targets_other %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
targets_other2$snoRNA_number <- sapply(targets_other2$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
grouped_targetschr15 <- targets_other2 %>% group_by(target_id, snoRNA_number)
summary_targetschr15 <- grouped_targetschr15 %>% summarize(count = n())
summary_targetschr15name <- merge(genes.table_other, summary_targetschr15, by="target_id", sort=TRUE, all = FALSE)
  #SNORD113
targets113_2 <- as.data.frame(targets113 %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
targets113_2$snoRNA_number <- sapply(targets113_2$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
targets113_2 <- as.data.frame(targets113_2 %>% separate(snoRNA_number, c("SNORD113", "SNORD113_number"), sep = "-"))
targets113_2$SNORD113_number[is.na(targets113_2$SNORD113_number)] <- targets113_2$snoENSEMBL_ID[is.na(targets113_2$SNORD113_number)]
targets113_2 <- targets113_2 %>% mutate(SNORD113_number = ifelse(SNORD113_number %in% NA, targets113_2$snoENSEMBL_ID, SNORD113_number))
targets113_2$snoRNA_number <- paste(targets113_2$SNORD113, targets113_2$SNORD113_number, sep = "-")
grouped_targets113 <- targets113_2 %>% group_by(target_id, snoRNA_number)
summary_targets113 <- grouped_targets113 %>% summarize(count = n())
summary_targets113name <- merge(genes.table_113, summary_targets113, by="target_id", sort=TRUE, all = FALSE)
summary_targets113name <- as.data.frame(summary_targets113name %>% separate(snoRNA_number, c("SNORD113", "SNORD113_number"), sep = "-"))
summary_targets113name <- summary_targets113name %>% mutate(SNORD113_number = as.numeric(gsub("[^0-9]", "", SNORD113_number))) %>% arrange(SNORD113_number)
summary_targets113name$snoRNA_number <- paste(summary_targets113name$SNORD113, summary_targets113name$SNORD113_number, sep = "-")
  #SNORD114
targets114_2 <- as.data.frame(targets114 %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
targets114_2$snoRNA_number <- sapply(targets114_2$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
grouped_targets114 <- targets114_2 %>% group_by(target_id, snoRNA_number)
summary_targets114 <- grouped_targets114 %>% summarize(count = n())
summary_targets114name <- merge(genes.table_114, summary_targets114, by="target_id", sort=TRUE, all = FALSE)
summary_targets114name <- as.data.frame(summary_targets114name %>% separate(snoRNA_number, c("SNORD114", "SNORD114_number"), sep = "-"))
summary_targets114name <- summary_targets114name %>% mutate(SNORD114_number = as.numeric(gsub("[^0-9]", "", SNORD114_number))) %>% arrange(SNORD114_number)
summary_targets114name$snoRNA_number <- paste(summary_targets114name$SNORD114, summary_targets114name$SNORD114_number, sep = "-")
#Plot
  #SNORD115
plotsummary115_A <- ggplot(data = summary_targets115name, aes(x=snoRNA_number, y=count, fill=external_gene_name)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD115 Summary: Gene Targeting Events by snoRNA Copy", x ="SNORD copy", y = "Number of Targeting Events") + scale_x_discrete(limits = levels(chr15sno115$snoRNA_number))
plotsummary115_A
plotsummary115_B <- ggplot(data = summary_targets115name, aes(x=external_gene_name, y=count, fill=snoRNA_number)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD115 Summary: Gene Targeting Events by snoRNA Copy", x ="Target Gene Name", y = "Number of Targeting Events")
plotsummary115_B
  #SNORD116
plotsummary116_A <- ggplot(data = summary_targets116name, aes(x=snoRNA_number, y=count, fill=external_gene_name)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD116 Summary: Gene Targeting Events by snoRNA Copy", x ="SNORD copy", y = "Number of Targeting Events") + scale_x_discrete(limits = levels(chr15sno116$snoRNA_number))
plotsummary116_A
plotsummary116_B <- ggplot(data = summary_targets116name, aes(x=external_gene_name, y=count, fill=snoRNA_number)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD116 Summary: Gene Targeting Events by snoRNA Copy", x ="Target Gene Name", y = "Number of Targeting Events")
plotsummary116_B
  #SNORDchr15
plotsummarychr15_A <- ggplot(data = summary_targetschr15name, aes(x=snoRNA_number, y=count, fill=external_gene_name)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="Other SNORDs on Chr15 Summary: Gene Targeting Events by snoRNA Copy", x ="SNORD copy", y = "Number of Targeting Events") + scale_x_discrete(limits = levels(snocountchr15name$snoRNA_number))
plotsummarychr15_A
plotsummarychr15_B <- ggplot(data = summary_targetschr15name, aes(x=external_gene_name, y=count, fill=snoRNA_number)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="Other SNORDs on Chr15 Summary: Gene Targeting Events by snoRNA Copy", x ="Target Gene Name", y = "Number of Targeting Events")
plotsummarychr15_B
  #SNORD113
plotsummary113_A <- ggplot(data = summary_targets113name, aes(x=snoRNA_number, y=count, fill=external_gene_name)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD113 Summary: Gene Targeting Events by snoRNA Copy", x ="SNORD copy", y = "Number of Targeting Events") + scale_x_discrete(limits = levels(summary_targets113name$snoRNA_number))
plotsummary113_A
plotsummary113_B <- ggplot(data = summary_targets113name, aes(x=external_gene_name, y=count, fill=snoRNA_number)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD113 Summary: Gene Targeting Events by snoRNA Copy", x ="Target Gene Name", y = "Number of Targeting Events")
plotsummary113_B
  #SNORD114
plotsummary114_A <- ggplot(data = summary_targets114name, aes(x=snoRNA_number, y=count, fill=external_gene_name)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD114 Summary: Gene Targeting Events by snoRNA Copy", x ="SNORD copy", y = "Number of Targeting Events") + scale_x_discrete(limits = levels(snocount114name$snoRNA_number))
plotsummary114_A
plotsummary114_B <- ggplot(data = summary_targets114name, aes(x=external_gene_name, y=count, fill=snoRNA_number)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD114 Summary: Gene Targeting Events by snoRNA Copy", x ="Target Gene Name", y = "Number of Targeting Events")
plotsummary114_B
#Save as pdf
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/TargetsandCopieschr15.pdf", width = 10, height = 8)
plotsummary115_A
plotsummary115_B
plotsummary116_A
plotsummary116_B
plotsummarychr15_A
plotsummarychr15_B
dev.off()
#Save controls as pdf
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/TargetsandCopiesControls.pdf", width = 10, height = 8)
plotsummary113_A
plotsummary113_B
plotsummary114_A
plotsummary114_B
dev.off()

##Calculate difference between average/median targets per gene by SNORD116 groups I-III vs SNORD115
#Label SNORD116 copies by grouping
summary_targets116name$orig_group <- "116-I"
summary_targets116name <- summary_targets116name %>% mutate(orig_group = ifelse(SNORD116_number %in% 10:24, "116-II", orig_group))
summary_targets116name <- summary_targets116name %>% mutate(orig_group = ifelse(SNORD116_number %in% 25:30, "116-III", orig_group))
#Subset by group
groupI <- summary_targets116name[which(summary_targets116name$orig_group == "116-I"),]
groupII <- summary_targets116name[which(summary_targets116name$orig_group == "116-II"),]
groupIII <- summary_targets116name[which(summary_targets116name$orig_group == "116-III"),]
#Calculate sum of counts for each group
calc_115 <- summary_targets115name %>% group_by(external_gene_name) %>% summarize(sum115 = sum(count))
calc_116_I <- groupI %>% group_by(external_gene_name) %>% summarize(sum116_I = sum(count))
calc_116_II <- groupII %>% group_by(external_gene_name) %>% summarize(sum116_II = sum(count))
calc_116_III <- groupIII %>% group_by(external_gene_name) %>% summarize(sum116_III = sum(count))
#Merge group data with SNORD115 data
merge_I <- merge(calc_115, calc_116_I, by="external_gene_name", sort = TRUE, all = TRUE)
merge_II <- merge(calc_115, calc_116_II, by="external_gene_name", sort = TRUE, all = TRUE)
merge_III <- merge(calc_115, calc_116_III, by="external_gene_name", sort = TRUE, all = TRUE)
#Set NA = 0 for counts
merge_I[is.na(merge_I)] <- 0
merge_II[is.na(merge_II)] <- 0
merge_III[is.na(merge_III)] <- 0
#Add 1 to all sums
enrich_I <- merge_I %>% mutate(sum115_plus1 = sum115 + 1)
enrich_I <- enrich_I %>% mutate(sum116_I_plus1 = sum116_I + 1)
enrich_II <- merge_II %>% mutate(sum115_plus1 = sum115 + 1)
enrich_II <- enrich_II %>% mutate(sum116_II_plus1 = sum116_II + 1)
enrich_III <- merge_III %>% mutate(sum115_plus1 = sum115 + 1)
enrich_III <- enrich_III %>% mutate(sum116_III_plus1 = sum116_III + 1)
#Divide columns by number of sno copies per group
enrich_I$avgpercopy_115 <- enrich_I$sum115_plus1 / 48
enrich_I$avgpercopy_116_I <- enrich_I$sum116_I_plus1 / 9
enrich_II$avgpercopy_115 <- enrich_II$sum115_plus1 / 48
enrich_II$avgpercopy_116_II <- enrich_II$sum116_II_plus1 / 15
enrich_III$avgpercopy_115 <- enrich_III$sum115_plus1 / 48
enrich_III$avgpercopy_116_III <- enrich_III$sum116_III_plus1 / 6
#Divide avg for SNORD116 per group by avg for SNORD115
enrich_I$enrichment <- enrich_I$avgpercopy_116_I / enrich_I$avgpercopy_115
enrich_II$enrichment <- enrich_II$avgpercopy_116_II / enrich_II$avgpercopy_115
enrich_III$enrichment <- enrich_III$avgpercopy_116_III / enrich_III$avgpercopy_115
#Order by enrichment
enrich_I <- enrich_I %>% arrange(desc(enrichment))
enrich_II <- enrich_II %>% arrange(desc(enrichment))
enrich_III <- enrich_III %>% arrange(desc(enrichment))
#Save tables
write.csv(enrich_I, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-IvsSNORD115_genetargets.csv")
write.csv(enrich_II, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-IIvsSNORD115_genetargets.csv")
write.csv(enrich_III, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-IIIvsSNORD115_genetargets.csv")
#Filter for the average group targeting for a gene being 1 or greater
topenrich_I <- enrich_I[which(enrich_I$avgpercopy_116_I > 1 & enrich_I$enrichment > 1),]
topenrich_II <- enrich_II[which(enrich_II$avgpercopy_116_II > 1 & enrich_II$enrichment > 1),]
topenrich_III <- enrich_III[which(enrich_III$avgpercopy_116_III > 1 & enrich_III$enrichment > 1),]
#Merge
enriched_all <- merge(topenrich_I, topenrich_II, by="external_gene_name", sort = TRUE, all = TRUE)
enriched_all <- merge(enriched_all, topenrich_III, by="external_gene_name", sort = TRUE, all = TRUE)
enriched_shared <- merge(topenrich_I, topenrich_II, by="external_gene_name", sort = TRUE, all = FALSE)
enriched_shared <- merge(enriched_shared, topenrich_III, by="external_gene_name", sort = TRUE, all = FALSE)
enriched_list <- enriched_all$external_gene_name

##Parse data per gene based on list made above
#Create new data frame listing all SNORD copies
all_snord115 <- as.data.frame(unique(summary_targets115name$SNORD115_number))
names(all_snord115)[1] <- "SNORD115_number"
all_snord116 <- as.data.frame(unique(summary_targets116name$SNORD116_number))
names(all_snord116)[1] <- "SNORD116_number"
all_chr15other <- as.data.frame(unique(summary_targetschr15name$snoRNA_number))
names(all_chr15other)[1] <- "snoRNA_number"
all_snord113 <- as.data.frame(unique(summary_targets113name$snoRNA_number))
names(all_snord113)[1] <- "snoRNA_number"
all_snord114 <- as.data.frame(unique(summary_targets114name$snoRNA_number))
names(all_snord114)[1] <- "snoRNA_number"
#Create empty variable
snord115 <- list()
snord116 <- list()
chr15other <- list()
snord113 <- list()
snord114 <- list()
snord115table <- list()
snord116table <- list()
chr15othertable <- list()
snord113table <- list()
snord114table <- list()
allSNORDtable <- list()
plot <- list()
#Create plots with p-values
library("ggpubr")
for (i in enriched_list)
{
  #SNORD115
    #Subset
  snord115[[i]] <- subset(summary_targets115name, external_gene_name == i)
    #Merge with all SNORD115 copies
  snord115[[i]] <- merge(snord115[[i]], all_snord115, by="SNORD115_number", sort = TRUE, all = TRUE)
    #Set NA counts to 0
  snord115[[i]]$count[is.na(snord115[[i]]$count)] <- 0
    #Fill in SNORD115 column
  snord115[[i]] <- snord115[[i]] %>% mutate(SNORD115 = ifelse(SNORD115 %in% NA, "SNORD115", SNORD115))
    #Subset columns for plotting
  snord115table[[i]] <- snord115[[i]][c(4:5)]
  names(snord115table[[i]])[1] <- "SNORD_group"
  #SNORD116
    #Subset
  snord116[[i]]<- subset(summary_targets116name, external_gene_name == i)
    #Merge with all SNORD116 copies
  snord116[[i]] <- merge(snord116[[i]], all_snord116, by="SNORD116_number", sort = TRUE, all = TRUE)
    #Set NA counts to 0
  snord116[[i]]$count[is.na(snord116[[i]]$count)] <- 0
    #Fill in SNORD116 group column
  snord116[[i]]$orig_group <- "116-I"
  snord116[[i]] <- snord116[[i]] %>% mutate(orig_group = ifelse(SNORD116_number %in% 10:24, "116-II", orig_group))
  snord116[[i]] <- snord116[[i]] %>% mutate(orig_group = ifelse(SNORD116_number %in% 25:30, "116-III", orig_group))
    #Subset columns for plotting
  snord116table[[i]] <- snord116[[i]][c(5,7)]
  names(snord116table[[i]])[2] <- "SNORD_group"
  #chr15 other
    #Subset
  chr15other[[i]] <- subset(summary_targetschr15name, external_gene_name == i)
    #Merge with all SNORD115 copies
  chr15other[[i]] <- merge(chr15other[[i]], all_chr15other, by="snoRNA_number", sort = TRUE, all = TRUE)
    #Set NA counts to 0
  chr15other[[i]]$count[is.na(chr15other[[i]]$count)] <- 0
    #Fill in SNORD116 group column
  chr15other[[i]]$snoRNA_number <- "chr15other"
    #Subset columns for plotting
  chr15othertable[[i]] <- chr15other[[i]][c(1,4)]
  names(chr15othertable[[i]])[1] <- "SNORD_group"
  #SNORD113
    #Subset
  snord113[[i]] <- subset(summary_targets113name, external_gene_name == i)
    #Merge with all SNORD113 copies
  snord113[[i]] <- merge(snord113[[i]], all_snord113, by="snoRNA_number", sort = TRUE, all = TRUE)
    #Set NA counts to 0
  snord113[[i]]$count[is.na(snord113[[i]]$count)] <- 0
    #Fill in SNORD113 column
  snord113[[i]] <- snord113[[i]] %>% mutate(SNORD113 = ifelse(SNORD113 %in% NA, "SNORD113", SNORD113))
    #Subset columns for plotting
  snord113table[[i]] <- snord113[[i]][c(4,6)]
  names(snord113table[[i]])[1] <- "SNORD_group"
  #SNORD114
    #Subset
  snord114[[i]] <- subset(summary_targets114name, external_gene_name == i)
    #Merge with all SNORD114 copies
  snord114[[i]] <- merge(snord114[[i]], all_snord114, by="snoRNA_number", sort = TRUE, all = TRUE)
    #Set NA counts to 0
  snord114[[i]]$count[is.na(snord114[[i]]$count)] <- 0
    #Fill in SNORD114 column
  snord114[[i]] <- snord114[[i]] %>% mutate(SNORD114 = ifelse(SNORD114 %in% NA, "SNORD114", SNORD114))
    #Subset columns for plotting
  snord114table[[i]] <- snord114[[i]][c(4,6)]
  names(snord114table[[i]])[1] <- "SNORD_group"
  #Bind lists for plotting
  allSNORDtable[[i]] <- rbind(snord115table[[i]], snord116table[[i]], chr15othertable[[i]], snord113table[[i]], snord114table[[i]])
  #Set comparisons & arguments
  comparisons <- list(c("SNORD115", "116-I"), c("SNORD115", "116-II"), c("SNORD115", "116-III"), c("SNORD114", "116-I"), c("SNORD114", "116-II"), c("SNORD114", "116-III"), c("SNORD113", "116-I"), c("SNORD113", "116-II"), c("SNORD113", "116-III"), c("chr15other", "116-I"), c("chr15other", "116-II"), c("chr15other", "116-III"), c("116-I", "116-II"), c("116-I", "116-III"), c("116-II", "116-III"))
  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "NS"))
  #Plot
  plot[[i]] <- ggplot(allSNORDtable[[i]], aes(x = SNORD_group, y = count)) + geom_boxplot(outlier.shape = NA)+ geom_point(position = position_jitter(w = 0.2, h = 0.12), alpha = 0.35, size = 2, shape = 16) +
    labs(title= i, x ="SNORD Group", y = "Number of Targeting Events") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", method.args = list(exact = FALSE), symnum.args = symnum.args, label = "p.signif")
}
#Save plots
  #Start a new PDF device
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/BoxplotsCountsPerGenePvalsALLCOMPS.pdf", width = 10, height = 8)
for (i in enriched_list)
{
  #Print the plot to the PDF device
  print(plot[[i]])
}
  #Close the PDF device
dev.off()

##Do permutations on sets of 42 random genes to run through snoGloBe
library("EDASeq")
library("purrr")
#Read in shared genes list that snoGloBe was run on & counts table from RNAseq
shared <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/ENSEMBLsharedsmlgDEL.txt", sep="\t", head=F)
counts <- read.csv("../Cotney_Lab/PWS_RNASeq/STAR/lgDEL/CombinedH9CT2_Results_DESeqSTAR.csv", header = TRUE, row.names = "Gene")
#Subset counts table for WT H9
counts_H9WT <- counts[c(26:31)]
#Make list of all genes from counts object
counts_H9WT_list <- row.names(counts_H9WT)
#Get summary statistics for all genes
allStats <- as.data.frame(getGeneLengthAndGCContent(counts_H9WT_list, "hsa", mode="biomart"))
#Save stats
write.table(allStats, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/allGeneStats.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  #To read in file
allStats <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/allGeneStats.tsv", sep = "\t", row.names = 1, header = TRUE)
#Subset list of shared genes from overall list
shared_stats <- allStats[shared$V1,]
#Remove unexpressed genes from overall list
counts_H9WT_nonzero <- counts_H9WT[which(rowSums(counts_H9WT) > 0),]
allnonZeroStats <- allStats[row.names(counts_H9WT_nonzero),]
#Remove shared genes from overall list
statToSample <- allnonZeroStats[!(row.names(allnonZeroStats) %in% row.names(shared_stats)),]
#Figure out average expression of our gene list in WT-H9
shared_counts <- counts_H9WT_nonzero[shared$V1,]
shared_counts$mean <- rowMeans(shared_counts)
#Calculate median GC content, gene length, expression
gc <-shared_stats$gc
  #round(median(gc), 2) = 0.49
len <- shared_stats$length
  #round(median(len), 0) = 5564
count <- shared_counts$mean
  #round(median(count), 0) = 551
#Remove mitochondrial genes from list
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes.table <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name", "chromosome_name"), values= row.names(statToSample), mart= mart)
  #This gets rid of 402 genes --> this was done using https://dec2017.archive.ensembl.org/ before it was taken down
noMT <- genes.table[which(genes.table$chromosome_name != "MT"),]
  #This gets rid of 36 genes
statToSample_noMT <- statToSample[noMT$ensembl_gene_id,]
#Make gene list to sample from
statToSample_list <- row.names(statToSample_noMT)
  #Write out list for downstream analysis
write.table(statToSample_list, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/statToSample_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#Set seed for reproducibility of results
set.seed(42)
#Set the desired number of lists
desired_lists <- 100
#Set a counter for the number of lists added to the resulting_genes object
counter <- 0
#Create an empty list to store the resulting genes
resulting_genes <- list()
#Create an empty list to store the stats
stat <- list()
#Loop until the desired number of lists have been added to the resulting_genes object
while (counter < desired_lists)
{
  #Generate a random list of genes
  random_genes <- sample(statToSample_list, size = 42, replace = FALSE)
  #Subset summary statistics for shared gene list
  random_genes_stats <- statToSample_noMT[random_genes,]
  random_genes_counts <- counts_H9WT_nonzero[random_genes,]
  random_genes_counts$mean <- rowMeans(random_genes_counts)
  #Run Wilcoxon Test to test for significance between lists
  results_gc <- wilcox.test(random_genes_stats$gc,gc, exact = FALSE)
  results_len <- wilcox.test(random_genes_stats$length,len, exact = FALSE)
  results_count <- wilcox.test(random_genes_counts$mean,count, exact = FALSE)
  #Test the attributes of the genes
  if (results_gc$p.value > 0.05 & results_len$p.value > 0.05 & results_count$p.value > 0.05)
  {
    #Store the genes that meet the parameters
    resulting_genes[[counter + 1]] <- random_genes
    #Increase the counter by 1
    counter <- counter + 1
    #Store stats
    stat[[counter + 1]] <- c(results_gc, results_len, results_count)
  }
}
  #Took about 40 mins to generate 100 lists
#Save list of lists & stats
library("rlist")
list.save(resulting_genes, '../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/permutations_100/listsForPermutations.yaml')
list.save(stats, '../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/permutations_100/WilcoxonTestStatsforPermLists.yaml')
  #To load list of lists into R
library("yaml")
resulting_genes <- read_yaml("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/permutations_100/listsForPermutations.yaml")
#Assign names to each list
for (i in 1:length(resulting_genes)) 
{
  names(resulting_genes)[i] <- paste0("random_gene_list", i)
}
#Save lists as separate txt files
for (i in 1:length(resulting_genes))
{
  write.table(resulting_genes[[i]], file=paste0("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/permutations_100/", names(resulting_genes)[[i]], ".txt"),  row.names = FALSE, col.names = FALSE, quote = FALSE)
}
  #Have to change EOL conversion to Unix in Notepad++ or on cluster

##After running snoGloBe in cluster, analyze results
#Import files
for (i in 1:length(resulting_genes))
{
  #Generate the file name
  file_name <- paste0("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/permutations_100/random_gene_list", i, "_snord116chr15")
  #Read the file and assign it to an object
  assign(paste0("perm", i), as.data.frame(read.table(file_name, sep="\t", head=T)))
}
#Calculate statistics for SNORD116 groups on shared list of 42 genes
  #Count by SNORD copy
snocount116 <- as.data.frame(targets116 %>% count(sno_id))
  #Split name of SNORD copy
snocount116name <- as.data.frame(snocount116 %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
  #Remove strand from SNORD copy name
snocount116name$snoRNA_number <- sapply(snocount116name$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
  #Separate copy number from "SNORD116"
snocount116name <- as.data.frame(snocount116name %>% separate(snoRNA_number, c("SNORD116", "SNORD116_number"), sep = "-"))
  #Remove copies not on chr15
snocount116name <- na.omit(snocount116name)
  #Create column for SNORD group & parse by sno copy number
snocount116name$orig_group <- "116-I"
snocount116name <- snocount116name %>% mutate(orig_group = ifelse(SNORD116_number %in% 10:24, "116-II", orig_group))
snocount116name <- snocount116name %>% mutate(orig_group = ifelse(SNORD116_number %in% 25:30, "116-III", orig_group))
  #Calculate statistics per group
statscount_116 <- snocount116name %>% group_by(orig_group) %>% summarize(sum116 = sum(n), median116 = median(n), mean116 = mean(n))
#Calculate statistics for SNORD116 groups on random gene lists
  #Create an empty list to store the objects
perm <- list()
  #Loop over the number of files
for (i in 1:length(resulting_genes))
{
  #Get the object name
  perm_name <- paste0("perm", i)
  #Add the object to the list
  perm[[i]] <- get(perm_name)
}
  #Create empty objects
permcount <- list()
permcountname <- list()
statscount116 <- list()
  #Run loop
for (i in 1:length(perm))
{
  permcount[[i]] <- as.data.frame(perm[[i]] %>% count(sno_id))
  permcountname[[i]] <- as.data.frame(permcount[[i]] %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
  permcountname[[i]]$snoRNA_number <- sapply(permcountname[[i]]$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
  permcountname[[i]] <- as.data.frame(permcountname[[i]] %>% separate(snoRNA_number, c("SNORD116", "SNORD116_number"), sep = "-"))
  permcountname[[i]]$orig_group <- "116-I"
  permcountname[[i]] <- permcountname[[i]] %>% mutate(orig_group = ifelse(SNORD116_number %in% 10:24, "116-II", orig_group))
  permcountname[[i]] <- permcountname[[i]] %>% mutate(orig_group = ifelse(SNORD116_number %in% 25:30, "116-III", orig_group))
  statscount116[[i]] <- permcountname[[i]] %>% group_by(orig_group) %>% summarize(sum116 = sum(n), median116 = median(n), mean116 = mean(n))
}
#Parse for mean, median, sum
  #Initialize vectors to store all of the test-stats
Perm.test.stat_meanI <- rep(0, 100)
Perm.test.stat_meanII <- rep(0, 100)
Perm.test.stat_meanIII <- rep(0, 100)
Perm.test.stat_medianI <- rep(0, 100)
Perm.test.stat_medianII <- rep(0, 100)
Perm.test.stat_medianIII <- rep(0, 100)
Perm.test.stat_sumI <- rep(0, 100)
Perm.test.stat_sumII <- rep(0, 100)
Perm.test.stat_sumIII <- rep(0, 100)
  #Run loop
for (i in 1:length(perm))
{
  Perm.test.stat_meanI[[i]] <- statscount116[[i]]$mean116[1]
  Perm.test.stat_meanII[[i]] <- statscount116[[i]]$mean116[2]
  Perm.test.stat_meanIII[[i]] <- statscount116[[i]]$mean116[3]
  Perm.test.stat_medianI[[i]] <- statscount116[[i]]$median116[1]
  Perm.test.stat_medianII[[i]] <- statscount116[[i]]$median116[2]
  Perm.test.stat_medianIII[[i]] <- statscount116[[i]]$median116[3]
  Perm.test.stat_sumI[[i]] <- statscount116[[i]]$sum116[1]
  Perm.test.stat_sumII[[i]] <- statscount116[[i]]$sum116[2]
  Perm.test.stat_sumIII[[i]] <- statscount116[[i]]$sum116[3]
}
  #Make results into dataframe to use with ggplot2
Perm.test.stat_meanIdf <- as.data.frame(Perm.test.stat_meanI)
Perm.test.stat_meanIIdf <- as.data.frame(Perm.test.stat_meanII)
Perm.test.stat_meanIIIdf <- as.data.frame(Perm.test.stat_meanIII)
Perm.test.stat_medianIdf <- as.data.frame(Perm.test.stat_medianI)
Perm.test.stat_medianIIdf <- as.data.frame(Perm.test.stat_medianII)
Perm.test.stat_medianIIIdf <- as.data.frame(Perm.test.stat_medianIII)
Perm.test.stat_sumIdf <- as.data.frame(Perm.test.stat_sumI)
Perm.test.stat_sumIIdf <- as.data.frame(Perm.test.stat_sumII)
Perm.test.stat_sumIIIdf <- as.data.frame(Perm.test.stat_sumIII)
  #Plot individual histograms 
hist_group1_mean <- ggplot(Perm.test.stat_meanIdf, aes(x=Perm.test.stat_meanI)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = median(Perm.test.stat_meanI)), col='red', size=.6) +
  geom_vline(aes(xintercept = statscount_116$mean116[1]), col='purple', size=1, linetype="dashed") +
  labs(x="Mean Targeting Events 116-I", y = "Frequency")
hist_group1_mean
hist_group2_mean <- ggplot(Perm.test.stat_meanIIdf, aes(x=Perm.test.stat_meanII)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = median(Perm.test.stat_meanII)), col='red', size=.6) +
  geom_vline(aes(xintercept = statscount_116$mean116[2]), col='purple', size=1, linetype="dashed") +
  labs(x="Mean Targeting Events 116-II", y = "Frequency")
hist_group2_mean
hist_group3_mean <- ggplot(Perm.test.stat_meanIIIdf, aes(x=Perm.test.stat_meanIII)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = median(Perm.test.stat_meanIII)), col='red', size=.6) +
  geom_vline(aes(xintercept = statscount_116$mean116[3]), col='purple', size=1, linetype="dashed") +
  labs(x="Mean Targeting Events 116-III", y = "Frequency")
hist_group3_mean
hist_group1_med <- ggplot(Perm.test.stat_medianIdf, aes(x=Perm.test.stat_medianI)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = median(Perm.test.stat_medianI)), col='red', size=.6) +
  geom_vline(aes(xintercept = statscount_116$median116[1]), col='purple', size=1, linetype="dashed") +
  labs(x="Median Targeting Events 116-I", y = "Frequency")
hist_group1_med
hist_group2_med <- ggplot(Perm.test.stat_medianIIdf, aes(x=Perm.test.stat_medianII)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = median(Perm.test.stat_medianII)), col='red', size=.6) +
  geom_vline(aes(xintercept = statscount_116$median116[2]), col='purple', size=1, linetype="dashed") +
  labs(x="Median Targeting Events 116-II", y = "Frequency")
hist_group2_med
hist_group3_med <- ggplot(Perm.test.stat_medianIIIdf, aes(x=Perm.test.stat_medianIII)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = median(Perm.test.stat_medianIII)), col='red', size=.6) +
  geom_vline(aes(xintercept = statscount_116$median116[3]), col='purple', size=1, linetype="dashed") +
  labs(x="Median Targeting Events 116-III", y = "Frequency")
hist_group3_med
hist_group1_sum <- ggplot(Perm.test.stat_sumIdf, aes(x=Perm.test.stat_sumI)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = median(Perm.test.stat_sumI)), col='red', size=.6) +
  geom_vline(aes(xintercept = statscount_116$sum116[1]), col='purple', size=1, linetype="dashed") +
  labs(x="Sum of Targeting Events 116-I", y = "Frequency")
hist_group1_sum
hist_group2_sum <- ggplot(Perm.test.stat_sumIIdf, aes(x=Perm.test.stat_sumII)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = median(Perm.test.stat_sumII)), col='red', size=.6) +
  geom_vline(aes(xintercept = statscount_116$sum116[2]), col='purple', size=1, linetype="dashed") +
  labs(x="Sum of Targeting Events 116-II", y = "Frequency")
hist_group2_sum
hist_group3_sum <- ggplot(Perm.test.stat_sumIIIdf, aes(x=Perm.test.stat_sumIII)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = median(Perm.test.stat_sumIII)), col='red', size=.6) +
  geom_vline(aes(xintercept = statscount_116$sum116[3]), col='purple', size=1, linetype="dashed") +
  labs(x="Sum of Targeting Events 116-III", y = "Frequency")
hist_group3_sum
  #To arrange plots
allhist <- ggarrange(hist_group1_mean, hist_group2_mean, hist_group3_mean, hist_group1_med, hist_group2_med, hist_group3_med, hist_group1_sum, hist_group2_sum, hist_group3_sum,
                     ncol = 3, nrow = 3)
allhist
  #Save plot
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/permutations_100/histogram_permutations.pdf", width = 10, height = 8)
allhist
dev.off()
#Parse for raw counts - ONLY did this for first set of permutations because the result can't be compared across groups (different # of copies/group)
  #Create empty variable
perm_rawcounts <- list()
  #Run loop to subset count & group columns
for (i in 1:10)
{
  perm_rawcounts[[i]] <- permcountname[[i]][,c(4:5)]
}
  #Bind into one dataframe
rawcounts_df <- bind_rows(perm_rawcounts)
  #Sort by group
rawcount_group1 <- rawcounts_df[which(rawcounts_df$orig_group == "116-I"),]
rawcount_group2 <- rawcounts_df[which(rawcounts_df$orig_group == "116-II"),]
rawcount_group3 <- rawcounts_df[which(rawcounts_df$orig_group == "116-III"),]
  #Do same with shared list, but make into list
sharedcount_group1 <- snocount116name[which(snocount116name$orig_group == "116-I"),]
sharedcount_group2 <- snocount116name[which(snocount116name$orig_group == "116-II"),]
sharedcount_group3 <- snocount116name[which(snocount116name$orig_group == "116-III"),]
  #Plot individual histograms
hist_group1_counts <- ggplot(rawcount_group1, aes(x=n)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = median(rawcount_group1$n)), col='red', size=.6) +
  geom_vline(aes(xintercept = median(sharedcount_group1$n)), col='purple', size=1, linetype="dashed") +
  labs(x="Raw Number Targeting Events for 116-I", y = "Frequency")
hist_group1_counts
hist_group2_counts <- ggplot(rawcount_group2, aes(x=n)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = median(rawcount_group2$n)), col='red', size=.6) +
  geom_vline(aes(xintercept = median(sharedcount_group2$n)), col='purple', size=1, linetype="dashed") +
  labs(x="Raw Number Targeting Events for 116-II", y = "Frequency")
hist_group2_counts
hist_group3_counts <- ggplot(rawcount_group3, aes(x=n)) + 
  geom_histogram(binwidth=1, color="black", fill="white") +
  geom_vline(aes(xintercept = median(rawcount_group3$n)), col='red', size=.6) +
  geom_vline(aes(xintercept = median(sharedcount_group3$n)), col='purple', size=1, linetype="dashed") +
  labs(x="Raw Number of Targeting Events 116-III", y = "Frequency")
hist_group3_counts
  #To arrange plots
rawhist <- ggarrange(hist_group1_counts, hist_group2_counts, hist_group3_counts,
                     ncol = 3, nrow = 1)
rawhist
  #Save plot
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/permutations_10/histogram_permutationsRawCounts.pdf", width = 10, height = 8)
rawhist
dev.off()

##Calculating % of predicted interactions across transcriptome
#Import files
bkgrnd_coverage <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/background_grouped_table.txt", sep="\t", head=F))
bkgrndShared_coverage <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/backgroundShared_grouped_table.txt", sep="\t", head=F))
bkgrndRandom_coverage <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/backgroundRandom_grouped_table.txt", sep="\t", head=F))
SNORD116_coverage <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/allSNORD116_grouped_table.txt", sep="\t", head=F))
SNORD116_coverageRandom <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/allSNORD116_randomPerms_groupedTableSample.txt", sep="\t", head=F))
SNORD115_coverage <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/snord115_grouped_table.txt", sep="\t", head=F))
otherchr15SNORDs_coverage <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/otherchr15SNORDs_grouped_table.txt", sep="\t", head=F))
#For background sets
  #Rename columns
colnames(bkgrnd_coverage) <- c("orig_category", "basepairs")
colnames(bkgrndShared_coverage) <- c("orig_category", "basepairs")
colnames(bkgrndRandom_coverage) <- c("orig_category", "basepairs")
  #Create column for category
bkgrnd_coverage$category <- c("Other", "3'UTR_only", "5'UTR_only", "Other", "CDS_only", "3'UTR+CDS", "5'UTR+CDS", "Other", "Junction", "Junction", "Junction", "Junction", "Junction", "Junction", "Junction", "Junction", "Intron", "Junction")
bkgrndShared_coverage$category <- c("Other", "3'UTR_only", "5'UTR_only", "CDS_only", "3'UTR+CDS", "5'UTR+CDS", "Junction", "Junction", "Junction", "Junction", "Junction", "Intron", "Junction")
bkgrndRandom_coverage$category <- c("Other", "3'UTR_only", "5'UTR_only", "Other", "CDS_only", "3'UTR+CDS", "5'UTR+CDS", "Other", "Junction", "Junction", "Junction", "Junction", "Junction", "Junction", "Junction", "Intron", "Junction")
    #To save
  write.table(bkgrnd_coverage, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/background_genomicFeatures_tablecategories.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(bkgrndShared_coverage, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/backgroundShared_genomicFeatures_tablecategories.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(bkgrndRandom_coverage, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/backgroundRandom_genomicFeatures_tablecategories.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  #Group
grouped_bkgrndcov <- bkgrnd_coverage %>% group_by(category) %>% summarize(basepairs = sum(basepairs))
grouped_bkgrndcovShared <- bkgrndShared_coverage %>% group_by(category) %>% summarize(basepairs = sum(basepairs))
grouped_bkgrndcovRandom <- bkgrndRandom_coverage %>% group_by(category) %>% summarize(basepairs = sum(basepairs))
  #Calculate percentage of exon, intron, and intron-exon junction of genes expressed in iNeurons transcriptome
summary_bkgrndcov <- grouped_bkgrndcov
summary_bkgrndcov$percent <- (summary_bkgrndcov$basepairs/(sum(summary_bkgrndcov$basepairs))*100)
summary_bkgrndcovShared <- grouped_bkgrndcovShared
summary_bkgrndcovShared$percent <- (summary_bkgrndcovShared$basepairs/(sum(summary_bkgrndcovShared$basepairs))*100)
summary_bkgrndcovRandom <- grouped_bkgrndcovRandom
summary_bkgrndcovRandom$percent <- (summary_bkgrndcovRandom$basepairs/(sum(summary_bkgrndcovRandom$basepairs))*100)
    #To save
  write.table(summary_bkgrndcov, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/background_genomicFeatures_calctable.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(summary_bkgrndcovShared, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/backgroundShared_genomicFeatures_calctable.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(summary_bkgrndcovRandom, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/backgroundRandom_genomicFeatures_calctable.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  #Modify to plot barplot
summary_bkgrndcov_bar <- summary_bkgrndcov
summary_bkgrndcov_bar$bar_category <- "Exon"
summary_bkgrndcov_bar <- summary_bkgrndcov_bar %>% mutate(bar_category = ifelse(category %in% "Intron", "Intron", bar_category))
summary_bkgrndcov_bar <- summary_bkgrndcov_bar %>% mutate(bar_category = ifelse(category %in% "Junction", "Junction", bar_category))
summary_bkgrndcov_bar$dataset <- "Background: iNeuron Transcriptome"
summary_bkgrndcovShared_bar <- summary_bkgrndcovShared
summary_bkgrndcovShared_bar$bar_category <- "Exon"
summary_bkgrndcovShared_bar <- summary_bkgrndcovShared_bar %>% mutate(bar_category = ifelse(category %in% "Intron", "Intron", bar_category))
summary_bkgrndcovShared_bar <- summary_bkgrndcovShared_bar %>% mutate(bar_category = ifelse(category %in% "Junction", "Junction", bar_category))
summary_bkgrndcovShared_bar$dataset <- "Background: Shared Genes"
summary_bkgrndcovRandom_bar <- summary_bkgrndcovRandom
summary_bkgrndcovRandom_bar$bar_category <- "Exon"
summary_bkgrndcovRandom_bar <- summary_bkgrndcovRandom_bar %>% mutate(bar_category = ifelse(category %in% "Intron", "Intron", bar_category))
summary_bkgrndcovRandom_bar <- summary_bkgrndcovRandom_bar %>% mutate(bar_category = ifelse(category %in% "Junction", "Junction", bar_category))
summary_bkgrndcovRandom_bar$dataset <- "Background: Random Permutations"
#For SNORD116 coverage on shared genes
  #Rename columns
colnames(SNORD116_coverage) <- c("orig_category", "basepairs")
  #Create column for SNORD group
SNORD116_coverage_bygroup <- SNORD116_coverage
SNORD116_coverage_bygroup$group <- SNORD116_coverage_bygroup$orig_category
SNORD116_coverage_bygroup <- as.data.frame(SNORD116_coverage_bygroup %>% separate(group, c("SNORDcopy"), sep = ","))
SNORD116_coverage_bygroup <- as.data.frame(SNORD116_coverage_bygroup %>% separate(SNORDcopy, c("SNORD116", "SNORD116_number"), sep = "-"))
  #Parse by sno copy number
SNORD116_coverage_bygroup$group <- "116-I"
SNORD116_coverage_bygroup <- SNORD116_coverage_bygroup %>% mutate(group = ifelse(SNORD116_number %in% 10:24, "116-II", group))
SNORD116_coverage_bygroup <- SNORD116_coverage_bygroup %>% mutate(group = ifelse(SNORD116_number %in% 25:30, "116-III", group))
  #Remove unnecessary columns
SNORD116_coverage_bygroup <- SNORD116_coverage_bygroup[-c(3)]
  #Get rid of copy number in first column
SNORD116_coverage_bygroup$orig_category <- sapply(SNORD116_coverage_bygroup$orig_category, function(x) sub("SNORD116-\\d+,", "SNORD116,", x))
  #Group by original category
SNORD116_coverage_bygroupALL <- SNORD116_coverage_bygroup %>% group_by(orig_category) %>% summarize(basepairs = sum(basepairs))
  #Create column for category
SNORD116_coverage_bygroupALL$category <- c("Other", "3'UTR_only", "5'UTR_only", "CDS_only", "5'UTR+CDS", "Junction", "Intron", "Junction")
    #To save
  write.table(SNORD116_coverage_bygroupALL, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/allSNORD116_genomicFeatures_tablecategories.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  #Group by new category
SNORD116cov_bycatALL <- SNORD116_coverage_bygroupALL %>% group_by(category) %>% summarize(basepairs = sum(basepairs))
  #Calculate percentage of exon, intron, and intron-exon junction of all SNORD116 binding events
summary_SNORD116cov_bycatALL <- SNORD116cov_bycatALL
summary_SNORD116cov_bycatALL$percent <- (summary_SNORD116cov_bycatALL$basepairs/(sum(summary_SNORD116cov_bycatALL$basepairs))*100)
    #To save
  write.table(summary_SNORD116cov_bycatALL, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/allSNORD116_genomicFeatures_calctable.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  #Modify to plot barplot
summary_SNORD116sharedALL_bar <- summary_SNORD116cov_bycatALL
summary_SNORD116sharedALL_bar$bar_category <- "Exon"
summary_SNORD116sharedALL_bar <- summary_SNORD116sharedALL_bar %>% mutate(bar_category = ifelse(category %in% "Intron", "Intron", bar_category))
summary_SNORD116sharedALL_bar <- summary_SNORD116sharedALL_bar %>% mutate(bar_category = ifelse(category %in% "Junction", "Junction", bar_category))
summary_SNORD116sharedALL_bar$dataset <- "SNORD116all_vs_sharedGenes"
  #Subset by group
SNORD116cov_groupI <- SNORD116_coverage_bygroup[which(SNORD116_coverage_bygroup$group == "116-I"),]
SNORD116cov_groupII <- SNORD116_coverage_bygroup[which(SNORD116_coverage_bygroup$group == "116-II"),]
SNORD116cov_groupIII <- SNORD116_coverage_bygroup[which(SNORD116_coverage_bygroup$group == "116-III"),]
SNORD116cov_groupI_III <- SNORD116_coverage_bygroup %>% filter(!SNORD116_coverage_bygroup$group %in% "116-II")
  #Group by original category
SNORD116cov_groupI_cat <- SNORD116cov_groupI %>% group_by(orig_category) %>% summarize(basepairs = sum(basepairs))
SNORD116cov_groupII_cat <- SNORD116cov_groupII %>% group_by(orig_category) %>% summarize(basepairs = sum(basepairs))
SNORD116cov_groupIII_cat <- SNORD116cov_groupIII %>% group_by(orig_category) %>% summarize(basepairs = sum(basepairs))
SNORD116cov_groupI_III_cat <- SNORD116cov_groupI_III %>% group_by(orig_category) %>% summarize(basepairs = sum(basepairs))
  #Create column for category
SNORD116cov_groupI_cat$category <- c("3'UTR_only", "5'UTR_only", "CDS_only", "Intron")
SNORD116cov_groupII_cat$category <- c("3'UTR_only", "5'UTR_only", "CDS_only", "Junction", "Intron", "Junction")
SNORD116cov_groupIII_cat$category <- c("Other", "3'UTR_only", "5'UTR_only", "CDS_only", "5'UTR+CDS", "Junction", "Intron", "Junction")
SNORD116cov_groupI_III_cat$category <- c("Other", "3'UTR_only", "5'UTR_only", "CDS_only", "5'UTR+CDS", "Junction", "Intron", "Junction")
    #To save
  write.table(SNORD116cov_groupI_cat, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-I_genomicFeatures_tablecategories.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(SNORD116cov_groupII_cat, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-II_genomicFeatures_tablecategories.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(SNORD116cov_groupIII_cat, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-III_genomicFeatures_tablecategories.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(SNORD116cov_groupI_III_cat, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-IandIII_genomicFeatures_tablecategories.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  #Group by new category
SNORD116cov_groupI_calc <- SNORD116cov_groupI_cat %>% group_by(category) %>% summarize(basepairs = sum(basepairs))
SNORD116cov_groupII_calc <- SNORD116cov_groupII_cat %>% group_by(category) %>% summarize(basepairs = sum(basepairs))
SNORD116cov_groupIII_calc <- SNORD116cov_groupIII_cat %>% group_by(category) %>% summarize(basepairs = sum(basepairs))
SNORD116cov_groupI_III_calc <- SNORD116cov_groupI_III_cat %>% group_by(category) %>% summarize(basepairs = sum(basepairs))
  #Calculate percentage of exon, intron, and intron-exon junction of all SNORD116 binding events
summary_SNORD116cov_groupI_calc <- SNORD116cov_groupI_calc
summary_SNORD116cov_groupI_calc$percent <- (summary_SNORD116cov_groupI_calc$basepairs/(sum(summary_SNORD116cov_groupI_calc$basepairs))*100)
summary_SNORD116cov_groupII_calc <- SNORD116cov_groupII_calc
summary_SNORD116cov_groupII_calc$percent <- (summary_SNORD116cov_groupII_calc$basepairs/(sum(summary_SNORD116cov_groupII_calc$basepairs))*100)
summary_SNORD116cov_groupIII_calc <- SNORD116cov_groupIII_calc
summary_SNORD116cov_groupIII_calc$percent <- (summary_SNORD116cov_groupIII_calc$basepairs/(sum(summary_SNORD116cov_groupIII_calc$basepairs))*100)
summary_SNORD116cov_groupI_III_calc <- SNORD116cov_groupI_III_calc
summary_SNORD116cov_groupI_III_calc$percent <- (summary_SNORD116cov_groupI_III_calc$basepairs/(sum(summary_SNORD116cov_groupI_III_calc$basepairs))*100)
    #To save
  write.table(summary_SNORD116cov_groupI_calc, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-I_genomicFeatures_calctable.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(summary_SNORD116cov_groupII_calc, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-II_genomicFeatures_calctable.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(summary_SNORD116cov_groupIII_calc, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-III_genomicFeatures_calctable.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(summary_SNORD116cov_groupI_III_calc, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-IandIII_genomicFeatures_calctable.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  #Modify to plot barplot
summary_SNORD116groupI_shared_bar <- summary_SNORD116cov_groupI_calc
summary_SNORD116groupI_shared_bar$bar_category <- "Exon"
summary_SNORD116groupI_shared_bar <- summary_SNORD116groupI_shared_bar %>% mutate(bar_category = ifelse(category %in% "Intron", "Intron", bar_category))
summary_SNORD116groupI_shared_bar <- summary_SNORD116groupI_shared_bar %>% mutate(bar_category = ifelse(category %in% "Junction", "Junction", bar_category))
summary_SNORD116groupI_shared_bar$dataset <- "SNORD116-I_vs_sharedGenes"
summary_SNORD116groupII_shared_bar <- summary_SNORD116cov_groupII_calc
summary_SNORD116groupII_shared_bar$bar_category <- "Exon"
summary_SNORD116groupII_shared_bar <- summary_SNORD116groupII_shared_bar %>% mutate(bar_category = ifelse(category %in% "Intron", "Intron", bar_category))
summary_SNORD116groupII_shared_bar <- summary_SNORD116groupII_shared_bar %>% mutate(bar_category = ifelse(category %in% "Junction", "Junction", bar_category))
summary_SNORD116groupII_shared_bar$dataset <- "SNORD116-II_vs_sharedGenes"
summary_SNORD116groupIII_shared_bar <- summary_SNORD116cov_groupIII_calc
summary_SNORD116groupIII_shared_bar$bar_category <- "Exon"
summary_SNORD116groupIII_shared_bar <- summary_SNORD116groupIII_shared_bar %>% mutate(bar_category = ifelse(category %in% "Intron", "Intron", bar_category))
summary_SNORD116groupIII_shared_bar <- summary_SNORD116groupIII_shared_bar %>% mutate(bar_category = ifelse(category %in% "Junction", "Junction", bar_category))
summary_SNORD116groupIII_shared_bar$dataset <- "SNORD116-III_vs_sharedGenes"
summary_SNORD116groupI_III_shared_bar <- summary_SNORD116cov_groupI_III_calc
summary_SNORD116groupI_III_shared_bar$bar_category <- "Exon"
summary_SNORD116groupI_III_shared_bar <- summary_SNORD116groupI_III_shared_bar %>% mutate(bar_category = ifelse(category %in% "Intron", "Intron", bar_category))
summary_SNORD116groupI_III_shared_bar <- summary_SNORD116groupI_III_shared_bar %>% mutate(bar_category = ifelse(category %in% "Junction", "Junction", bar_category))
summary_SNORD116groupI_III_shared_bar$dataset <- "SNORD116-IandIII_vs_sharedGenes"
#For SNORD116 coverage on random gene lists
  #Rename columns
colnames(SNORD116_coverageRandom) <- c("orig_category", "basepairs", "permutation")
  #Create column for SNORD group
SNORD116_covRandom_bygroup <- SNORD116_coverageRandom
SNORD116_covRandom_bygroup$group <- SNORD116_covRandom_bygroup$orig_category
SNORD116_covRandom_bygroup <- as.data.frame(SNORD116_covRandom_bygroup %>% separate(group, c("SNORDcopy"), sep = ","))
SNORD116_covRandom_bygroup <- as.data.frame(SNORD116_covRandom_bygroup %>% separate(SNORDcopy, c("SNORD116", "SNORD116_number"), sep = "-"))
  #Parse by sno copy number
SNORD116_covRandom_bygroup$group <- "116-I"
SNORD116_covRandom_bygroup <- SNORD116_covRandom_bygroup %>% mutate(group = ifelse(SNORD116_number %in% 10:24, "116-II", group))
SNORD116_covRandom_bygroup <- SNORD116_covRandom_bygroup %>% mutate(group = ifelse(SNORD116_number %in% 25:30, "116-III", group))
  #Remove unnecessary column
SNORD116_covRandom_bygroup <- SNORD116_covRandom_bygroup[-c(4)]
  #Get rid of copy number in first column
SNORD116_covRandom_bygroup$orig_category <- sapply(SNORD116_covRandom_bygroup$orig_category, function(x) sub("SNORD116-\\d+,", "SNORD116,", x))
  #Group by original category
SNORD116_covRandom_bygroupALL <- SNORD116_covRandom_bygroup %>% group_by(orig_category) %>% summarize(basepairs = sum(basepairs))
  #Create column for category
SNORD116_covRandom_bygroupALL$category <- c("Other", "3'UTR_only", "5'UTR_only", "Other", "CDS_only", "3'UTR+CDS", "5'UTR+CDS", "Other", "Junction", "Junction", "Junction", "Junction", "Junction", "Junction", "Intron", "Junction")
    #To save
  write.table(SNORD116_covRandom_bygroupALL, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/allSNORD116vsRandom_genomicFeatures_tablecategories.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  #Group by new category
SNORD116covRandom_bycatALL <- SNORD116_covRandom_bygroupALL %>% group_by(category) %>% summarize(basepairs = sum(basepairs))
  #Calculate percentage of exon, intron, and intron-exon junction of all SNORD116 binding events
summary_SNORD116covRandom_bycatALL <- SNORD116covRandom_bycatALL
summary_SNORD116covRandom_bycatALL$percent <- (summary_SNORD116covRandom_bycatALL$basepairs/(sum(summary_SNORD116covRandom_bycatALL$basepairs))*100)
    #To save
  write.table(summary_SNORD116covRandom_bycatALL, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/allSNORD116vRandom_genomicFeatures_calctable.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  #Modify to plot barplot
summary_SNORD116randomALL_bar <- summary_SNORD116covRandom_bycatALL
summary_SNORD116randomALL_bar$bar_category <- "Exon"
summary_SNORD116randomALL_bar <- summary_SNORD116randomALL_bar %>% mutate(bar_category = ifelse(category %in% "Intron", "Intron", bar_category))
summary_SNORD116randomALL_bar <- summary_SNORD116randomALL_bar %>% mutate(bar_category = ifelse(category %in% "Junction", "Junction", bar_category))
summary_SNORD116randomALL_bar$dataset <- "SNORD116all_vs_randomPerm"
  #Subset by group
SNORD116covRandom_groupI <- SNORD116_covRandom_bygroup[which(SNORD116_covRandom_bygroup$group == "116-I"),]
SNORD116covRandom_groupII <- SNORD116_covRandom_bygroup[which(SNORD116_covRandom_bygroup$group == "116-II"),]
SNORD116covRandom_groupIII <- SNORD116_covRandom_bygroup[which(SNORD116_covRandom_bygroup$group == "116-III"),]
SNORD116covRandom_groupI_III <- SNORD116_covRandom_bygroup %>% filter(!SNORD116_covRandom_bygroup$group %in% "116-II")
  #Group by original category
SNORD116covRandom_groupI_cat <- SNORD116covRandom_groupI %>% group_by(orig_category) %>% summarize(basepairs = sum(basepairs))
SNORD116covRandom_groupII_cat <- SNORD116covRandom_groupII %>% group_by(orig_category) %>% summarize(basepairs = sum(basepairs))
SNORD116covRandom_groupIII_cat <- SNORD116covRandom_groupIII %>% group_by(orig_category) %>% summarize(basepairs = sum(basepairs))
SNORD116covRandom_groupI_III_cat <- SNORD116covRandom_groupI_III %>% group_by(orig_category) %>% summarize(basepairs = sum(basepairs))
  #Create column for category
SNORD116covRandom_groupI_cat$category <- c("Other", "3'UTR_only", "5'UTR_only", "Other", "CDS_only", "3'UTR+CDS", "5'UTR+CDS", "Junction", "Junction", "Junction", "Junction", "Intron", "Junction")
SNORD116covRandom_groupII_cat$category <- c("Other", "3'UTR_only", "5'UTR_only", "CDS_only", "3'UTR+CDS", "5'UTR+CDS", "Other", "Junction", "Junction", "Junction", "Junction", "Junction", "Intron", "Junction")
SNORD116covRandom_groupIII_cat$category <- c("Other", "3'UTR_only", "5'UTR_only", "CDS_only", "3'UTR+CDS", "5'UTR+CDS", "Junction", "Junction", "Junction", "Junction", "Junction", "Intron", "Junction")
SNORD116covRandom_groupI_III_cat$category <- c("Other", "3'UTR_only", "5'UTR_only", "Other", "CDS_only", "3'UTR+CDS", "5'UTR+CDS", "Junction", "Junction", "Junction", "Junction", "Junction", "Intron", "Junction")
    #To save
  write.table(SNORD116covRandom_groupI_cat, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-IvRandom_genomicFeatures_tablecategories.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(SNORD116covRandom_groupII_cat, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-IIvRandom_genomicFeatures_tablecategories.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(SNORD116covRandom_groupIII_cat, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-IIIvRandom_genomicFeatures_tablecategories.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(SNORD116covRandom_groupI_III_cat, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-IandIIIvRandom_genomicFeatures_tablecategories.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  #Group by new category
SNORD116covRandom_groupI_calc <- SNORD116covRandom_groupI_cat %>% group_by(category) %>% summarize(basepairs = sum(basepairs))
SNORD116covRandom_groupII_calc <- SNORD116covRandom_groupII_cat %>% group_by(category) %>% summarize(basepairs = sum(basepairs))
SNORD116covRandom_groupIII_calc <- SNORD116covRandom_groupIII_cat %>% group_by(category) %>% summarize(basepairs = sum(basepairs))
SNORD116covRandom_groupI_III_calc <- SNORD116covRandom_groupI_III_cat %>% group_by(category) %>% summarize(basepairs = sum(basepairs))
  #Calculate percentage of exon, intron, and intron-exon junction of all SNORD116 binding events
summary_SNORD116covRandom_groupI_calc <- SNORD116covRandom_groupI_calc
summary_SNORD116covRandom_groupI_calc$percent <- (summary_SNORD116covRandom_groupI_calc$basepairs/(sum(summary_SNORD116covRandom_groupI_calc$basepairs))*100)
summary_SNORD116covRandom_groupII_calc <- SNORD116covRandom_groupII_calc
summary_SNORD116covRandom_groupII_calc$percent <- (summary_SNORD116covRandom_groupII_calc$basepairs/(sum(summary_SNORD116covRandom_groupII_calc$basepairs))*100)
summary_SNORD116covRandom_groupIII_calc <- SNORD116covRandom_groupIII_calc
summary_SNORD116covRandom_groupIII_calc$percent <- (summary_SNORD116covRandom_groupIII_calc$basepairs/(sum(summary_SNORD116covRandom_groupIII_calc$basepairs))*100)
summary_SNORD116covRandom_groupI_III_calc <- SNORD116covRandom_groupI_III_calc
summary_SNORD116covRandom_groupI_III_calc$percent <- (summary_SNORD116covRandom_groupI_III_calc$basepairs/(sum(summary_SNORD116covRandom_groupI_III_calc$basepairs))*100)
    #To save
  write.table(summary_SNORD116covRandom_groupI_calc, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-IvRandom_genomicFeatures_calctable.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(summary_SNORD116covRandom_groupII_calc, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-IIvRandom_genomicFeatures_calctable.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(summary_SNORD116covRandom_groupIII_calc, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-IIIvRandom_genomicFeatures_calctable.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(summary_SNORD116covRandom_groupI_III_calc, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116-IandIIIvRandom_genomicFeatures_calctable.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  #Modify to plot barplot
summary_SNORD116groupI_random_bar <- summary_SNORD116covRandom_groupI_calc
summary_SNORD116groupI_random_bar$bar_category <- "Exon"
summary_SNORD116groupI_random_bar <- summary_SNORD116groupI_random_bar %>% mutate(bar_category = ifelse(category %in% "Intron", "Intron", bar_category))
summary_SNORD116groupI_random_bar <- summary_SNORD116groupI_random_bar %>% mutate(bar_category = ifelse(category %in% "Junction", "Junction", bar_category))
summary_SNORD116groupI_random_bar$dataset <- "SNORD116-I_vs_randomPerm"
summary_SNORD116groupII_random_bar <- summary_SNORD116covRandom_groupII_calc
summary_SNORD116groupII_random_bar$bar_category <- "Exon"
summary_SNORD116groupII_random_bar <- summary_SNORD116groupII_random_bar %>% mutate(bar_category = ifelse(category %in% "Intron", "Intron", bar_category))
summary_SNORD116groupII_random_bar <- summary_SNORD116groupII_random_bar %>% mutate(bar_category = ifelse(category %in% "Junction", "Junction", bar_category))
summary_SNORD116groupII_random_bar$dataset <- "SNORD116-II_vs_randomPerm"
summary_SNORD116groupIII_random_bar <- summary_SNORD116covRandom_groupIII_calc
summary_SNORD116groupIII_random_bar$bar_category <- "Exon"
summary_SNORD116groupIII_random_bar <- summary_SNORD116groupIII_random_bar %>% mutate(bar_category = ifelse(category %in% "Intron", "Intron", bar_category))
summary_SNORD116groupIII_random_bar <- summary_SNORD116groupIII_random_bar %>% mutate(bar_category = ifelse(category %in% "Junction", "Junction", bar_category))
summary_SNORD116groupIII_random_bar$dataset <- "SNORD116-III_vs_randomPerm"
summary_SNORD116groupI_III_random_bar <- summary_SNORD116covRandom_groupI_III_calc
summary_SNORD116groupI_III_random_bar$bar_category <- "Exon"
summary_SNORD116groupI_III_random_bar <- summary_SNORD116groupI_III_random_bar %>% mutate(bar_category = ifelse(category %in% "Intron", "Intron", bar_category))
summary_SNORD116groupI_III_random_bar <- summary_SNORD116groupI_III_random_bar %>% mutate(bar_category = ifelse(category %in% "Junction", "Junction", bar_category))
summary_SNORD116groupI_III_random_bar$dataset <- "SNORD116-IandIII_vs_randomPerm"
#For SNORD115 coverage on shared genes
  #Rename columns
colnames(SNORD115_coverage) <- c("orig_category", "basepairs")
  #Create column for category
SNORD115_coverage$category <- c("3'UTR_only", "5'UTR_only", "CDS_only", "5'UTR+CDS", "Intron")
    #To save
  write.table(SNORD115_coverage, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/allSNORD115_genomicFeatures_tablecategories.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  #Calculate percentage of exon, intron, and intron-exon junction of all SNORD115 binding events
summary_SNORD115cov <- SNORD115_coverage
summary_SNORD115cov$percent <- (summary_SNORD115cov$basepairs/(sum(summary_SNORD115cov$basepairs))*100)
    #To save
  write.table(summary_SNORD115cov, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/allSNORD115_genomicFeatures_calctable.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  #Modify to plot barplot
summary_SNORD115cov_bar <- summary_SNORD115cov
summary_SNORD115cov_bar$bar_category <- "Exon"
summary_SNORD115cov_bar <- summary_SNORD115cov_bar %>% mutate(bar_category = ifelse(category %in% "Intron", "Intron", bar_category))
summary_SNORD115cov_bar$dataset <- "SNORD115all_vs_sharedGenes"
  #Remove unnecessary column
summary_SNORD115cov_bar <- summary_SNORD115cov_bar[-c(1)]
#For other chr15 SNORD copies coverage on shared genes
  #Rename columns
colnames(otherchr15SNORDs_coverage) <- c("orig_category", "basepairs")
  #Create column for category
otherchr15SNORDs_coverage$category <- c("3'UTR_only", "5'UTR_only", "CDS_only", "5'UTR+CDS", "Junction", "Intron", "Junction")
    #To save
  write.table(otherchr15SNORDs_coverage, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/otherchr15SNORDs_genomicFeatures_tablecategories.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  #Group by new category
otherchr15SNORDs_coverage_cat <- otherchr15SNORDs_coverage %>% group_by(category) %>% summarize(basepairs = sum(basepairs))
  #Calculate percentage of exon, intron, and intron-exon junction of all other chr15 SNORD binding events
summary_otherchr15SNORDs <- otherchr15SNORDs_coverage_cat
summary_otherchr15SNORDs$percent <- (summary_otherchr15SNORDs$basepairs/(sum(summary_otherchr15SNORDs$basepairs))*100)
    #To save
  write.table(summary_otherchr15SNORDs, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/otherchr15SNORDs_genomicFeatures_calctable.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
  #Modify to plot barplot
summary_otherchr15SNORDs_bar <- summary_otherchr15SNORDs
summary_otherchr15SNORDs_bar$bar_category <- "Exon"
summary_otherchr15SNORDs_bar <- summary_otherchr15SNORDs_bar %>% mutate(bar_category = ifelse(category %in% "Intron", "Intron", bar_category))
summary_otherchr15SNORDs_bar <- summary_otherchr15SNORDs_bar %>% mutate(bar_category = ifelse(category %in% "Junction", "Junction", bar_category))
summary_otherchr15SNORDs_bar$dataset <- "otherchr15SNORDs_vs_sharedGenes"
#Plotting
  #Make barplots
barplot_df <- rbind(summary_bkgrndcov_bar, summary_bkgrndcovRandom_bar, summary_bkgrndcovShared_bar, summary_SNORD116sharedALL_bar, summary_SNORD116groupI_shared_bar, summary_SNORD116groupII_shared_bar, summary_SNORD116groupIII_shared_bar, summary_SNORD116randomALL_bar, summary_SNORD116groupI_random_bar, summary_SNORD116groupII_random_bar, summary_SNORD116groupIII_random_bar)
barplot <- ggplot(barplot_df, aes(fill=bar_category, y=percent, x=dataset)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = c(Exon = "#8CD04A", Intron="#9C9C9C", Junction="#76094A")) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust=1, size=7), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_discrete(limits = c("Background: iNeuron Transcriptome", "Background: Shared Genes", "Background: Random Permutations", "SNORD116all_vs_sharedGenes", "SNORD116-I_vs_sharedGenes", "SNORD116-II_vs_sharedGenes", "SNORD116-III_vs_sharedGenes", "SNORD116all_vs_randomPerm", "SNORD116-I_vs_randomPerm", "SNORD116-II_vs_randomPerm", "SNORD116-III_vs_randomPerm"))
barplot
barplot_IandIII_df <- rbind(summary_SNORD116groupI_III_shared_bar, summary_SNORD116groupI_III_random_bar)
barplot_IandIII <- ggplot(barplot_IandIII_df, aes(fill=bar_category, y=percent, x=dataset)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = c(Exon = "#8CD04A", Intron="#9C9C9C", Junction="#76094A")) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust=1, size=7), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_discrete(limits = c("SNORD116-IandIII_vs_sharedGenes", "SNORD116-IandIII_vs_randomPerm"))
barplot_IandIII
barplot_SNORD115_other_df <- rbind(summary_SNORD115cov_bar, summary_otherchr15SNORDs_bar)
barplot_SNORD115_other <- ggplot(barplot_SNORD115_other_df, aes(fill=bar_category, y=percent, x=dataset)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = c(Exon = "#8CD04A", Intron="#9C9C9C", Junction="#76094A")) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust=1, size=7), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_discrete(limits = c("SNORD115all_vs_sharedGenes", "otherchr15SNORDs_vs_sharedGenes"))
barplot_SNORD115_other
  #Load required libraries for donut plotting
library("ggiraph")
library("plyr")
library("ggiraphExtra")
  #Parse for exon only
donutplot_df <- barplot_df[which(barplot_df$bar_category == "Exon"),]
donutplot_IandIII_df <- barplot_IandIII_df[which(barplot_IandIII_df$bar_category == "Exon"),]
donutplot_SNORD115_other_df <- barplot_SNORD115_other_df[which(barplot_SNORD115_other_df$bar_category == "Exon"),]
  #Parse by dataset
donut_bckground <- donutplot_df[which(donutplot_df$dataset == "Background: iNeuron Transcriptome"),]
donut_bckgroundShared <- donutplot_df[which(donutplot_df$dataset == "Background: Shared Genes"),]
donut_bckgroundRandom <- donutplot_df[which(donutplot_df$dataset == "Background: Random Permutations"),]
donut_116all_shared <- donutplot_df[which(donutplot_df$dataset == "SNORD116all_vs_sharedGenes"),]
donut_116groupI_shared <- donutplot_df[which(donutplot_df$dataset == "SNORD116-I_vs_sharedGenes"),]
donut_116groupII_shared <- donutplot_df[which(donutplot_df$dataset == "SNORD116-II_vs_sharedGenes"),]
donut_116groupIII_shared <- donutplot_df[which(donutplot_df$dataset == "SNORD116-III_vs_sharedGenes"),]
donut_116all_random <- donutplot_df[which(donutplot_df$dataset == "SNORD116all_vs_randomPerm"),]
donut_116groupI_random <- donutplot_df[which(donutplot_df$dataset == "SNORD116-I_vs_randomPerm"),]
donut_116groupII_random <- donutplot_df[which(donutplot_df$dataset == "SNORD116-II_vs_randomPerm"),]
donut_116groupIII_random <- donutplot_df[which(donutplot_df$dataset == "SNORD116-III_vs_randomPerm"),]
donut_116groupI_III_shared <- donutplot_IandIII_df[which(donutplot_IandIII_df$dataset == "SNORD116-IandIII_vs_sharedGenes"),]
donut_116groupI_III_random <- donutplot_IandIII_df[which(donutplot_IandIII_df$dataset == "SNORD116-IandIII_vs_randomPerm"),]
donut_SNORD115_shared <- donutplot_SNORD115_other_df[which(donutplot_SNORD115_other_df$dataset == "SNORD115all_vs_sharedGenes"),]
donut_otherchr15_shared <- donutplot_SNORD115_other_df[which(donutplot_SNORD115_other_df$dataset == "otherchr15SNORDs_vs_sharedGenes"),]
  #Calculate new percent based on exons only (this was just done to double check, will not actually be used to plot)
donut_bckground$percent <- (donut_bckground$basepairs/(sum(donut_bckground$basepairs))*100)
donut_bckgroundShared$percent <- (donut_bckgroundShared$basepairs/(sum(donut_bckgroundShared$basepairs))*100)
donut_bckgroundRandom$percent <- (donut_bckgroundRandom$basepairs/(sum(donut_bckgroundRandom$basepairs))*100)
donut_116all_shared$percent <- (donut_116all_shared$basepairs/(sum(donut_116all_shared$basepairs))*100)
donut_116groupI_shared$percent <- (donut_116groupI_shared$basepairs/(sum(donut_116groupI_shared$basepairs))*100)
donut_116groupII_shared$percent <- (donut_116groupII_shared$basepairs/(sum(donut_116groupII_shared$basepairs))*100)
donut_116groupIII_shared$percent <- (donut_116groupIII_shared$basepairs/(sum(donut_116groupIII_shared$basepairs))*100)
donut_116all_random$percent <- (donut_116all_random$basepairs/(sum(donut_116all_random$basepairs))*100)
donut_116groupI_random$percent <- (donut_116groupI_random$basepairs/(sum(donut_116groupI_random$basepairs))*100)
donut_116groupII_random$percent <- (donut_116groupII_random$basepairs/(sum(donut_116groupII_random$basepairs))*100)
donut_116groupIII_random$percent <- (donut_116groupIII_random$basepairs/(sum(donut_116groupIII_random$basepairs))*100)
donut_116groupI_III_shared$percent <- (donut_116groupI_III_shared$basepairs/(sum(donut_116groupI_III_shared$basepairs))*100)
donut_116groupI_III_random$percent <- (donut_116groupI_III_random$basepairs/(sum(donut_116groupI_III_random$basepairs))*100)
donut_SNORD115_shared$percent <- (donut_SNORD115_shared$basepairs/(sum(donut_SNORD115_shared$basepairs))*100)
donut_SNORD115_shared <- donut_SNORD115_shared %>% arrange(category)
donut_otherchr15_shared$percent <- (donut_otherchr15_shared$basepairs/(sum(donut_otherchr15_shared$basepairs))*100)
  #Plot donuts
plotdonut_bckground <- ggDonut(donut_bckground, aes(donuts=category,count=basepairs), title = "Background: iNeuron Transcriptome", labelsize = 2)
plotdonut_bckground
plotdonut_bckgroundShared <- ggDonut(donut_bckgroundShared, aes(donuts=category,count=basepairs), title = "Background: Shared Genes", labelsize = 2)
plotdonut_bckgroundShared
plotdonut_bckgroundRandom <- ggDonut(donut_bckgroundRandom, aes(donuts=category,count=basepairs), title = "Background: Random Genes", labelsize = 2)
plotdonut_bckgroundRandom
plotdonut_116all_shared <- ggDonut(donut_116all_shared, aes(donuts=category,count=basepairs), title = "SNORD116all vs sharedGenes", labelsize = 2)
plotdonut_116all_shared
plotdonut_116groupI_shared <- ggDonut(donut_116groupI_shared, aes(donuts=category,count=basepairs), title = "SNORD116-I vs sharedGenes", labelsize = 2)
plotdonut_116groupI_shared
plotdonut_116groupII_shared <- ggDonut(donut_116groupII_shared, aes(donuts=category,count=basepairs), title = "SNORD116-II vs sharedGenes", labelsize = 2)
plotdonut_116groupII_shared
plotdonut_116groupIII_shared <- ggDonut(donut_116groupIII_shared, aes(donuts=category,count=basepairs), title = "SNORD116-III vs sharedGenes", labelsize = 2)
plotdonut_116groupIII_shared
plotdonut_116all_random <- ggDonut(donut_116all_random, aes(donuts=category,count=basepairs), title = "SNORD116all vs randomPerm", labelsize = 2)
plotdonut_116all_random
plotdonut_116groupI_random <- ggDonut(donut_116groupI_random, aes(donuts=category,count=basepairs), title = "SNORD116-I vs randomPerm", labelsize = 2)
plotdonut_116groupI_random
plotdonut_116groupII_random <- ggDonut(donut_116groupII_random, aes(donuts=category,count=basepairs), title = "SNORD116-II vs randomPerm", labelsize = 2)
plotdonut_116groupII_random
plotdonut_116groupIII_random <- ggDonut(donut_116groupIII_random, aes(donuts=category,count=basepairs), title = "SNORD116-III vs randomPerm", labelsize = 2)
plotdonut_116groupIII_random
plotdonut_116groupI_III_shared <- ggDonut(donut_116groupI_III_shared, aes(donuts=category,count=basepairs), title = "SNORD116-I and III vs sharedGenes", labelsize = 2)
plotdonut_116groupI_III_shared
plotdonut_116groupI_III_random <- ggDonut(donut_116groupI_III_random, aes(donuts=category,count=basepairs), title = "SNORD116-I and III vs randomPerm", labelsize = 2)
plotdonut_116groupI_III_random
plotdonut_SNORD115_shared <- ggDonut(donut_SNORD115_shared, aes(donuts=category,count=basepairs), title = "SNORD115 vs sharedGenes", labelsize = 2)
plotdonut_SNORD115_shared
plotdonut_otherchr15_shared <- ggDonut(donut_otherchr15_shared, aes(donuts=category,count=basepairs), title = "other chr15 SNORDs vs sharedGenes", labelsize = 2)
plotdonut_otherchr15_shared
  #To arrange plots
alldonut <- ggarrange(plotdonut_bckground, plotdonut_bckgroundShared, plotdonut_bckgroundRandom, plotdonut_116all_shared, plotdonut_116all_random, plotdonut_116groupI_shared, plotdonut_116groupII_shared, plotdonut_116groupIII_shared, plotdonut_116groupI_random, plotdonut_116groupII_random, plotdonut_116groupIII_random,
                     ncol = 4, nrow = 3)
alldonut
plotsIandIII <- ggarrange(barplot_IandIII, plotdonut_116groupI_III_shared, plotdonut_116groupI_III_random,
                          ncol = 3, nrow = 1)
plotsIandIII
plots115andother <- ggarrange(barplot_SNORD115_other, plotdonut_SNORD115_shared, plotdonut_otherchr15_shared,
                             ncol = 3, nrow = 1)
plots115andother
#Save as pdfs
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/genomicFeatureCov_plots.pdf", width = 10, height = 8)
barplot
alldonut
dev.off()
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/genomicFeatureCov_IandIII_plots.pdf", width = 10, height = 8)
plotsIandIII
dev.off()
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/genomicFeatureCov_SNORD115andother_plots.pdf", width = 10, height = 8)
plots115andother
dev.off()

##Distribution of the predicted region of interaction of snoRNAs
#Import file
SNORDseq <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/chr15allSNORDs_seq.txt", sep="\t", head=F))
#Rename columns
colnames(SNORDseq) <- c("snoENSEMBL_ID", "snoRNA_number", "sequence", "length")
#Creating plot of number of binding events
  #Merge files
dist115 <- merge(targets115name, SNORDseq, by="snoENSEMBL_ID", sort = TRUE, all = FALSE)
dist115 <- merge(dist115, genes.table_115, by="target_id", sort = TRUE, all = FALSE)
dist116 <- merge(targets116name, SNORDseq, by="snoENSEMBL_ID", sort = TRUE, all = FALSE)
dist116 <- merge(dist116, genes.table_116, by="target_id", sort = TRUE, all = FALSE)
dist_other <- merge(targets_other2, SNORDseq, by="snoENSEMBL_ID", sort = TRUE, all = FALSE)
dist_other <- merge(dist_other, genes.table_other, by="target_id", sort = TRUE, all = FALSE)
  #Calculate center of binding event
dist115$center_sno_bind <- ((dist115$sno_window_end - dist115$sno_window_start)/2) + dist115$sno_window_start
dist116$center_sno_bind <- ((dist116$sno_window_end - dist116$sno_window_start)/2) + dist116$sno_window_start
dist_other$center_sno_bind <- ((dist_other$sno_window_end - dist_other$sno_window_start)/2) + dist_other$sno_window_start
  #Make ratio column
dist115$rel_bind <- (dist115$center_sno_bind/dist115$length)
dist116$rel_bind <- (dist116$center_sno_bind/dist116$length)
dist_other$rel_bind <- (dist_other$center_sno_bind/dist_other$length)
  #Make binding event length column
dist115$bind_len <- (dist115$sno_window_end - dist115$sno_window_start)
dist116$bind_len <- (dist116$sno_window_end - dist116$sno_window_start)
dist_other$bind_len <- (dist_other$sno_window_end - dist_other$sno_window_start)
  #Create column for SNORD116 groups
dist116_group <- as.data.frame(dist116 %>% separate(snoRNA_number.x, c("SNORD116", "snoRNA_number"), sep = "-"))
dist116_group$group <- "116-I"
dist116_group <- dist116_group %>% mutate(group = ifelse(snoRNA_number %in% 10:24, "116-II", group))
dist116_group <- dist116_group %>% mutate(group = ifelse(snoRNA_number %in% 25:30, "116-III", group))
  #Subset by group
dist116_groupI <- dist116_group[which(dist116_group$group == "116-I"),]
dist116_groupII <- dist116_group[which(dist116_group$group == "116-II"),]
dist116_groupIII <- dist116_group[which(dist116_group$group == "116-III"),]
dist116_groupI_III <- dist116_group %>% filter(!dist116_group$group %in% "116-II")
  #Group & count occurances
grouped_dist115 <- dist115 %>% group_by(rel_bind)
summary_dist115 <- grouped_dist115 %>% summarize(count = n())
grouped_dist116 <- dist116 %>% group_by(rel_bind)
summary_dist116 <- grouped_dist116 %>% summarize(count = n())
grouped_distother <- dist_other %>% group_by(rel_bind)
summary_distother <- grouped_distother %>% summarize(count = n())
grouped_dist116_groupI <- dist116_groupI %>% group_by(rel_bind)
summary_dist116_groupI <- grouped_dist116_groupI %>% summarize(count = n())
grouped_dist116_groupII <- dist116_groupII %>% group_by(rel_bind)
summary_dist116_groupII <- grouped_dist116_groupII %>% summarize(count = n())
grouped_dist116_groupIII <- dist116_groupIII %>% group_by(rel_bind)
summary_dist116_groupIII <- grouped_dist116_groupIII %>% summarize(count = n())
grouped_dist116_groupI_III <- dist116_groupI_III %>% group_by(rel_bind)
summary_dist116_groupI_III <- grouped_dist116_groupI_III %>% summarize(count = n())
  #Plot
plot115dist <- ggplot(summary_dist115, aes(x=rel_bind)) + geom_line(aes(y=count)) + 
  ggtitle("All SNORD115 Copies", subtitle = "Distribution of Predicted Interactions") + scale_x_continuous(name = "Relative position in snoRNA", breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0, 1)) + ylab("# of predicted interactions") +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plot115dist
plot116dist <- ggplot(summary_dist116, aes(x=rel_bind)) + geom_line(aes(y=count)) +
  ggtitle("All SNORD116 Copies", subtitle = "Distribution of Predicted Interactions") + scale_x_continuous(name = "Relative position in snoRNA", breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0, 1)) + ylab("# of predicted interactions") +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plot116dist
plototherdist <- ggplot(summary_distother, aes(x=rel_bind)) + geom_line(aes(y=count)) +
  ggtitle("Other SNORD Copies on chr15", subtitle = "Distribution of Predicted Interactions") + scale_x_continuous(name = "Relative position in snoRNA", breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0, 1)) + ylab("# of predicted interactions") +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plototherdist
plot116dist_groupI <- ggplot(summary_dist116_groupI, aes(x=rel_bind)) + geom_line(aes(y=count)) +
  ggtitle("SNORD116-I", subtitle = "Distribution of Predicted Interactions") + scale_x_continuous(name = "Relative position in snoRNA", breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0, 1)) + ylab("# of predicted interactions") +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plot116dist_groupI
plot116dist_groupII <- ggplot(summary_dist116_groupII, aes(x=rel_bind)) + geom_line(aes(y=count)) +
  ggtitle("SNORD116-II", subtitle = "Distribution of Predicted Interactions") + scale_x_continuous(name = "Relative position in snoRNA", breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0, 1)) + ylab("# of predicted interactions") +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plot116dist_groupII
plot116dist_groupIII <- ggplot(summary_dist116_groupIII, aes(x=rel_bind)) + geom_line(aes(y=count)) +
  ggtitle("SNORD116-III", subtitle = "Distribution of Predicted Interactions") + scale_x_continuous(name = "Relative position in snoRNA", breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0, 1)) + ylab("# of predicted interactions") +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plot116dist_groupIII
plot116dist_groupI_III <- ggplot(summary_dist116_groupI_III, aes(x=rel_bind)) + geom_line(aes(y=count)) +
  ggtitle("SNORD116-I & III", subtitle = "Distribution of Predicted Interactions") + scale_x_continuous(name = "Relative position in snoRNA", breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0, 1)) + ylab("# of predicted interactions") +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plot116dist_groupI_III
  #Arrange plots
alldist <- ggarrange(plot115dist, plot116dist, plototherdist, plot116dist_groupI, plot116dist_groupII, plot116dist_groupIII, plot116dist_groupI_III,
                     ncol = 4, nrow = 2)
alldist
  #Save as pdf
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/alldistacrosssnoRNA_plots.pdf", width = 10, height = 8)
alldist
plot116dist
plot115dist
plototherdist
plot116dist_groupI
plot116dist_groupII
plot116dist_groupIII
plot116dist_groupI_III
dev.off()
#Creating heatmap for location of snoRNA elements
  #Parse sequence file
snord115seq <- SNORDseq[grepl("SNORD115", SNORDseq$snoRNA_number), ]
snord116seq <- SNORDseq[grepl("SNORD116", SNORDseq$snoRNA_number), ]
otherSNORDseq <- SNORDseq[1:5,]
  #Analyze SNORD116
    #Run loop to parse for C-boxes
library("stringr")
snord116_cbox <- list()
query <- "ATGATGA|ATGACGA|ACGATGA"
for (i in 1:30)
{
  x <- snord116seq[i,3]
  snord116_cbox[[i]] <- as.data.frame(str_locate_all(x, query))
  snord116_cbox[[i]]$length <- nchar(x)
  snord116_cbox[[i]]$rel_start <- snord116_cbox[[i]]$start/snord116_cbox[[i]]$length
  snord116_cbox[[i]]$rel_end <- snord116_cbox[[i]]$end/snord116_cbox[[i]]$length
  snord116_cbox[[i]]$snoRNA_number <- snord116seq[i,2]
}
snord116_cbox_all <- bind_rows(snord116_cbox)
    #Run loop to parse for C'-boxes
snord116_cprimebox <- list()
query <- "ATGAGTG|ACGAGTG|ATGAATG"
for (i in 1:30)
{
  x <- snord116seq[i,3]
  snord116_cprimebox[[i]] <- as.data.frame(str_locate_all(x, query))
  snord116_cprimebox[[i]]$length <- nchar(x)
  snord116_cprimebox[[i]]$rel_start <- snord116_cprimebox[[i]]$start/snord116_cprimebox[[i]]$length
  snord116_cprimebox[[i]]$rel_end <- snord116_cprimebox[[i]]$end/snord116_cprimebox[[i]]$length
  snord116_cprimebox[[i]]$snoRNA_number <- snord116seq[i,2]
}
snord116_cprimebox_all <- bind_rows(snord116_cprimebox)
    #Run loop to parse for D'-boxes
snord116_dprimebox <- list()
query <- "AAGCTGA|AAGGTGA|AAGTTGA|AATCTGA"
for (i in 1:30)
{
  x <- snord116seq[i,3]
  snord116_dprimebox[[i]] <- as.data.frame(str_locate_all(x, query))
  snord116_dprimebox[[i]]$d_start <- snord116_dprimebox[[i]]$start + 3
  snord116_dprimebox[[i]]$length <- nchar(x)
  snord116_dprimebox[[i]]$rel_start <- snord116_dprimebox[[i]]$d_start/snord116_dprimebox[[i]]$length
  snord116_dprimebox[[i]]$rel_end <- snord116_dprimebox[[i]]$end/snord116_dprimebox[[i]]$length
  snord116_dprimebox[[i]]$snoRNA_number <- snord116seq[i,2]
}
snord116_dprimebox_all <- bind_rows(snord116_dprimebox)
    #Run loop to parse for D-boxes
snord116_dbox <- list()
query <- "AACTGAGG|aactgagg|AGCTGAGG"
for (i in 1:30)
{
  x <- snord116seq[i,3]
  snord116_dbox[[i]] <- as.data.frame(str_locate_all(x, query))
  snord116_dbox[[i]]$d_start <- snord116_dbox[[i]]$start + 2
  snord116_dbox[[i]]$d_end <- snord116_dbox[[i]]$end - 2
  snord116_dbox[[i]]$length <- nchar(x)
  snord116_dbox[[i]]$rel_start <- snord116_dbox[[i]]$d_start/snord116_dbox[[i]]$length
  snord116_dbox[[i]]$rel_end <- snord116_dbox[[i]]$d_end/snord116_dbox[[i]]$length
  snord116_dbox[[i]]$snoRNA_number <- snord116seq[i,2]
}
snord116_dbox_all <- bind_rows(snord116_dbox)
  #Analyze SNORD115
    #Run loop to parse for C-boxes
snord115_cbox <- list()
query <- "GATGATGAGAA|AATGATGAGAA|AATGAGAACCT|AGTGTCGAGAA|GATTATGAGAA|AGTGATGAGAA|GGTGATGAGAA|AATATGGAGAT|AATAATGAAAT|AACGTAAATAT|AATGATGAGAT"
for (i in 1:48)
{
  x <- snord115seq[i,3]
  snord115_cbox[[i]] <- as.data.frame(str_locate_all(x, query))
  snord115_cbox[[i]]$c_start <- snord115_cbox[[i]]$start + 1
  snord115_cbox[[i]]$c_end <- snord115_cbox[[i]]$end - 3
  snord115_cbox[[i]]$length <- nchar(x)
  snord115_cbox[[i]]$rel_start <- snord115_cbox[[i]]$c_start/snord115_cbox[[i]]$length
  snord115_cbox[[i]]$rel_end <- snord115_cbox[[i]]$c_end/snord115_cbox[[i]]$length
  snord115_cbox[[i]]$snoRNA_number <- snord115seq[i,2]
}
snord115_cbox_all <- bind_rows(snord115_cbox)
    #Run loop to parse for C'-boxes
snord115_cprimebox <- list()
query <- "GGTGATGACTT|GGTGATAACTT|GATGATGACTT|GGTGGTGACTT|GGTGATTATTT|GAAGATGACAT|AATAATGCCTT|CGTAATAATGT|AATGATGACGT"
for (i in 1:48)
{
  x <- snord115seq[i,3]
  snord115_cprimebox[[i]] <- as.data.frame(str_locate_all(x, query))
  snord115_cprimebox[[i]]$c_start <- snord115_cprimebox[[i]]$start + 1
  snord115_cprimebox[[i]]$c_end <- snord115_cprimebox[[i]]$end - 3
  snord115_cprimebox[[i]]$length <- nchar(x)
  snord115_cprimebox[[i]]$rel_start <- snord115_cprimebox[[i]]$c_start/snord115_cprimebox[[i]]$length
  snord115_cprimebox[[i]]$rel_end <- snord115_cprimebox[[i]]$c_end/snord115_cprimebox[[i]]$length
  snord115_cprimebox[[i]]$snoRNA_number <- snord115seq[i,2]
}
snord115_cprimebox_all <- bind_rows(snord115_cprimebox)
    #Run loop to parse for D'-boxes
snord115_dprimebox <- list()
query <- "CCTGAAG|CTTGAAG|TCTGAAG|GTTGAAG|CCTGAAA|TTCGACA|TCTGATT|TCAGGAG"
for (i in 1:48)
{
  x <- snord115seq[i,3]
  snord115_dprimebox[[i]] <- as.data.frame(str_locate_all(x, query))
  snord115_dprimebox[[i]]$d_start <- snord115_dprimebox[[i]]$start + 1
  snord115_dprimebox[[i]]$d_end <- snord115_dprimebox[[i]]$end - 2
  snord115_dprimebox[[i]]$length <- nchar(x)
  snord115_dprimebox[[i]]$rel_start <- snord115_dprimebox[[i]]$d_start/snord115_dprimebox[[i]]$length
  snord115_dprimebox[[i]]$rel_end <- snord115_dprimebox[[i]]$d_end/snord115_dprimebox[[i]]$length
  snord115_dprimebox[[i]]$snoRNA_number <- snord115seq[i,2]
}
snord115_dprimebox_all <- bind_rows(snord115_dprimebox)
    #Run loop to parse for D-boxes
snord115_dbox <- list()
query <- "GCTGAGG|GCTGAGT|ATTGAAG|ACTTAGG|ACTGAGG|GCTGAAG|GTGGAGA"
for (i in 1:48)
{
  x <- snord115seq[i,3]
  snord115_dbox[[i]] <- as.data.frame(str_locate_all(x, query))
  snord115_dbox[[i]]$d_start <- snord115_dbox[[i]]$start + 1
  snord115_dbox[[i]]$d_end <- snord115_dbox[[i]]$end - 2
  snord115_dbox[[i]]$length <- nchar(x)
  snord115_dbox[[i]]$rel_start <- snord115_dbox[[i]]$d_start/snord115_dbox[[i]]$length
  snord115_dbox[[i]]$rel_end <- snord115_dbox[[i]]$d_end/snord115_dbox[[i]]$length
  snord115_dbox[[i]]$snoRNA_number <- snord115seq[i,2]
}
snord115_dbox_all <- bind_rows(snord115_dbox)
  #Analyze other SNORDs
    #Run loop to parse for C-boxes
snordother_cbox <- list()
query <- "GTGATGAG|ATGATGAC|ATGATGAG"
for (i in 1:5)
{
  x <- otherSNORDseq[i,3]
  snordother_cbox[[i]] <- as.data.frame(str_locate_all(x, query))
  snordother_cbox[[i]]$c_end <- snordother_cbox[[i]]$end - 1
  snordother_cbox[[i]]$length <- nchar(x)
  snordother_cbox[[i]]$rel_start <- snordother_cbox[[i]]$start/snordother_cbox[[i]]$length
  snordother_cbox[[i]]$rel_end <- snordother_cbox[[i]]$c_end/snordother_cbox[[i]]$length
  snordother_cbox[[i]]$snoRNA_number <- otherSNORDseq[i,2]
}
snordother_cbox_all <- bind_rows(snordother_cbox)
    #Run loop to parse for C'-boxes
snordother_cprimebox <- list()
query <- "ATGATGAA|ATAATGAT|TGGATGAC|ATGCTGAG"
for (i in 1:5)
{
  x <- otherSNORDseq[i,3]
  snordother_cprimebox[[i]] <- as.data.frame(str_locate_all(x, query))
  snordother_cprimebox[[i]]$c_end <- snordother_cprimebox[[i]]$end - 1
  snordother_cprimebox[[i]]$length <- nchar(x)
  snordother_cprimebox[[i]]$rel_start <- snordother_cprimebox[[i]]$start/snordother_cprimebox[[i]]$length
  snordother_cprimebox[[i]]$rel_end <- snordother_cprimebox[[i]]$c_end/snordother_cprimebox[[i]]$length
  snordother_cprimebox[[i]]$snoRNA_number <- otherSNORDseq[i,2]
}
snordother_cprimebox_all <- bind_rows(snordother_cprimebox)
    #Run loop to parse for D & D'-boxes
snordother_dboxes <- list()
query <- "ACTGAG|TCTGAA|CTTGAA|CGTGAG|TCTGAG"
for (i in 1:5)
{
  x <- otherSNORDseq[i,3]
  snordother_dboxes[[i]] <- as.data.frame(str_locate_all(x, query))
  snordother_dboxes[[i]]$d_start <- snordother_dboxes[[i]]$start + 1
  snordother_dboxes[[i]]$d_end <- snordother_dboxes[[i]]$end - 1
  snordother_dboxes[[i]]$length <- nchar(x)
  snordother_dboxes[[i]]$rel_start <- snordother_dboxes[[i]]$d_start/snordother_dboxes[[i]]$length
  snordother_dboxes[[i]]$rel_end <- snordother_dboxes[[i]]$d_end/snordother_dboxes[[i]]$length
  snordother_dboxes[[i]]$snoRNA_number <- otherSNORDseq[i,2]
}
snordother_dboxes_all <- bind_rows(snordother_dboxes)
  #Parse groups
groupI_cbox <- snord116_cbox_all[1:9,]
groupI_cprimebox <- snord116_cprimebox_all[1:9,]
groupI_dbox <- snord116_dbox_all[1:9,]
groupI_dprimebox <- snord116_dprimebox_all[1:9,]
groupII_cbox <- snord116_cbox_all[10:24,]
groupII_cprimebox <- snord116_cprimebox_all[10:24,]
groupII_dbox <- snord116_dbox_all[10:24,]
groupII_dprimebox <- snord116_dprimebox_all[10:24,]
groupIII_cbox <- snord116_cbox_all[25:30,]
groupIII_cprimebox <- snord116_cprimebox_all[25:30,]
groupIII_dbox <- snord116_dbox_all[25:30,]
groupIII_dprimebox <- snord116_dprimebox_all[25:30,]
  #Graph
library("grid")
pdf('../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/snord_elements.pdf',w=10, h=1)
grid.newpage()
for (rel_start in snord116_cbox_all)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            # 1A is 90% tranparency: https://stackoverflow.com/questions/23201134/transparent-argb-hex-value
            #width may need to be a variable that is the length of the k-mer
            width = unit((7/snord116_cbox_all$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#008C001A'))
}
for (rel_start in snord116_cprimebox_all)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((7/snord116_cprimebox_all$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#008C001A'))
}
for (rel_start in snord116_dprimebox_all)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((4/snord116_dprimebox_all$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#4B00921A'))
}
for (rel_start in snord116_dbox_all)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((4/snord116_dbox_all$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#4B00921A'))
}
grid.newpage()
for (rel_start in snord115_cbox_all)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((7/snord115_cbox_all$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#008C001A'))
}
for (rel_start in snord115_cprimebox_all)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((7/snord115_cprimebox_all$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#008C001A'))
}
for (rel_start in snord115_dprimebox_all)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((4/snord115_dprimebox_all$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#4B00921A'))
}
for (rel_start in snord115_dbox_all)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((4/snord115_dbox_all$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#4B00921A'))
}
grid.newpage()
for (rel_start in snordother_cbox_all)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((7/snordother_cbox_all$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#008C001A'))
}
for (rel_start in snordother_cprimebox_all)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((7/snordother_cprimebox_all$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#008C001A'))
}
for (rel_start in snordother_dboxes_all)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((4/snordother_dboxes_all$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#4B00921A'))
}
dev.off()
pdf('../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/snord116_elements_bygroup.pdf',w=10, h=1)
grid.newpage()
for (rel_start in groupI_cbox)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((7/groupI_cbox$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#008C001A'))
}
for (rel_start in groupI_cprimebox)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((7/groupI_cprimebox$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#008C001A'))
}
for (rel_start in groupI_dprimebox)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((4/groupI_dprimebox$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#4B00921A'))
}
for (rel_start in groupI_dbox)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((4/groupI_dbox$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#4B00921A'))
}
grid.newpage()
for (rel_start in groupII_cbox)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((7/groupII_cbox$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#008C001A'))
}
for (rel_start in groupII_cprimebox)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((7/groupII_cprimebox$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#008C001A'))
}
for (rel_start in groupII_dprimebox)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((4/groupII_dprimebox$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#4B00921A'))
}
for (rel_start in groupII_dbox)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((4/groupII_dbox$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#4B00921A'))
}
grid.newpage()
for (rel_start in groupIII_cbox)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((7/groupIII_cbox$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#008C001A'))
}
for (rel_start in groupIII_cprimebox)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((7/groupIII_cprimebox$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#008C001A'))
}
for (rel_start in groupIII_dprimebox)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"),
            width = unit((4/groupIII_dprimebox$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#4B00921A'))
}
for (rel_start in groupIII_dbox)
{
  grid.rect(x = unit(rel_start, "npc"), y = unit(0.5, "npc"), 
            width = unit((4/groupIII_dbox$length), "npc"), height = unit(1, "npc"), just = 'left', gp = gpar(col = NA, fill = '#4B00921A'))
}
dev.off()

##Comparison of SNORDs known to target rRNAs vs list of rRNAs & our shared gene list
#Import files
rRNAtargets116 <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/snord116_rRNA", sep="\t", head=T))
rRNAtargets_knownrRNASNORDs <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/rRNASNORDs_rRNA", sep="\t", head=T))
targets_knownrRNASNORDs <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/rRNASNORDs_degshared", sep="\t", head=T))
#Remove non-chr15 SNORD116 from file
rRNAtargets116 <- rRNAtargets116[-(2),]
#Change ESEMBL ID's to external gene name
mart <- useMart("ensembl", host = "https://apr2018.archive.ensembl.org/", dataset = "hsapiens_gene_ensembl")
  #SNORD116 on rRNAs
R1 <- data.table(rRNAtargets116$target_id)
colnames(R1) <- c("target_id")
genes.table_1 <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"), values = R1$target_id, mart = mart)
names(genes.table_1)[1] <- "target_id"
  #rRNA SNORDs on rRNAs
R2 <- data.table(rRNAtargets_knownrRNASNORDs$target_id)
colnames(R2) <- c("target_id")
genes.table_2 <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"), values = R2$target_id, mart = mart)
names(genes.table_2)[1] <- "target_id"
  #rRNA SNORDs on shared genes
R3 <- data.table(targets_knownrRNASNORDs$target_id)
colnames(R3) <- c("target_id")
genes.table_3 <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"), values = R3$target_id, mart = mart)
names(genes.table_3)[1] <- "target_id"
#Parse by gene target & plot
  #SNORD116 on rRNAs
rRNAcount116 <- as.data.frame(rRNAtargets116 %>% count(target_id))
rRNAcount116name <- merge(genes.table_1, rRNAcount116, by="target_id", sort = TRUE, all = TRUE)
rRNAplot116 <- ggplot(data = rRNAcount116name, aes(x=external_gene_name, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD116: rRNA Targeting Events", x ="Target Gene Name", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5)
rRNAplot116
  #rRNA SNORDs on rRNAs
rRNAcountrRNASNORDs <- as.data.frame(rRNAtargets_knownrRNASNORDs %>% count(target_id))
rRNAcountrRNASNORDsname <- merge(genes.table_2, rRNAcountrRNASNORDs, by="target_id", sort = TRUE, all = TRUE)
rRNAplotrRNASNORDs <- ggplot(data = rRNAcountrRNASNORDsname, aes(x=external_gene_name, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="Known rRNA SNORDs: rRNA Targeting Events", x ="Target Gene Name", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5)
rRNAplotrRNASNORDs
  #rRNA SNORDs on shared genes
genecountrRNASNORDs <- as.data.frame(targets_knownrRNASNORDs %>% count(target_id))
genecountrRNASNORDsname <- merge(genes.table_3, genecountrRNASNORDs, by="target_id", sort = TRUE, all = TRUE)
geneplotrRNASNORDs <- ggplot(data = genecountrRNASNORDsname, aes(x=external_gene_name, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="Known rRNA SNORDs: Shared Gene Targeting Events", x ="Target Gene Name", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5)
geneplotrRNASNORDs
#Save as pdf
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/rRNAtargetingSNORDs_byGene.pdf", width = 10, height = 8)
rRNAplot116
rRNAplotrRNASNORDs
geneplotrRNASNORDs
dev.off()
#Parse by sno copy & plot
  #SNORD116 on rRNAs
rRNAsnocount116 <- as.data.frame(rRNAtargets116 %>% count(sno_id))
rRNAsnocount116name <- as.data.frame(rRNAsnocount116 %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
rRNAsnocount116name$snoRNA_number <- sapply(rRNAsnocount116name$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
rRNAsnoplot116 <- ggplot(data = rRNAsnocount116name, aes(x=snoRNA_number, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD116: rRNA Targeting Events by snoRNA Copy", x ="snoRNA Copy", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5)
rRNAsnoplot116
  #rRNA SNORDs on rRNAs
rRNAsnocountrRNASNORDs <- as.data.frame(rRNAtargets_knownrRNASNORDs %>% count(sno_id))
rRNAsnocountrRNASNORDsname <- as.data.frame(rRNAsnocountrRNASNORDs %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
rRNAsnocountrRNASNORDsname$snoRNA_number <- sapply(rRNAsnocountrRNASNORDsname$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
rRNAsnoplotrRNASNORDs <- ggplot(data = rRNAsnocountrRNASNORDsname, aes(x=snoRNA_number, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="Known rRNA SNORDs: rRNA Targeting Events by snoRNA Copy", x ="snoRNA Copy", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5)
rRNAsnoplotrRNASNORDs
  #rRNA SNORDs on shared genes
snocountrRNASNORDs <- as.data.frame(targets_knownrRNASNORDs %>% count(sno_id))
snocountrRNASNORDsname <- snocountrRNASNORDs
snocountrRNASNORDsname$sno_id <- sapply(snocountrRNASNORDsname$sno_id, function(x) sub("\\(.*\\)", "", x))
snoplotrRNASNORDs <- ggplot(data = snocountrRNASNORDsname, aes(x=sno_id, y=n)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="Known rRNA SNORDs: Shared Gene Targeting Events by snoRNA Copy", x ="snoRNA Copy", y = "Number of Targeting Events") + geom_text(
      aes(label=n), vjust=1.1, color="white", size=3.5)
snoplotrRNASNORDs
#Save as pdf
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/rRNAtargetingSNORDs_byCopy.pdf", width = 10, height = 8)
rRNAsnoplot116
rRNAsnoplotrRNASNORDs
snoplotrRNASNORDs
dev.off()
#Subset by sno copy & gene
  #SNORD116 on rRNAs
rRNAtargets116name <- as.data.frame(rRNAtargets116 %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
rRNAtargets116name$snoRNA_number <- sapply(rRNAtargets116name$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
grouped_rRNAtargets116 <- rRNAtargets116name %>% group_by(target_id, snoRNA_number)
summary_rRNAtargets116 <- grouped_rRNAtargets116 %>% summarize(count = n())
summary_rRNAtargets116name <- merge(genes.table_1, summary_rRNAtargets116, by="target_id", sort = TRUE, all = TRUE)
plotsummaryrRNAtargets116_A <- ggplot(data = summary_rRNAtargets116name, aes(x=snoRNA_number, y=count, fill=external_gene_name)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD116 Summary: rRNA Targeting Events by snoRNA Copy", x ="SNORD copy", y = "Number of Targeting Events")
plotsummaryrRNAtargets116_A
plotsummaryrRNAtargets116_B <- ggplot(data = summary_rRNAtargets116name, aes(x=external_gene_name, y=count, fill=snoRNA_number)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="SNORD116 Summary: rRNA Targeting Events by snoRNA Copy", x ="Target Gene Name", y = "Number of Targeting Events")
plotsummaryrRNAtargets116_B
  #rRNA SNORDs on rRNAs
rRNAtargetsrRNASNORDsname <- rRNAtargets_knownrRNASNORDs
rRNAtargetsrRNASNORDsname$sno_id <- sapply(rRNAtargetsrRNASNORDsname$sno_id, function(x) sub("\\(.*\\)", "", x))
grouped_rRNAtargetsrRNASNORDs <- rRNAtargetsrRNASNORDsname %>% group_by(target_id, sno_id)
summary_rRNAtargetsrRNASNORDs <- grouped_rRNAtargetsrRNASNORDs %>% summarize(count = n())
summary_rRNAtargetsrRNASNORDs <- merge(genes.table_2, summary_rRNAtargetsrRNASNORDs, by="target_id", sort = TRUE, all = TRUE)
plotsummaryrRNAtargetsrRNASNORDs_A <- ggplot(data = summary_rRNAtargetsrRNASNORDs, aes(x=sno_id, y=count, fill=external_gene_name)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="Known rRNA SNORDs Summary: rRNA Targeting Events by snoRNA Copy", x ="SNORD copy", y = "Number of Targeting Events")
plotsummaryrRNAtargetsrRNASNORDs_A
plotsummaryrRNAtargetsrRNASNORDs_B <- ggplot(data = summary_rRNAtargetsrRNASNORDs, aes(x=external_gene_name, y=count, fill=sno_id)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="Known rRNA SNORDs Summary: rRNA Targeting Events by snoRNA Copy", x ="Target Gene Name", y = "Number of Targeting Events")
plotsummaryrRNAtargetsrRNASNORDs_B
  #rRNA SNORDs on shared genes
targetsrRNASNORDsname <- targets_knownrRNASNORDs
targetsrRNASNORDsname$sno_id <- sapply(targetsrRNASNORDsname$sno_id, function(x) sub("\\(.*\\)", "", x))
grouped_targetsrRNASNORDs <- targetsrRNASNORDsname %>% group_by(target_id, sno_id)
summary_targetsrRNASNORDs <- grouped_targetsrRNASNORDs %>% summarize(count = n())
summary_targetsrRNASNORDs <- merge(genes.table_3, summary_targetsrRNASNORDs, by="target_id", sort = TRUE, all = TRUE)
plotsummarytargetsrRNASNORDs_A <- ggplot(data = summary_targetsrRNASNORDs, aes(x=sno_id, y=count, fill=external_gene_name)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="Known rRNA SNORDs Summary: Shared Gene Targeting Events by snoRNA Copy", x ="SNORD copy", y = "Number of Targeting Events")
plotsummarytargetsrRNASNORDs_A
plotsummarytargetsrRNASNORDs_B <- ggplot(data = summary_targetsrRNASNORDs, aes(x=external_gene_name, y=count, fill=sno_id)) + geom_bar(stat="identity") + theme(
  axis.text.x = element_text(angle = 60, hjust=1, size=7)) + labs(
    title="Known rRNA SNORDs Summary: Shared Gene Targeting Events by snoRNA Copy", x ="Target Gene Name", y = "Number of Targeting Events")
plotsummarytargetsrRNASNORDs_B
#Save as pdf
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/rRNAtargetingSNORDs_byCopyandGene.pdf", width = 10, height = 8)
plotsummaryrRNAtargets116_A
plotsummaryrRNAtargets116_B
plotsummaryrRNAtargetsrRNASNORDs_A
plotsummaryrRNAtargetsrRNASNORDs_B
plotsummarytargetsrRNASNORDs_A
plotsummarytargetsrRNASNORDs_B
dev.off()

##Test set with SNORD22 vs 3 genes from snoGloBe paper
#Import files
targets_test <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD22_testGenes", sep="\t", head=T))
#Multiply count column by 100 to generate a "score" for the bed file
bedtest <- targets_test
bedtest$score <- bedtest$count * 100
#Separate sno column
bedtest <- as.data.frame(bedtest %>% separate(sno_id, c("snoENSEMBL_ID", "snoRNA_number"), sep = "_"))
bedtest$snoRNA_number <- sapply(bedtest$snoRNA_number, function(x) sub("\\(.*\\)", "", x))
#Drop unnecessary columns & reorder
bedtest <- bedtest[c(1:3,8,14,6)]
#Write out file
write.table(bedtest, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/snord22_testtargets_snoglobe_hg38.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  #In notepad, edited file to include trackline & uploaded into browser (https://genome.ucsc.edu/s/rbgilmore/snoGloBe_targets_RBG)

##Does the differential expression value correlate with the number of targeting events?
#Read in results object
counts_smDEL <- read.csv("../Cotney_Lab/PWS_RNASeq/STAR/smDEL/CombinedH9CT2_Results_DESeqSTAR.csv", header = TRUE, row.names = "Gene")
#Subset log2foldChange for H9
fc <- as.data.frame(counts_smDEL[,c(1,3)])
#Calculate absolute value of the log2foldChange
fc$abslog2FC_H9 <- abs(fc$log2FoldChange_H9)
names(fc)[1] <- "target_id"
#Calculate stats
stats_115 <- summary_targets115name %>% group_by(target_id) %>% summarize(sum115 = sum(count), median115 = median(count), mean115 = mean(count))
stats_116_I <- groupI %>% group_by(target_id) %>% summarize(sum116_I = sum(count), median116_I = median(count), mean116_I = mean(count))
stats_116_II <- groupII %>% group_by(target_id) %>% summarize(sum116_II = sum(count), median116_II = median(count), mean116_II = mean(count))
stats_116_III <- groupIII %>% group_by(target_id) %>% summarize(sum116_III = sum(count), median116_III = median(count), mean116_III = mean(count))
#Merge all
stats_116_ALL <- merge(stats_116_I, stats_116_II, by="target_id", sort = TRUE, all = TRUE)
stats_116_ALL <- merge(stats_116_ALL, stats_116_III, by="target_id", sort = TRUE, all = TRUE)
#Merge with foldchange object
stats_116_ALL <- merge(fc, stats_116_ALL, by="target_id", sort = TRUE, all = FALSE)
#Add gene names
stats_116_ALL <- merge(genes.table_116, stats_116_ALL, by="target_id", sort = TRUE, all = FALSE)
#Merge with SNORD115
stats_115_116_ALL <- merge(stats_116_ALL, stats_115, by="target_id", sort = TRUE, all = TRUE)
#Subset "enriched" genes from above & remove SNHG14 (low expression)
statstoPlot <- stats_115_116_ALL %>% filter(stats_115_116_ALL$external_gene_name %in% enriched_list)
statstoPlot <- statstoPlot[-(19),]
#Plot
library("devtools")
plot_group1_mean <- ggplot(statstoPlot, aes(x=mean116_I, y=abslog2FC_H9)) +
  geom_point(color = "blue") + geom_smooth(method = lm, colour="blue", fill = "lightblue", alpha = 0.25) + stat_cor(aes(label = paste(..p.label..)), label.x = 0, label.y = 2.5) + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~")), label.x = 0, label.y = 3) + geom_text(label=statstoPlot$external_gene_name, colour = "darkblue") +
  ggtitle("SNORD116-I") + theme(plot.title = element_text(hjust = 0.5)) + labs(x = "Mean Counts", y = "abs(log2foldchange)")
plot_group2_mean <- ggplot(statstoPlot, aes(x=mean116_II, y=abslog2FC_H9)) +
  geom_point(color = "red") + geom_smooth(method = lm, colour="red", fill="pink", alpha = 0.35) + stat_cor(aes(label = paste(..p.label..)), label.x = 0, label.y = 2.5) + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~")), label.x = 0, label.y = 3) + geom_text(label=statstoPlot$external_gene_name, colour = "darkred") +
  ggtitle("SNORD116-II") + theme(plot.title = element_text(hjust = 0.5)) + labs(x = "Mean Counts", y = "abs(log2foldchange)")
plot_group3_mean <- ggplot(statstoPlot, aes(x=mean116_III, y=abslog2FC_H9)) +
  geom_point(color = "green") + geom_smooth(method = lm, colour="green", fill="lightgreen", alpha = 0.25) + stat_cor(aes(label = paste(..p.label..)), label.x = 0, label.y = 2.5) + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~")), label.x = 0, label.y = 3) + geom_text(label=statstoPlot$external_gene_name, colour = "darkgreen") +
  ggtitle("SNORD116-III") + theme(plot.title = element_text(hjust = 0.5)) + labs(x = "Mean Counts", y = "abs(log2foldchange)")
plot_115_mean <- ggplot(statstoPlot, aes(x=mean115, y=abslog2FC_H9)) +
  geom_point(color = "gray") + geom_smooth(method = lm, colour="gray", fill="lightgray", alpha = 0.4) + stat_cor(aes(label = paste(..p.label..)), label.x = 0, label.y = 2.5) + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~")), label.x = 0, label.y = 3) + geom_text(label=statstoPlot$external_gene_name, colour = "darkgray") +
  ggtitle("SNORD115") + theme(plot.title = element_text(hjust = 0.5)) + labs(x = "Mean Counts", y = "abs(log2foldchange)")
plot_group1_med <- ggplot(statstoPlot, aes(x=median116_I, y=abslog2FC_H9)) +
  geom_point(color = "blue") + geom_smooth(method = lm, colour="blue", fill = "lightblue", alpha = 0.25) + stat_cor(aes(label = paste(..p.label..)), label.x = 0, label.y = 2.5) + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~")), label.x = 0, label.y = 3) + geom_text(label=statstoPlot$external_gene_name, colour = "darkblue") +
  labs(x = "Median Counts", y = "abs(log2foldchange)")
plot_group2_med <- ggplot(statstoPlot, aes(x=median116_II, y=abslog2FC_H9)) +
  geom_point(color = "red") + geom_smooth(method = lm, colour="red", fill="pink", alpha = 0.35) + stat_cor(aes(label = paste(..p.label..)), label.x = 0, label.y = 2.5) + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~")), label.x = 0, label.y = 3) + geom_text(label=statstoPlot$external_gene_name, colour = "darkred") +
  labs(x = "Median Counts", y = "abs(log2foldchange)")
plot_group3_med <- ggplot(statstoPlot, aes(x=median116_III, y=abslog2FC_H9)) +
  geom_point(color = "green") + geom_smooth(method = lm, colour="green", fill="lightgreen", alpha = 0.25) + stat_cor(aes(label = paste(..p.label..)), label.x = 0, label.y = 2.5) + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~")), label.x = 0, label.y = 3) + geom_text(label=statstoPlot$external_gene_name, colour = "darkgreen") +
  labs(x = "Median Counts", y = "abs(log2foldchange)")
plot_115_med <- ggplot(statstoPlot, aes(x=median115, y=abslog2FC_H9)) +
  geom_point(color = "gray") + geom_smooth(method = lm, colour="gray", fill="lightgray", alpha = 0.4) + stat_cor(aes(label = paste(..p.label..)), label.x = 0, label.y = 2.5) + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~")), label.x = 0, label.y = 3) + geom_text(label=statstoPlot$external_gene_name, colour = "darkgray") +
  labs(x = "Median Counts", y = "abs(log2foldchange)")
plot_group1_sum <- ggplot(statstoPlot, aes(x=sum116_I, y=abslog2FC_H9)) +
  geom_point(color = "blue") + geom_smooth(method = lm, colour="blue", fill = "lightblue", alpha = 0.25) + stat_cor(aes(label = paste(..p.label..)), label.x = 0, label.y = 2.5) + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~")), label.x = 0, label.y = 3) + geom_text(label=statstoPlot$external_gene_name, colour = "darkblue") +
  labs(x = "Sum Counts", y = "abs(log2foldchange)")
plot_group2_sum <- ggplot(statstoPlot, aes(x=sum116_II, y=abslog2FC_H9)) +
  geom_point(color = "red") + geom_smooth(method = lm, colour="red", fill="pink", alpha = 0.35) + stat_cor(aes(label = paste(..p.label..)), label.x = 0, label.y = 2.5) + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~")), label.x = 0, label.y = 3) + geom_text(label=statstoPlot$external_gene_name, colour = "darkred") +
  labs(x = "Sum Counts", y = "abs(log2foldchange)")
plot_group3_sum <- ggplot(statstoPlot, aes(x=sum116_III, y=abslog2FC_H9)) +
  geom_point(color = "green") + geom_smooth(method = lm, colour="green", fill="lightgreen", alpha = 0.25) + stat_cor(aes(label = paste(..p.label..)), label.x = 0, label.y = 2.5) + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~")), label.x = 0, label.y = 3) + geom_text(label=statstoPlot$external_gene_name, colour = "darkgreen") +
  labs(x = "Sum Counts", y = "abs(log2foldchange)")
plot_115_sum <- ggplot(statstoPlot, aes(x=sum115, y=abslog2FC_H9)) +
  geom_point(color = "gray") + geom_smooth(method = lm, colour="gray", fill="lightgray", alpha = 0.4) + stat_cor(aes(label = paste(..p.label..)), label.x = 0, label.y = 2.5) + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~")), label.x = 0, label.y = 3) + geom_text(label=statstoPlot$external_gene_name, colour = "darkgray") +
  labs(x = "Sum Counts", y = "abs(log2foldchange)")
#To arrange plots
allplot <- ggarrange(plot_group1_mean, plot_group2_mean, plot_group3_mean, plot_115_mean, plot_group1_med, plot_group2_med, plot_group3_med, plot_115_med, plot_group1_sum, plot_group2_sum, plot_group3_sum, plot_115_sum,
          ncol = 4, nrow = 3)
allplot
#Save plot
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/scatterPlots_targetingeventsvsfoldchange.pdf", width = 10, height = 8)
allplot
dev.off()
  #No, there is no significant correlation.

##Create meta-gene plot to analyze distribution of SNORD116 binding events across a gene body
#Import HOMER output
meta116 <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/metaGeneProfile_SNORD116.txt", sep="\t", head=T))
meta115 <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/metaGeneProfile_SNORD115.txt", sep="\t", head=T))
#Rename columns
colnames(meta116) <- c("MetaGene_Profile", "Peaks_per_bp_per_gene")
colnames(meta115) <- c("MetaGene_Profile", "Peaks_per_bp_per_gene")
#Plot
plot116meta<- ggplot(meta116, aes(x=MetaGene_Profile, y=Peaks_per_bp_per_gene)) +
  geom_line() + ggtitle("Meta-Gene Profile", subtitle = "Distribution of Predicted SNORD116 Interactions") +
  scale_x_continuous(name = "Meta-Gene Profile", breaks = seq(from = -5000, to = 25000, by = 5000)) + ylab("Peaks per bp per gene") +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plot116meta
plot115meta<- ggplot(meta115, aes(x=MetaGene_Profile, y=Peaks_per_bp_per_gene)) +
  geom_line() + ggtitle("Meta-Gene Profile", subtitle = "Distribution of Predicted SNORD115 Interactions") +
  scale_x_continuous(name = "Meta-Gene Profile", breaks = seq(from = -5000, to = 25000, by = 5000)) + ylab("Peaks per bp per gene") +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plot115meta
#Make overlaid plot
  #Combine data.frames
metageneComb <- merge(meta116, meta115, by="MetaGene_Profile")
  #Rename columns
colnames(metageneComb) <- c("MetaGene_Profile", "SNORD116", "SNORD115")
  #Plot
plotCombined<- ggplot(metageneComb, aes(x=MetaGene_Profile)) +
  geom_line(aes(y = SNORD116), color = "black") +  geom_line(aes(y = SNORD115), color = "darkgray") +
  ggtitle("Meta-Gene Profile", subtitle = "Distribution of Predicted Interactions") +
  scale_x_continuous(name = "Meta-Gene Profile", breaks = seq(from = -5000, to = 25000, by = 5000)) + ylab("Peaks per bp per gene") +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plotCombined

##Create alternative meta-gene plot to analyze distribution of SNORD116 binding events across a gene body
#Import results from kallisto
transcriptsALL <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/transcript_TPMs_allWT_gene.tsv", sep = "\t", header = TRUE)
#Change column names
colnames(transcriptsALL) <- c("transcript_id","gene_id","havana_gene","havana_transcript","transcript_name","gene_name","length","transcript_type","length","CT2-ngn2-s1_counts","CT2-ngn2-s2_counts","CT2-ngn2-s3_counts","CT2-ngn2-s4_counts","CT2-ngn2-s5_counts","CT2-ngn2-s6_counts","H9-NGN-1_counts","H9-NGN-2_counts","H9-NGN-3_counts","H9-NGN-4_counts","H9-NGN-5_counts","H9-NGN-6_counts")
#Create column for summing all counts
transcriptsALL$sum_counts <- rowSums(transcriptsALL[10:21])
#Create a column for calculating the average counts
transcriptsALL$mean_counts <- rowMeans(transcriptsALL[10:21])
#Read in distance measure file for SNORD116 targets
SNORD116all.dist <- read.delim ("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/chr15SNORD116targets.dist.measures.txt", header = T)
SNORD116Iall.dist<- read.delim ("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116groupItargets.dist.measures.txt", header = T)
SNORD116IIall.dist<- read.delim ("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116groupIItargets.dist.measures.txt", header = T)
SNORD116IIIall.dist<- read.delim ("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116groupIIItargets.dist.measures.txt", header = T)
SNORD115all.dist<- read.delim ("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/chr15SNORD115targets.dist.measures.txt", header = T)
#Pull out transcripts detected in metaPlotR analysis
temp_transcripts116 <- transcriptsALL[(transcriptsALL$transcript_id) %in% SNORD116all.dist$refseqID, ]
temp_transcripts115 <- transcriptsALL[(transcriptsALL$transcript_id) %in% SNORD115all.dist$refseqID, ]
#Sort by gene_id and then highest counts
temp_transcripts116 <- temp_transcripts116[order(temp_transcripts116$gene_id, -temp_transcripts116$mean_counts),]
temp_transcripts115 <- temp_transcripts115[order(temp_transcripts115$gene_id, -temp_transcripts115$mean_counts),]
#Filter by removing duplicate gene_id's
temp_transcripts116HI <- temp_transcripts116[!duplicated(temp_transcripts116$gene_id),]
temp_transcripts115HI <- temp_transcripts115[!duplicated(temp_transcripts115$gene_id),]
#Filter transcripts from distance measure file that most expressed in our data
SNORD116exp.dist <- SNORD116all.dist[(SNORD116all.dist$refseqID) %in% temp_transcripts116HI$transcript_id, ]
SNORD116Iexp.dist <- SNORD116Iall.dist[(SNORD116Iall.dist$refseqID) %in% temp_transcripts116HI$transcript_id, ]
SNORD116IIexp.dist <- SNORD116IIall.dist[(SNORD116IIall.dist$refseqID) %in% temp_transcripts116HI$transcript_id, ]
SNORD116IIIexp.dist <- SNORD116IIIall.dist[(SNORD116IIIall.dist$refseqID) %in% temp_transcripts116HI$transcript_id, ]
SNORD115exp.dist <- SNORD115all.dist[(SNORD115all.dist$refseqID) %in% temp_transcripts115HI$transcript_id, ]
#Re-scale the widths of the 5'UTR and 3'UTR relative to the CDS (which is set constant to a width of 1 unit)
utr5.SF <- median(SNORD116exp.dist$utr5_size, na.rm = T)/median(SNORD116exp.dist$cds_size, na.rm = T)
utr3.SF <- median(SNORD116exp.dist$utr3_size, na.rm = T)/median(SNORD116exp.dist$cds_size, na.rm = T)
utr5gpI.SF <- median(SNORD116Iexp.dist$utr5_size, na.rm = T)/median(SNORD116Iexp.dist$cds_size, na.rm = T)
utr3gpI.SF <- median(SNORD116Iexp.dist$utr3_size, na.rm = T)/median(SNORD116Iexp.dist$cds_size, na.rm = T)
utr5gpII.SF <- median(SNORD116IIexp.dist$utr5_size, na.rm = T)/median(SNORD116IIexp.dist$cds_size, na.rm = T)
utr3gpII.SF <- median(SNORD116IIexp.dist$utr3_size, na.rm = T)/median(SNORD116IIexp.dist$cds_size, na.rm = T)
utr5gpIII.SF <- median(SNORD116IIIexp.dist$utr5_size, na.rm = T)/median(SNORD116IIIexp.dist$cds_size, na.rm = T)
utr3gpIII.SF <- median(SNORD116IIIexp.dist$utr3_size, na.rm = T)/median(SNORD116IIIexp.dist$cds_size, na.rm = T)
utr5_115.SF <- median(SNORD115exp.dist$utr5_size, na.rm = T)/median(SNORD115exp.dist$cds_size, na.rm = T)
utr3_115.SF <- median(SNORD115exp.dist$utr3_size, na.rm = T)/median(SNORD115exp.dist$cds_size, na.rm = T)
#Assign the regions to new dataframes
utr5.SNORD116.dist <- SNORD116exp.dist[SNORD116exp.dist$rel_location < 1, ]
cds.SNORD116.dist <- SNORD116exp.dist [SNORD116exp.dist$rel_location < 2 & SNORD116exp.dist$rel_location >= 1, ]
utr3.SNORD116.dist <- SNORD116exp.dist[SNORD116exp.dist$rel_location >= 2, ]
utr5.SNORD116I.dist <- SNORD116Iexp.dist[SNORD116Iexp.dist$rel_location < 1, ]
cds.SNORD116I.dist <- SNORD116Iexp.dist [SNORD116Iexp.dist$rel_location < 2 & SNORD116Iexp.dist$rel_location >= 1, ]
utr3.SNORD116I.dist <- SNORD116Iexp.dist[SNORD116Iexp.dist$rel_location >= 2, ]
utr5.SNORD116II.dist <- SNORD116IIexp.dist[SNORD116IIexp.dist$rel_location < 1, ]
cds.SNORD116II.dist <- SNORD116IIexp.dist [SNORD116IIexp.dist$rel_location < 2 & SNORD116IIexp.dist$rel_location >= 1, ]
utr3.SNORD116II.dist <- SNORD116IIexp.dist[SNORD116IIexp.dist$rel_location >= 2, ]
utr5.SNORD116III.dist <- SNORD116IIIexp.dist[SNORD116IIIexp.dist$rel_location < 1, ]
cds.SNORD116III.dist <- SNORD116IIIexp.dist [SNORD116IIIexp.dist$rel_location < 2 & SNORD116IIIexp.dist$rel_location >= 1, ]
utr3.SNORD116III.dist <- SNORD116IIIexp.dist[SNORD116IIIexp.dist$rel_location >= 2, ]
utr5.SNORD115.dist <- SNORD115exp.dist[SNORD115exp.dist$rel_location < 1, ]
cds.SNORD115.dist <- SNORD115exp.dist [SNORD115exp.dist$rel_location < 2 & SNORD115exp.dist$rel_location >= 1, ]
utr3.SNORD115.dist <- SNORD115exp.dist[SNORD115exp.dist$rel_location >= 2, ]
#Rescale 5'UTR and 3'UTR
library("scales")
utr5.SNORD116.dist$rel_location <- rescale(utr5.SNORD116.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0,1))
utr3.SNORD116.dist$rel_location <- rescale(utr3.SNORD116.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2,3))
utr5.SNORD116I.dist$rel_location <- rescale(utr5.SNORD116I.dist$rel_location, to = c(1-utr5gpI.SF, 1), from = c(0,1))
utr3.SNORD116I.dist$rel_location <- rescale(utr3.SNORD116I.dist$rel_location, to = c(2, 2+utr3gpI.SF), from = c(2,3))
utr5.SNORD116II.dist$rel_location <- rescale(utr5.SNORD116II.dist$rel_location, to = c(1-utr5gpII.SF, 1), from = c(0,1))
utr3.SNORD116II.dist$rel_location <- rescale(utr3.SNORD116II.dist$rel_location, to = c(2, 2+utr3gpII.SF), from = c(2,3))
utr5.SNORD116III.dist$rel_location <- rescale(utr5.SNORD116III.dist$rel_location, to = c(1-utr5gpIII.SF, 1), from = c(0,1))
utr3.SNORD116III.dist$rel_location <- rescale(utr3.SNORD116III.dist$rel_location, to = c(2, 2+utr3gpIII.SF), from = c(2,3))
utr5.SNORD115.dist$rel_location <- rescale(utr5.SNORD115.dist$rel_location, to = c(1-utr5_115.SF, 1), from = c(0,1))
utr3.SNORD115.dist$rel_location <- rescale(utr3.SNORD115.dist$rel_location, to = c(2, 2+utr3_115.SF), from = c(2,3))
#Combine
SNORD116.metagene.coord <- c(utr5.SNORD116.dist$rel_location, cds.SNORD116.dist$rel_location, utr3.SNORD116.dist$rel_location)
SNORD116I.metagene.coord <- c(utr5.SNORD116I.dist$rel_location, cds.SNORD116I.dist$rel_location, utr3.SNORD116I.dist$rel_location)
SNORD116II.metagene.coord <- c(utr5.SNORD116II.dist$rel_location, cds.SNORD116II.dist$rel_location, utr3.SNORD116II.dist$rel_location)
SNORD116III.metagene.coord <- c(utr5.SNORD116III.dist$rel_location, cds.SNORD116III.dist$rel_location, utr3.SNORD116III.dist$rel_location)
SNORD115.metagene.coord <- c(utr5.SNORD115.dist$rel_location, cds.SNORD115.dist$rel_location, utr3.SNORD115.dist$rel_location)
#Draw combined plots
metagene.cord <- c(SNORD116.metagene.coord, SNORD116I.metagene.coord, SNORD116II.metagene.coord, SNORD116III.metagene.coord)
groups116 <- c(rep("SNORD116", length(SNORD116.metagene.coord)), 
         rep("SNORD116-I", length(SNORD116I.metagene.coord)),
         rep("SNORD116-II", length(SNORD116II.metagene.coord)),
         rep("SNORD116-III", length(SNORD116III.metagene.coord))) 
comb <- data.frame(metagene.cord, groups116)
cols <- c("SNORD116" = "black", "SNORD116-I" = "#AA4499", "SNORD116-II" = "#CC6677", "SNORD116-III" = "#882255")
combPlot <- ggplot(comb) + geom_density(aes(x = metagene.cord, colour = groups116)) + theme_bw() + 
  scale_colour_manual(values = cols) +
  geom_vline(xintercept = 1:2, col = "#DDCC77")
combPlot
#Save plot
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/metaGenePlotSNORD116groups.pdf", width = 10, height = 8)
combPlot
dev.off()
#Single plot
meta115df <- as.data.frame(SNORD115.metagene.coord)
pSNORD115 <- ggplot(meta115df) + geom_density(aes(x = SNORD115.metagene.coord)) + theme_bw() + 
  scale_x_continuous(limits = c(0.5, 3.5), breaks = c(1, 2, 3)) +
  geom_vline(xintercept = 1:2, col = "#DDCC77")
pSNORD115
#Save plot
pdf("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/metaGenePlotSNORD115.pdf", width = 10, height = 8)
pSNORD115
dev.off()

#SCRATCH WORK#
##Run Mann-Whitney U test to test significance
#Create empty variables
test.stat.group1 <- list()
test.stat.group2 <- list()
test.stat.group3 <- list()
results_group1 <- list()
results_group2 <- list()
results_group3 <- list()
#Run for loop
for (i in enriched_list)
{
  #Set variables
  y <- snord115[[i]]$count
  x1 <- snord116[[i]][snord116[[i]]$orig_group == "116-I",5]
  results_group1[[i]] <- wilcox.test(x1,y, exact = FALSE)
  test.stat.group1[[i]] <- results_group1[[i]]$p.value
  x2 <- snord116[[i]][snord116[[i]]$orig_group == "116-II",5]
  results_group2[[i]] <- wilcox.test(x2,y, exact = FALSE)
  test.stat.group2[[i]] <- results_group2[[i]]$p.value
  x3 <- snord116[[i]][snord116[[i]]$orig_group == "116-III",5]
  results_group3[[i]] <- wilcox.test(x3,y, exact = FALSE)
  test.stat.group3[[i]] <- results_group3[[i]]$p.value
}
#Create list of transcriptome iNeurons for analysis
genes.table_all <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name", "chromosome_name"), values= row.names(allnonZeroStats), mart= mart)
allnonZeroStats_noMT <- allnonZeroStats[which(genes.table_all$chromosome_name != "MT"),]
write.table(row.names(allnonZeroStats_noMT), file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/allnonZeroStats_noMT_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Create list of binding events for parsing files
targets116chr15 <- targets116 %>% filter(!targets116$sno_id %in% c("ENSG00000202498_SNORD116(-)", "ENSG00000212553_SNORD116(-)", "ENSG00000252985_SNORD116(+)"))
targets116_list <- unique(targets116chr15$target_id)
#To save
write.table(targets116_list, file = "../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116_targets_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#To load
targets116_list <- as.data.frame(read.table("../Cotney_Lab/PWS_RNASeq/STAR/snoGloBe/SNORD116_targets_list.txt", sep="\t", head=F))