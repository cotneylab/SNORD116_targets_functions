##---Determining shared genes between lg & sm DEL lines---##
#Import files
lgDEL_up <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/lgDEL/UpGene_3analyses_lgDEL.txt", sep="\t", head=F)
smDEL_up <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/smDEL/UpGene_3analyses_smDEL.txt", sep="\t", head=F)
lgDEL_down <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/lgDEL/DownGene_3analyses_lgDEL.txt", sep="\t", head=F)
smDEL_down<- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/smDEL/DownGene_3analyses_smDEL.txt", sep="\t", head=F)
smDEL <- rbind(smDEL_up, smDEL_down)
lgDEL <- rbind(lgDEL_up, lgDEL_down)
#Merge to find shared genes
sharedlgsmDEL <- merge(lgDEL, smDEL, by="V1", sort=TRUE, all = FALSE)
sharedlgsmDEL_up <- merge(lgDEL_up, smDEL_up, by="V1", sort=TRUE, all = FALSE)
sharedlgsmDEL_down <- merge(lgDEL_down, smDEL_down, by="V1", sort=TRUE, all = FALSE)
  #Export list of shared genes
  write.table(sharedlgsmDEL, "../Cotney_Lab/PWS_RNASeq/STAR/ENSEMBLsharedsmlgDEL.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(sharedlgsmDEL_up, "../Cotney_Lab/PWS_RNASeq/STAR/ENSEMBLsharedsmlgDEL_up.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(sharedlgsmDEL_down, "../Cotney_Lab/PWS_RNASeq/STAR/ENSEMBLsharedsmlgDEL_down.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)

#Import files
lgDEL_name <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/lgDEL/allSigGeneName_3analyses_lgDEL.txt", sep="\t", head=F)
smDEL_name <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/smDEL/allSigGeneName_3analyses_smDEL.txt", sep="\t", head=F)
#Merge to find shared genes
SharedlgsmDEL_name <- merge(lgDEL_name, smDEL_name, by="V1", sort=TRUE, all = FALSE)
  #Export list of shared genes
  write.table(SharedlgsmDEL_name, "../Cotney_Lab/PWS_RNASeq/STAR/ENSEMBLsharedsmlgDEL.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)

##---Complex Upset---##
#Change first column name
colnames(lgDEL)[1] <- "Gene"
colnames(lgDEL_up)[1] <- "Gene"
colnames(lgDEL_down)[1] <- "Gene"
colnames(smDEL)[1] <- "Gene"
colnames(smDEL_up)[1] <- "Gene"
colnames(smDEL_down)[1] <- "Gene"
#Make column for TRUE/FALSE
lgDEL$lgDEL <- TRUE
smDEL$smDEL <- TRUE
lgDEL_up$lgDEL <- TRUE
smDEL_up$smDEL <- TRUE
lgDEL_down$lgDEL <- TRUE
smDEL_down$smDEL <- TRUE
#Merge lists
  #All dysregulated genes
all_upset_df <- merge(lgDEL, smDEL, by = 'Gene', all = TRUE)
  #Upregulated genes
up_upset_df <- merge(lgDEL_up, smDEL_up, by = 'Gene', all = TRUE)
  #Downregulated genes
down_upset_df <- merge(lgDEL_down, smDEL_down, by = 'Gene', all = TRUE)
#Set NA values to FALSE
all_upset_df[is.na(all_upset_df)] <- FALSE
down_upset_df[is.na(down_upset_df)] <- FALSE
up_upset_df[is.na(up_upset_df)] <- FALSE
#Plot
library("ggplot2")
library("ggupset")
library("ComplexUpset")
all_plot <- upset(all_upset_df, colnames(all_upset_df)[2:3], name = "DE Genes", width_ratio = 0.1, keep_empty_groups=TRUE,
                  matrix=(intersection_matrix(geom=geom_point(shape='square', size=3.5), segment=geom_segment(linetype='dotted'), outline_color=list(active='#FE6100', inactive='#929292'))
                          + scale_color_manual(values=c('TRUE'='#FFB000', 'FALSE'='#B5B4B4'), labels=c('TRUE'='yes', 'FALSE'='no'), breaks=c('TRUE', 'FALSE'), name='Intersection')
                          + scale_y_discrete(position='right')
                          + annotate(geom='text', label='Look here →', x='Comedy-Drama', y='Drama', size=5, hjust=1)),
                  queries=list(upset_query(intersect=c('smDEL', 'lgDEL'), color='#DC267F', fill='#DC267F', only_components=c('intersections_matrix', 'Intersection size')))) + 
  ggtitle('All Shared DE Genes (Disregarding Directionality)')
all_plot
up_plot <- upset(up_upset_df, colnames(up_upset_df)[2:3], name = "DE Genes", width_ratio = 0.1, keep_empty_groups=TRUE,
                 matrix=(intersection_matrix(geom=geom_point(shape='square', size=3.5), segment=geom_segment(linetype='dotted'), outline_color=list(active='#FE6100', inactive='#929292'))
                         + scale_color_manual(values=c('TRUE'='#FFB000', 'FALSE'='#B5B4B4'), labels=c('TRUE'='yes', 'FALSE'='no'), breaks=c('TRUE', 'FALSE'), name='Intersection')
                         + scale_y_discrete(position='right')
                         + annotate(geom='text', label='Look here →', x='Comedy-Drama', y='Drama', size=5, hjust=1)),
                 queries=list(upset_query(intersect=c('smDEL', 'lgDEL'), color='#DC267F', fill='#DC267F', only_components=c('intersections_matrix', 'Intersection size')))) + 
  ggtitle('Upregulated Shared DE Genes')
up_plot
down_plot <- upset(down_upset_df, colnames(down_upset_df)[2:3], name = "DE Genes", width_ratio = 0.1, keep_empty_groups=TRUE,
                   matrix=(intersection_matrix(geom=geom_point(shape='square', size=3.5), segment=geom_segment(linetype='dotted'), outline_color=list(active='#FE6100', inactive='#929292'))
                           + scale_color_manual(values=c('TRUE'='#FFB000', 'FALSE'='#B5B4B4'), labels=c('TRUE'='yes', 'FALSE'='no'), breaks=c('TRUE', 'FALSE'), name='Intersection')
                           + scale_y_discrete(position='right')
                           + annotate(geom='text', label='Look here →', x='Comedy-Drama', y='Drama', size=5, hjust=1)),
                   queries=list(upset_query(intersect=c('smDEL', 'lgDEL'), color='#DC267F', fill='#DC267F', only_components=c('intersections_matrix', 'Intersection size')))) + 
  ggtitle('Downregulated Shared DE Genes')
down_plot
#Save plots as pdf
pdf("../Cotney_Lab/PWS_RNASeq/STAR/lgVSsmDEL_Upset.pdf", width = 10, height = 8)
all_plot
up_plot
down_plot
dev.off()

##---Gene Ontology (on lg/sm DEL only---##
#Parse for lg/smDEL only
lgDEL_only <- all_upset_df[which(all_upset_df$lgDEL == 'TRUE' & all_upset_df$smDEL == 'FALSE'),]
smDEL_only <- all_upset_df[which(all_upset_df$lgDEL == 'FALSE' & all_upset_df$smDEL == 'TRUE'),]
#Make into lists
lgDEL_only <- lgDEL_only$Gene
smDEL_only <- smDEL_only$Gene
  #Save lists
  write.table(lgDEL_only, "../Cotney_Lab/PWS_RNASeq/STAR/lgDEL_only.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(smDEL_only, "../Cotney_Lab/PWS_RNASeq/STAR/smDEL_only.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
##Ran gene lists through Metascape
  ##Followed up with disgenet2r below!

##Following up on PWS ontology returned for lgDEL only dataset - figure out what genes are in group
  #See tutorial page: https://www.disgenet.org/static/disgenet2r/disgenet2r.html
library("devtools")
library("disgenet2r")
#Retrieve API key (must be registered for account first - register at https://www.disgenet.org/signup/)
disgenet_api_key <- get_disgenet_api_key(email = "rgilmore@uchc.edu", password = "")
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)
#Pull info for ontology term
PWS_disgenet <- disease2gene(disease = "C0032897", database = "ALL")
#Make results into data frame
PWS_disgenet_df <- as.data.frame(PWS_disgenet@qresult)
#Use biomart to find gene symbol for lgDEL only list
library("biomaRt")
mart <- useMart("ensembl", host = "https://apr2018.archive.ensembl.org/", dataset = "hsapiens_gene_ensembl")
genes.table <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene","hgnc_symbol"), values = lgDEL_only, mart = mart)
#Parse PWS ontology to find genes from lgDEL only contained within ontology
library("dplyr")
PWS_lgDEL_res <- PWS_disgenet_df %>% filter(PWS_disgenet_df$gene_symbol %in% genes.table$hgnc_symbol)
PWS_lgDEL_genelist <- PWS_lgDEL_res$gene_symbol
#Figure out where in the genome those genes are located
PWS_lgDEL_geneInfo <- getBM(filters = "hgnc_symbol", attributes= c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "description"), values = PWS_lgDEL_genelist, mart = mart)
  #Save table
  write.csv(PWS_lgDEL_geneInfo, file = "../Cotney_Lab/PWS_RNASeq/STAR/lgDELonly_PWS_Disgenet_terms.csv")
#Figure out where the other genes are (all in category)
PWS_C0032897_geneInfo <- getBM(filters = "hgnc_symbol", attributes= c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "description"), values = PWS_disgenet_df$gene_symbol, mart = mart)
#Get rid of weird chromosomes
PWS_C0032897_geneInfo <- PWS_C0032897_geneInfo %>% filter(!grepl('CHR', chromosome_name))
  #Save table
  write.csv(PWS_C0032897_geneInfo, file = "../Cotney_Lab/PWS_RNASeq/STAR/allPWS_Disgenet_terms.csv")

##---Gene Ontology on shared gene list---##
#Use biomart to find gene symbol
library("biomaRt")
mart <- useMart("ensembl", host = "https://apr2018.archive.ensembl.org/", dataset = "hsapiens_gene_ensembl")
genes.table_shared <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene","hgnc_symbol"), values = sharedlgsmDEL, mart = mart)
#Run disgenet2r
shared_DO <- disease_enrichment(entities = genes.table_shared$hgnc_symbol, vocabulary = "HGNC", database = "ALL")
#Create dataframe with results
summary_DO <- data.frame(shared_DO@qresult)
  #Save results table
  write.table(summary_DO, sep = "\t", file = "../Cotney_Lab/PWS_RNASeq/STAR/sharedGenes_DisGeNET_results.tsv", quote = FALSE)
#Plot
  #Calculate fold change to use on x-axis of dotplot
DOall_FC <- summary_DO
DOall_FC$Ratio <- sapply(DOall_FC$Ratio, function(x) eval(parse(text=x)))
DOall_FC$BgRatio <- sapply(DOall_FC$BgRatio, function(x) eval(parse(text=x)))
DOall_FC$FoldEnrichment <- DOall_FC$Ratio/DOall_FC$BgRatio
DOall_FC$log2FoldEnrichment <- log2(DOall_FC$FoldEnrichment)
DOall_FC$neglog10Padj <- -log(DOall_FC$FDR)
  #Order by p.adj
ordDOall_FC <- order(DOall_FC$neglog10Padj, decreasing=TRUE)
sortDOall_FC <- DOall_FC[ordDOall_FC, ]
  #Only use top 25 terms
sortDOall_FC <- sortDOall_FC[c(0:25),]
  #Order by logFoldEnrichment
sortDOall_FC <- sortDOall_FC[order(sortDOall_FC$log2FoldEnrichment),]
sortDOall_FC$Description <- factor(sortDOall_FC$Description, levels = sortDOall_FC$Description)
  #Make dotplots
library("ggplot2")
dotDOgg <- ggplot(data = sortDOall_FC, aes(x = log2FoldEnrichment, y = Description, 
                                              color = `neglog10Padj`, size = Count)) + 
  geom_point() +
  scale_colour_gradient2(low = "#E2E2E2", mid = "#37DC50", high = "#007B13", midpoint = 5.2) +
  theme_bw() + 
  labs(title="DisGeNET Enrichment for Shared Genes", y="", x="log2FoldEnrichment")
dotDOgg
  #Save plot as pdf
pdf("../Cotney_Lab/PWS_RNASeq/STAR/sharedGenes_DisGeNET_plot.pdf", width = 10, height = 8)
dotDOgg
dev.off()

##---To find list of TF's in shared genes---##
#Import file of all TFs: official list of human TFs from Lambert et al (PMID:29425488)
TFall <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/DatabaseExtract_v_1.01.txt", sep="\t", head=T)
#Delete first column
TFall <- TFall[-c(1)]
#Take only proteins thought to be TFs
TFall <- TFall[which(TFall$Is.TF. == "Yes"),]
#Make gene name rownames
rownames(TFall) <- TFall$Ensembl.ID
#Import shared gene list (if necessary)
SharedlgsmDEL_ID <- read.delim("../Cotney_Lab/PWS_RNASeq/STAR/ENSEMBLsharedsmlgDEL.txt", sep="\t", head=F)
rownames(SharedlgsmDEL_ID) <- SharedlgsmDEL_ID$V1
colnames(SharedlgsmDEL_ID) <- c("Ensembl.ID")
#Order list of all TFs & parse by shared gene list
TFall <- TFall[order(TFall$Ensembl.ID),]
sharedTF <- merge(TFall, SharedlgsmDEL_ID, by="Ensembl.ID", sort=TRUE, all = FALSE)
