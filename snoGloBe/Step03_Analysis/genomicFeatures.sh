#Go into working directory
cd /home/FCAM/rgilmore/analysis/snoGloBe/genomicFeatures
#Load bedtools
module load bedtools/2.29.0
#Filter gtf file based on lists of genes
	#Background of iNeuron transcriptome (used to sample for random permutations)
cat /home/FCAM/jcotney/ANALYSIS/snoglobe/Homo_sapiens.GRCh38.88.gtf | grep -f allnonZeroStats_noMT_list.txt > Homo_sapiens.GRCh38.88_sampledGenes.gtf
	#Background of 42 shared genes
cat /home/FCAM/jcotney/ANALYSIS/snoglobe/Homo_sapiens.GRCh38.88.gtf | grep -f ../ENSEMBLsharedsmlgDel.txt > Homo_sapiens.GRCh38.88_sharedGenes.gtf
	#Background of random gene lists
cd /home/FCAM/rgilmore/analysis/snoGloBe/permutations
cat random_gene_list*.txt > all_randomGenes.txt
cat /home/FCAM/jcotney/ANALYSIS/snoglobe/Homo_sapiens.GRCh38.88.gtf | grep -f all_randomGenes.txt > Homo_sapiens.GRCh38.88_randomGenes.gtf

## For TRANSCRIPT files ##
cd /home/FCAM/rgilmore/analysis/snoGloBe/genomicFeatures
#iNeuron
	#Create sorted bed file
cat Homo_sapiens.GRCh38.88_sampledGenes.gtf | grep -v "transcript_support_level \"NA\"" | grep -v "transcript_support_level \"5\"" | grep -v "transcript_support_level \"4\"" | grep "tag \"basic\"" | awk 'OFS="\t" {if ($3=="transcript") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' | awk '{print "chr"$0}' | sed -e 's/chrKI/chrUn_KI/g' -e 's/chrGL/chrUn_GL/g'| grep -v chrUn | grep -v chrMT | sort -k1,1 -k2,2n > sampledTranscript.sorted.bed
    #Merge
bedtools merge -s -i sampledTranscript.sorted.bed | sort -k1,1 -k2,2n > sampledTranscript_merge.sorted.bed
	#To create file for uploading to UCSC
		#sampledTranscript_merge.sorted.bed is good as is
cat sampledTranscript.sorted.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"0"\t"$6}' > sampledTranscript_hg38_UCSC.bed
#sharedGenes
	#Create sorted bed file
cat Homo_sapiens.GRCh38.88_sharedGenes.gtf | grep -v "transcript_support_level \"NA\"" | grep -v "transcript_support_level \"5\"" | grep -v "transcript_support_level \"4\"" | grep "tag \"basic\"" | awk 'OFS="\t" {if ($3=="transcript") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' | awk '{print "chr"$0}' | sed -e 's/chrKI/chrUn_KI/g' -e 's/chrGL/chrUn_GL/g'| grep -v chrUn | grep -v chrMT | sort -k1,1 -k2,2n > sharedTranscript.sorted.bed
    #Merge
bedtools merge -s -i sharedTranscript.sorted.bed | sort -k1,1 -k2,2n > sharedTranscript_merge.sorted.bed
#randomGenes
cd /home/FCAM/rgilmore/analysis/snoGloBe/permutations
	#Create sorted bed file
cat Homo_sapiens.GRCh38.88_randomGenes.gtf | grep -v "transcript_support_level \"NA\"" | grep -v "transcript_support_level \"5\"" | grep -v "transcript_support_level \"4\"" | grep "tag \"basic\"" | awk 'OFS="\t" {if ($3=="transcript") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' | awk '{print "chr"$0}' | sed -e 's/chrKI/chrUn_KI/g' -e 's/chrGL/chrUn_GL/g'| grep -v chrUn | grep -v chrMT | sort -k1,1 -k2,2n > randomTranscript.sorted.bed
    #Merge
bedtools merge -s -i randomTranscript.sorted.bed | sort -k1,1 -k2,2n > randomTranscript_merge.sorted.bed

## For EXON files ##
cd /home/FCAM/rgilmore/analysis/snoGloBe/genomicFeatures
#iNeuron
	#Create sorted bed file
cat Homo_sapiens.GRCh38.88_sampledGenes.gtf | grep -v "transcript_support_level \"NA\"" | grep -v "transcript_support_level \"5\"" | grep -v "transcript_support_level \"4\"" | grep "tag \"basic\"" | awk 'OFS="\t" {if ($3=="exon") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' | awk '{print "chr"$0}' | sed -e 's/chrKI/chrUn_KI/g' -e 's/chrGL/chrUn_GL/g'| grep -v chrUn | grep -v chrMT | sort -k1,1 -k2,2n > sampledExon.sorted.bed
    #Merge
bedtools merge -s -i sampledExon.sorted.bed | sort -k1,1 -k2,2n > sampledExon_merge.sorted.bed
	#To create file for uploading to UCSC
		#sampledExon_merge.sorted.bed is good as is
cat sampledExon.sorted.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"0"\t"$6}' > sampledExon_hg38_UCSC.sorted.bed
#sharedGenes
	#Create sorted bed file
cat Homo_sapiens.GRCh38.88_sharedGenes.gtf | grep -v "transcript_support_level \"NA\"" | grep -v "transcript_support_level \"5\"" | grep -v "transcript_support_level \"4\"" | grep "tag \"basic\"" | awk 'OFS="\t" {if ($3=="exon") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' | awk '{print "chr"$0}' | sed -e 's/chrKI/chrUn_KI/g' -e 's/chrGL/chrUn_GL/g'| grep -v chrUn | grep -v chrMT | sort -k1,1 -k2,2n > sharedExon.sorted.bed
    #Merge
bedtools merge -s -i sharedExon.sorted.bed | sort -k1,1 -k2,2n > sharedExon_merge.sorted.bed
#randomGenes
cd /home/FCAM/rgilmore/analysis/snoGloBe/permutations
	#Create sorted bed file
cat Homo_sapiens.GRCh38.88_randomGenes.gtf | grep -v "transcript_support_level \"NA\"" | grep -v "transcript_support_level \"5\"" | grep -v "transcript_support_level \"4\"" | grep "tag \"basic\"" | awk 'OFS="\t" {if ($3=="exon") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' | awk '{print "chr"$0}' | sed -e 's/chrKI/chrUn_KI/g' -e 's/chrGL/chrUn_GL/g'| grep -v chrUn | grep -v chrMT | sort -k1,1 -k2,2n > randomExon.sorted.bed
    #Merge
bedtools merge -s -i randomExon.sorted.bed | sort -k1,1 -k2,2n > randomExon_merge.sorted.bed

## For INTRON file ##
cd /home/FCAM/rgilmore/analysis/snoGloBe/genomicFeatures
	#Sort chrom sizes file
sort -k1,1 /home/FCAM/jcotney/GENOME/hg38/dna/hg38.chrom.sizes > hg38_sorted.chrom.sizes
#iNeuron
	#First take complement of exon file
bedtools complement -L -i sampledExon.sorted.bed -g hg38_sorted.chrom.sizes | sort -k1,1 -k2,2n > complement_sampledExon.sorted.bed
    #Intersect
bedtools intersect -a complement_sampledExon.sorted.bed -b sampledTranscript_merge.sorted.bed | sort -k1,1 -k2,2n > sampledIntron_fromTranscript.sorted.bed
#sharedGenes
	#First take complement of exon file
bedtools complement -L -i sharedExon.sorted.bed -g hg38_sorted.chrom.sizes | sort -k1,1 -k2,2n > complement_sharedExon.sorted.bed
    #Intersect
bedtools intersect -a complement_sharedExon.sorted.bed -b sharedTranscript_merge.sorted.bed | sort -k1,1 -k2,2n > sharedIntron_fromTranscript.sorted.bed
#randomGenes
cd /home/FCAM/rgilmore/analysis/snoGloBe/permutations
	#First take complement of exon file
bedtools complement -L -i randomExon.sorted.bed -g ../genomicFeatures/hg38_sorted.chrom.sizes | sort -k1,1 -k2,2n > complement_randomExon.sorted.bed
    #Intersect
bedtools intersect -a complement_randomExon.sorted.bed -b randomTranscript_merge.sorted.bed | sort -k1,1 -k2,2n > randomIntron_fromTranscript.sorted.bed

## For JUNCTION file ##
cd /home/FCAM/rgilmore/analysis/snoGloBe/genomicFeatures
#iNeuron
	#Create a file for intron-exon junctions
cat sampledIntron_fromTranscript.sorted.bed | awk 'OFS="\t" { print $1,$2-1,$2+1"\n"$1,$3-1,$3+1}' | sort -k1,1 -k2,2n > sampledJunction.sorted.bed
#sharedGenes
	#Create a file for intron-exon junctions
cat sharedIntron_fromTranscript.sorted.bed | awk 'OFS="\t" { print $1,$2-1,$2+1"\n"$1,$3-1,$3+1}' | sort -k1,1 -k2,2n > sharedJunction.sorted.bed
#randomGenes
cd /home/FCAM/rgilmore/analysis/snoGloBe/permutations
	#Create a file for intron-exon junctions
cat randomIntron_fromTranscript.sorted.bed | awk 'OFS="\t" { print $1,$2-1,$2+1"\n"$1,$3-1,$3+1}' | sort -k1,1 -k2,2n > randomJunction.sorted.bed

## For all CDS ##
cd /home/FCAM/rgilmore/analysis/snoGloBe/genomicFeatures
#iNeuron
	#Create sorted bed file
cat Homo_sapiens.GRCh38.88_sampledGenes.gtf | grep -v "transcript_support_level \"NA\"" | grep -v "transcript_support_level \"5\"" | grep -v "transcript_support_level \"4\"" | grep "tag \"basic\"" | awk 'OFS="\t" {if ($3=="CDS") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' | awk '{print "chr"$0}' | sed -e 's/chrKI/chrUn_KI/g' -e 's/chrGL/chrUn_GL/g'| grep -v chrUn | grep -v chrMT | sort -k1,1 -k2,2n > sampledCDSall.sorted.bed
    #To create file for uploading to UCSC
cat sampledCDSall.sorted.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"0"\t"$6}' > sampledCDSall_hg38_UCSC.bed
#sharedGenes
	#Create sorted bed file
cat Homo_sapiens.GRCh38.88_sharedGenes.gtf | grep -v "transcript_support_level \"NA\"" | grep -v "transcript_support_level \"5\"" | grep -v "transcript_support_level \"4\"" | grep "tag \"basic\"" | awk 'OFS="\t" {if ($3=="CDS") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' | awk '{print "chr"$0}' | sed -e 's/chrKI/chrUn_KI/g' -e 's/chrGL/chrUn_GL/g'| grep -v chrUn | grep -v chrMT | sort -k1,1 -k2,2n > sharedCDSall.sorted.bed
#randomGenes
cd /home/FCAM/rgilmore/analysis/snoGloBe/permutations
	#Create sorted bed file
cat Homo_sapiens.GRCh38.88_randomGenes.gtf | grep -v "transcript_support_level \"NA\"" | grep -v "transcript_support_level \"5\"" | grep -v "transcript_support_level \"4\"" | grep "tag \"basic\"" | awk 'OFS="\t" {if ($3=="CDS") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' | awk '{print "chr"$0}' | sed -e 's/chrKI/chrUn_KI/g' -e 's/chrGL/chrUn_GL/g'| grep -v chrUn | grep -v chrMT | sort -k1,1 -k2,2n > randomCDSall.sorted.bed

## For all 5'UTR ##
cd /home/FCAM/rgilmore/analysis/snoGloBe/genomicFeatures
#iNeuron
	#Create sorted bed file
cat Homo_sapiens.GRCh38.88_sampledGenes.gtf | grep -v "transcript_support_level \"NA\"" | grep -v "transcript_support_level \"5\"" | grep -v "transcript_support_level \"4\"" | grep "tag \"basic\"" | awk 'OFS="\t" {if ($3=="five_prime_utr") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' | awk '{print "chr"$0}' | sed -e 's/chrKI/chrUn_KI/g' -e 's/chrGL/chrUn_GL/g'| grep -v chrUn | grep -v chrMT | sort -k1,1 -k2,2n > sampled5UTRall.sorted.bed
    #To create file for uploading to UCSC
cat sampled5UTRall.sorted.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"0"\t"$6}' > sampled5UTRall_hg38_UCSC.bed
#sharedGenes
	#Create sorted bed file
cat Homo_sapiens.GRCh38.88_sharedGenes.gtf | grep -v "transcript_support_level \"NA\"" | grep -v "transcript_support_level \"5\"" | grep -v "transcript_support_level \"4\"" | grep "tag \"basic\"" | awk 'OFS="\t" {if ($3=="five_prime_utr") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' | awk '{print "chr"$0}' | sed -e 's/chrKI/chrUn_KI/g' -e 's/chrGL/chrUn_GL/g'| grep -v chrUn | grep -v chrMT | sort -k1,1 -k2,2n > shared5UTRall.sorted.bed
#randomGenes
cd /home/FCAM/rgilmore/analysis/snoGloBe/permutations
	#Create sorted bed file
cat Homo_sapiens.GRCh38.88_randomGenes.gtf | grep -v "transcript_support_level \"NA\"" | grep -v "transcript_support_level \"5\"" | grep -v "transcript_support_level \"4\"" | grep "tag \"basic\"" | awk 'OFS="\t" {if ($3=="five_prime_utr") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' | awk '{print "chr"$0}' | sed -e 's/chrKI/chrUn_KI/g' -e 's/chrGL/chrUn_GL/g'| grep -v chrUn | grep -v chrMT | sort -k1,1 -k2,2n > random5UTRall.sorted.bed

## For all 3'UTR ##
cd /home/FCAM/rgilmore/analysis/snoGloBe/genomicFeatures
#iNeuron
	#Create sorted bed file
cat Homo_sapiens.GRCh38.88_sampledGenes.gtf | grep -v "transcript_support_level \"NA\"" | grep -v "transcript_support_level \"5\"" | grep -v "transcript_support_level \"4\"" | grep "tag \"basic\"" | awk 'OFS="\t" {if ($3=="three_prime_utr") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' | awk '{print "chr"$0}' | sed -e 's/chrKI/chrUn_KI/g' -e 's/chrGL/chrUn_GL/g'| grep -v chrUn | grep -v chrMT | sort -k1,1 -k2,2n > sampled3UTRall.sorted.bed
    #To create file for uploading to UCSC
cat sampled3UTRall.sorted.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"0"\t"$6}' > sampled3UTRall_hg38_UCSC.bed
#sharedGenes
	#Create sorted bed file
cat Homo_sapiens.GRCh38.88_sharedGenes.gtf | grep -v "transcript_support_level \"NA\"" | grep -v "transcript_support_level \"5\"" | grep -v "transcript_support_level \"4\"" | grep "tag \"basic\"" | awk 'OFS="\t" {if ($3=="three_prime_utr") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' | awk '{print "chr"$0}' | sed -e 's/chrKI/chrUn_KI/g' -e 's/chrGL/chrUn_GL/g'| grep -v chrUn | grep -v chrMT | sort -k1,1 -k2,2n > shared3UTRall.sorted.bed
#randomGenes
cd /home/FCAM/rgilmore/analysis/snoGloBe/permutations
	#Create sorted bed file
cat Homo_sapiens.GRCh38.88_randomGenes.gtf | grep -v "transcript_support_level \"NA\"" | grep -v "transcript_support_level \"5\"" | grep -v "transcript_support_level \"4\"" | grep "tag \"basic\"" | awk 'OFS="\t" {if ($3=="three_prime_utr") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' | awk '{print "chr"$0}' | sed -e 's/chrKI/chrUn_KI/g' -e 's/chrGL/chrUn_GL/g'| grep -v chrUn | grep -v chrMT | sort -k1,1 -k2,2n > random3UTRall.sorted.bed

## Create bed files for SNORD116 targets vs random permutations ##
cd /home/FCAM/rgilmore/analysis/snoGloBe/permutations
cat sample_list.txt | while read i
do
 cat ${i}_snord116chr15 | sed '1d' | awk 'OFS="\t" {split($7, a, "_"); print $1,$2,$3,a[2],$5*100,$6}' | awk 'OFS="\t" {gsub(/[(+)]/, "", $4); print $0}' > ${i}_snord116targets_hg38.bed
done

## Parse SNORD116 bed file per copy ##
cd /home/FCAM/rgilmore/analysis/snoGloBe/genomicFeatures
	#For each SNORD116 copy vs shared genes
for i in {1..30}
do
	cat ../snord116targets_snoglobe_hg38.bed | awk -v i="$i" 'OFS="\t" {if ($4=="SNORD116-"i) {print $0}}' > snord116-${i}_hg38.bed 
done
	#For each SNORD116 copy vs random permutations
cd /home/FCAM/rgilmore/analysis/snoGloBe/permutations
cat sample_list.txt | while read i
do
	for n in {1..30}
	do
		cat ${i}_snord116targets_hg38.bed | awk -v n="$n" 'OFS="\t" {if ($4=="SNORD116-"n) {print $0}}' > snord116-${n}_${i}_hg38.bed
	done
done

## Parse SNORD115 bed file for chr15 only (& remove header) ##
cat ../snord115targets_snoglobe_hg38.bed | awk 'OFS="\t" {if ($4 ~ /SNORD115-/) {print $0}}' > snord115_allchr15_hg38.bed
	
## Run multiinter ##
cd /home/FCAM/rgilmore/analysis/snoGloBe/genomicFeatures
	#For each SNORD116 copy vs shared genes
for i in {1..30}
do
	bedtools multiinter -header -g hg38_sorted.chrom.sizes -names SNORD116-${i} exon intron junction CDS 5UTR 3UTR -i snord116-${i}_hg38.bed sampledExon.sorted.bed sampledIntron_fromTranscript.sorted.bed sampledJunction.sorted.bed sampledCDSall.sorted.bed sampled5UTRall.sorted.bed sampled3UTRall.sorted.bed > snord116-${i}_multiinter.txt
done
	#For background set
		#iNeuron
	bedtools multiinter -header -g hg38_sorted.chrom.sizes -names exon intron junction CDS 5UTR 3UTR -i sampledExon.sorted.bed sampledIntron_fromTranscript.sorted.bed sampledJunction.sorted.bed sampledCDSall.sorted.bed sampled5UTRall.sorted.bed sampled3UTRall.sorted.bed > background_multiinter.txt
		#sharedGenes
	bedtools multiinter -header -g hg38_sorted.chrom.sizes -names exon intron junction CDS 5UTR 3UTR -i sharedExon.sorted.bed sharedIntron_fromTranscript.sorted.bed sharedJunction.sorted.bed sharedCDSall.sorted.bed shared5UTRall.sorted.bed shared3UTRall.sorted.bed > backgroundShared_multiinter.txt
		#randomGenes
	cd /home/FCAM/rgilmore/analysis/snoGloBe/permutations
	bedtools multiinter -header -g ../genomicFeatures/hg38_sorted.chrom.sizes -names exon intron junction CDS 5UTR 3UTR -i randomExon.sorted.bed randomIntron_fromTranscript.sorted.bed randomJunction.sorted.bed randomCDSall.sorted.bed random5UTRall.sorted.bed random3UTRall.sorted.bed > backgroundRandom_multiinter.txt
	#For each SNORD116 copy vs random permutations (takes awhile - run as sbatch)
cat ../permutations/sample_list.txt | while read i
do
	for n in {1..30}
	do
		bedtools multiinter -header -g hg38_sorted.chrom.sizes -names SNORD116-${n} exon intron junction CDS 5UTR 3UTR -i ../permutations/snord116-${n}_${i}_hg38.bed sampledExon.sorted.bed sampledIntron_fromTranscript.sorted.bed sampledJunction.sorted.bed sampledCDSall.sorted.bed sampled5UTRall.sorted.bed sampled3UTRall.sorted.bed > ../permutations/snord116-${n}_${i}_multiinter.txt
	done
done
	#For all SNORD115 copies vs shared genes
bedtools multiinter -header -g hg38_sorted.chrom.sizes -names SNORD115 exon intron junction CDS 5UTR 3UTR -i snord115_allchr15_hg38.bed sampledExon.sorted.bed sampledIntron_fromTranscript.sorted.bed sampledJunction.sorted.bed sampledCDSall.sorted.bed sampled5UTRall.sorted.bed sampled3UTRall.sorted.bed > snord115_multiinter.txt
	#For all other chr15 SNORD copies vs shared genes
bedtools multiinter -header -g hg38_sorted.chrom.sizes -names otherchr15SNORDs exon intron junction CDS 5UTR 3UTR -i otherchr15SNORDstargets_snoglobe_hg38.bed sampledExon.sorted.bed sampledIntron_fromTranscript.sorted.bed sampledJunction.sorted.bed sampledCDSall.sorted.bed sampled5UTRall.sorted.bed sampled3UTRall.sorted.bed > otherchr15SNORDs_multiinter.txt

## Filter output tables ##
cd /home/FCAM/rgilmore/analysis/snoGloBe/genomicFeatures
	#For each SNORD116 copy vs shared genes
for i in {1..30}
do
	cat snord116-${i}_multiinter.txt |  awk 'OFS="\t" {if ($4>1 && $5 ~ /SNORD116/) {print $0,$3-$2}}' > snord116-${i}_multiinter_size.txt
done
	#For background set
		#iNeuron
	grep -v chrom background_multiinter.txt |  awk 'OFS="\t" {print $0,$3-$2}' > background_multiinter_size.txt
		#sharedGenes
	grep -v chrom backgroundShared_multiinter.txt |  awk 'OFS="\t" {print $0,$3-$2}' > backgroundShared_multiinter_size.txt
		#randomGenes
	cd /home/FCAM/rgilmore/analysis/snoGloBe/permutations
	grep -v chrom backgroundRandom_multiinter.txt |  awk 'OFS="\t" {print $0,$3-$2}' > backgroundRandom_multiinter_size.txt
	#For each SNORD116 copy vs random permutations (run as sbatch)
cd /home/FCAM/rgilmore/analysis/snoGloBe/permutations
cat sample_list.txt | while read i
do
	for n in {1..30}
	do
		cat snord116-${n}_${i}_multiinter.txt |  awk 'OFS="\t" {if ($4>1 && $5 ~ /SNORD116/) {print $0,$3-$2}}' > snord116-${n}_${i}_multiinter_size.txt
	done
done
	#For all SNORD115 copies vs shared genes
cd /home/FCAM/rgilmore/analysis/snoGloBe/genomicFeatures
cat snord115_multiinter.txt |  awk 'OFS="\t" {if ($4>1 && $5 ~ /SNORD115/) {print $0,$3-$2}}' > snord115_multiinter_size.txt
	#For all other chr15 SNORD copies vs shared genes
cat otherchr15SNORDs_multiinter.txt |  awk 'OFS="\t" {if ($4>1 && $5 ~ /SNORD/) {print $0,$3-$2}}' > otherchr15SNORDs_multiinter_size.txt

## Group output & combine into one table for reading into R ##
cd /home/FCAM/rgilmore/analysis/snoGloBe/genomicFeatures
	#For each SNORD116 copy vs shared genes
for i in {1..30}
do
	bedtools groupby -i snord116-${i}_multiinter_size.txt -g 5 -c 13 -o sum | sort -k1,1 | bedtools groupby -i stdin -g 1 -c 2 -o sum > snord116-${i}_grouped_table.txt
done
cat snord116-*_grouped_table.txt > allSNORD116_grouped_table.txt
	#For background set
		#iNeuron
	bedtools groupby -i background_multiinter_size.txt -g 5 -c 12 -o sum | sort -k1,1 | bedtools groupby -i stdin -g 1 -c 2 -o sum > background_grouped_table.txt
		#sharedGenes
	bedtools groupby -i backgroundShared_multiinter_size.txt -g 5 -c 12 -o sum | sort -k1,1 | bedtools groupby -i stdin -g 1 -c 2 -o sum > backgroundShared_grouped_table.txt
		#randomGenes
	cd /home/FCAM/rgilmore/analysis/snoGloBe/permutations
	bedtools groupby -i backgroundRandom_multiinter_size.txt -g 5 -c 12 -o sum | sort -k1,1 | bedtools groupby -i stdin -g 1 -c 2 -o sum > backgroundRandom_grouped_table.txt
	#For each SNORD116 copy vs random permutations (run as sbatch)
cd /home/FCAM/rgilmore/analysis/snoGloBe/permutations
cat sample_list.txt | while read i
do
	for n in {1..30}
	do
		bedtools groupby -i snord116-${n}_${i}_multiinter_size.txt -g 5 -c 13 -o sum | sort -k1,1 | bedtools groupby -i stdin -g 1 -c 2 -o sum > snord116-${n}_${i}_groupedTable.txt
	done
done
cat sample_list.txt | while read i
do
	for n in {1..30}
	do
		cat snord116-${n}_${i}_groupedTable.txt | awk -v i="$i" 'OFS="\t" {print $0, i}' > snord116-${n}_${i}_groupedTableSample.txt
	done
done
cat snord116-*_random_gene_list*_groupedTableSample.txt > allSNORD116_randomPerms_groupedTableSample.txt
	#For all SNORD115 copies vs shared genes
cd /home/FCAM/rgilmore/analysis/snoGloBe/genomicFeatures
bedtools groupby -i snord115_multiinter_size.txt -g 5 -c 13 -o sum | sort -k1,1 | bedtools groupby -i stdin -g 1 -c 2 -o sum > snord115_grouped_table.txt
	#For all other chr15 SNORD copies vs shared genes
bedtools groupby -i otherchr15SNORDs_multiinter_size.txt -g 5 -c 13 -o sum | sort -k1,1 | bedtools groupby -i stdin -g 1 -c 2 -o sum > otherchr15SNORDs_grouped_table.txt

## Move output table(s) to R ##