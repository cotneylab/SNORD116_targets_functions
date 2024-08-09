#! /usr/bin/bash
#SBATCH --job-name=Venn_PWSexp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mem=64G
#SBATCH --mail-user=rgilmore@uchc.edu
#SBATCH -o sbatch_%x_%j.out
#SBATCH -e sbatch_%x_%j.err

##DESCRIPTION: This is script is a loop which runs tests of overlap to text significance of intersects found in ggVennDiagram.
##USAGE: sbatch venn_intersect_test_rbg.sbatch
##INSTRUCTION: User must upload exported gene lists from R to cluster & modify names of files.

hostname; date

#Go to directory for analysis
cd /home/FCAM/rgilmore/analysis/RNAseq/PWS_RNAseq/STAR/PWS_Exp/VennAnalysis

for ((n=1; n<=10000; n++))
do
	#Shuffle lists
	shuf AllH9Genes.txt > tempH9.txt
	shuf AllCT2Genes.txt > tempCT2.txt
	#Create gene lists using same number of genes as in analysis
	head -n 1578 tempH9.txt > A.txt
	head -n 2671 tempCT2.txt > B.txt
	tail -n 2317 tempCT2.txt > C.txt
	tail -n 1422 tempH9.txt > D.txt
	#Sort files
	sort -k1,1 A.txt > A_sort.txt
	sort -k1,1 B.txt > B_sort.txt
	sort -k1,1 C.txt > C_sort.txt
	sort -k1,1 D.txt > D_sort.txt
	#Join files and count for overlaps
	join A_sort.txt B_sort.txt | wc -l >> A_B.txt
	join C_sort.txt D_sort.txt | wc -l >> C_D.txt
done
