#! /usr/bin/bash
#SBATCH --job-name=PermTest_human2mouse
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

##DESCRIPTION: This is script is a loop which runs tests of overlap to text significance of overlaps.
##USAGE: sbatch permutation_test_rbg.sbatch
##INSTRUCTION: User must upload exported gene lists from R to cluster & modify names of files.

hostname; date

#Go to directory for analysis
cd /home/FCAM/rgilmore/analysis/RNAseq/PWS_RNAseq/STAR/human2mouse

#Sort list of DEGs
sort -k1,1 smDEL_ortho_list.txt > smDELortho_sort.txt

for ((n=1; n<=1000; n++))
do
	#Shuffle lists, create gene lists using same number of genes as in analysis, sort files, join files and count for overlaps.
	shuf all_tested_mouse_ortho.txt | head -n 5348 - | sort -k1,1 - | join - smDELortho_sort.txt | wc -l >> permtest.txt
done
