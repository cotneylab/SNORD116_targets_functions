#! /usr/bin/bash
#SBATCH --job-name=PermTest_smDELwlgDEL
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

##DESCRIPTION: This is script is a loop which runs tests of overlap to text significance of intersects.
##USAGE: sbatch permutation_test_rbg.sbatch
##INSTRUCTION: User must upload exported gene lists from R to cluster & modify names of files.

hostname; date

#Go to directory for analysis
cd /home/FCAM/rgilmore/analysis/RNAseq/STAR/PWS3/permutationTest

#Sort lists
sort -k1,1 allSigGene_3analyses_smDEL.txt > smDELsigGene_sort.txt
sort -k1,1 allSigGene_3analyses_lgDEL.txt > lgDELsigGene_sort.txt

for ((n=1; n<=1000; n++))
do
	#Shuffle lists, create gene lists using same number of genes as in analysis, sort files, join files and count for overlaps.
	shuf smDELsigGene_permtest.txt | head -n 232 - | sort -k1,1 - | join - lgDELsigGene_sort.txt | wc -l >> smtolg.txt
	shuf lgDELsigGene_permtest.txt | head -n 691 - | sort -k1,1 - | join - smDELsigGene_sort.txt | wc -l >> lgtosm.txt
done
