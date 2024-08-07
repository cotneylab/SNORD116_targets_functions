#! /usr/bin/bash
#SBATCH --job-name=BedtoBigwig_PWS
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mem=32G
#SBATCH --mail-user=rgilmore@uchc.edu
#SBATCH -o sbatch_%x_%j.out
#SBATCH -e sbatch_%x_%j.err

##DESCRIPTION: This is script for sorting bedgraph files & making them into bigwigs for visualization on UCSC browser.
##USAGE: sbatch BedtoBigwig_rbg.sbatch
##INSTRUCTION: User must define variables up to end of user input.

#Create variable for list of files to be processed
SampleList="CT2-ngn2-s1 CT2-ngn2-s2 CT2-ngn2-s3 CT2-ngn2-s4 CT2-ngn2-s5 CT2-ngn2-s6 CT2-snord116-ic-ngn2-s1 CT2-snord116-ic-ngn2-s2 CT2-snord116-ic-ngn2-s3 CT2-snord116-ic-ngn2-s4 CT2-snord116-ic-ngn2-s5 CT2-snord116-ic-ngn2-s6 CT2-snord116-ngn2-s1 CT2-snord116-ngn2-s2 CT2-snord116-ngn2-s3 CT2-snord116-ngn2-s4 CT2-snord116-ngn2-s5 CT2-snord116-ngn2-s6 H9-snord116-ic-ngn2-s1 H9-snord116-ic-ngn2-s2 H9-snord116-ic-ngn2-s5"
#Designate where the unsorted bedgraphs are
inPATH="/home/FCAM/rgilmore/analysis/RNAseq/STAR/PWS_Exp/BAM_sorted"
#Designate where bigwigs go
outPATH="/tgc/TGCore_User_Data/WebData/cotney/hubs/gilmore/PWS_RNAseq"
#Designate the chromosome file (& path to it)
chrFile="/home/FCAM/jcotney/GENOME/Homo_sapiens/UCSC/hg38/Sequence/STARIndex/chrNameLength.txt"

##################
##END USER INPUT##
##################

if [ ! -d "$outPATH" ]; then
mkdir -p $outPATH
fi

#Load modulefile
module load kent-tools/1.04.00

for sample in ${SampleList}
do
	#List hostname & date for troubleshooting
	hostname
	date
	#Print to screen for script update and name confirmation
        echo Sorting $sample bedgraph
	#Sort bedgraph
		sort -k1,1 -k2,2n ${inPATH}/${sample}_Signal.Unique.str1.out.bg > ${inPATH}/${sample}.Unique.str1.sorted.out.bg
		sort -k1,1 -k2,2n ${inPATH}/${sample}_Signal.Unique.str2.out.bg > ${inPATH}/${sample}.Unique.str2.sorted.out.bg
	#Convert sorted bedgraph to bigwig
		echo Converting $sample bedgraph to bigwig
		bedGraphToBigWig ${inPATH}/${sample}.Unique.str1.sorted.out.bg ${chrFile} ${outPATH}/${sample}.str1.bw
		bedGraphToBigWig ${inPATH}/${sample}.Unique.str2.sorted.out.bg ${chrFile} ${outPATH}/${sample}.str2.bw
done
