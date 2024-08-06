#!/bin/bash
#SBATCH --job-name=fileConcat
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=32G
#SBATCH --mail-type=BEGIN,END
#SBATCH mail-user=rgilmore@uchc.edu
#SBATCH -o sbatch_%x_%j.out
#SBATCH -e sbatch_%x_%j.err

##DESCRIPTION: This is a loop to concatenate files from sequencing experiments.
##USAGE: sbatch fileConcatenate_rbg.sbatch
##INSTRUCTION: User must define variables up to end of user input.

#Create variable for list of files to be processed
SampleList="CT2-ngn2-s1 CT2-ngn2-s2 CT2-ngn2-s3 CT2-ngn2-s4 CT2-ngn2-s5 CT2-ngn2-s6 CT2-snord116-ic-ngn2-s1 CT2-snord116-ic-ngn2-s2 CT2-snord116-ic-ngn2-s3 CT2-snord116-ic-ngn2-s4 CT2-snord116-ic-ngn2-s5 CT2-snord116-ic-ngn2-s6 CT2-snord116-ngn2-s1 CT2-snord116-ngn2-s2 CT2-snord116-ngn2-s3 CT2-snord116-ngn2-s4 CT2-snord116-ngn2-s5 CT2-snord116-ngn2-s6 H9-snord116-ic-ngn2-s1 H9-snord116-ic-ngn2-s2 H9-snord116-ic-ngn2-s3 H9-snord116-ic-ngn2-s5"
#Set path to basespace directory
basePATH="/home/FCAM/rgilmore/basespace/"
#Set location of folder with symbolic links to files to be analyzed
filePATH="/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS2/ToBasespace"
#Designate where output data goes
outPATH="/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS2/fastq"

##################
##END USER INPUT##
##################

#Load module for fastqc
module load basemount/0.13.3.1573
module load pigz/2.2.3

if [ ! -d "$outPATH" ]; then
mkdir -p $outPATH
fi

#Mount to basespace
basemount ${basePATH}

for sample in ${SampleList}
do
	#List hostname & date for troubleshooting
		hostname
		date
	#Concatenate
		zcat ${filePATH}/${sample}*_R1* | pigz -p 16 > ${outPATH}/${sample}_R1.fastq.gz
		zcat ${filePATH}/${sample}*_R2* | pigz -p 16 > ${outPATH}/${sample}_R2.fastq.gz
done

#When finished, unmount basemount
basemount --unmount ${basePATH}

#Add wait command then a "Y" for lazy unmount question