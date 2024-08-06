#!/bin/bash
#SBATCH --job-name=fastqc_PWS2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,END
#SBATCH mail-user=rgilmore@uchc.edu
#SBATCH -o sbatch_%x_%j.out
#SBATCH -e sbatch_%x_%j.err

##DESCRIPTION: This is a loop to run fastqc.
##USAGE: sbatch fastqc_rbg.sbatch
##INSTRUCTION: User must define variables up to end of user input.

#Create variable for list of files to be processed
SampleList="CT2-ngn2-s1 CT2-ngn2-s2 CT2-ngn2-s3 CT2-ngn2-s4 CT2-ngn2-s5 CT2-ngn2-s6 CT2-snord116-ic-ngn2-s1 CT2-snord116-ic-ngn2-s2 CT2-snord116-ic-ngn2-s3 CT2-snord116-ic-ngn2-s4 CT2-snord116-ic-ngn2-s5 CT2-snord116-ic-ngn2-s6 CT2-snord116-ngn2-s1 CT2-snord116-ngn2-s2 CT2-snord116-ngn2-s3 CT2-snord116-ngn2-s4 CT2-snord116-ngn2-s5 CT2-snord116-ngn2-s6 H9-snord116-ic-ngn2-s1 H9-snord116-ic-ngn2-s2 H9-snord116-ic-ngn2-s3 H9-snord116-ic-ngn2-s5"
#Set location of files to be analyzed
filePATH="/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS2/fastq"
#Designate where output data goes
outPATH="/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS2/fastqc"

##################
##END USER INPUT##
##################

#Load module for fastqc
module load fastqc/0.11.7

if [ ! -d "$outPATH" ]; then
mkdir -p $outPATH
fi

for sample in ${SampleList}
do
	#List hostname & date for troubleshooting
		hostname
		date
	#Run fastqc
		fastqc -o ${outPATH} -t 8 ${filePATH}/${sample}*.fastq.gz
done
