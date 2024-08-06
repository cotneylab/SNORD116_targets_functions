#! /usr/bin/bash
#SBATCH --job-name=HISAT2_align_WTvSNORDmut
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mem=64G
#SBATCH --mail-user=rgilmore@uchc.edu
#SBATCH -o sbatch_%x_%j.out
#SBATCH -e sbatch_%x_%j.err

##DESCRIPTION: This is script for aligning RNAseq reads.
##USAGE: sbatch HISAT2_align.sbatch
##INSTRUCTION: User must define variables up to end of user input.

#Create variable for list of files to be processed
SampleList="CT2-ngn2-s1 CT2-ngn2-s2 CT2-ngn2-s3 CT2-ngn2-s4 CT2-ngn2-s5 CT2-ngn2-s6 CT2-snord116-ic-ngn2-s1 CT2-snord116-ic-ngn2-s2 CT2-snord116-ic-ngn2-s3 CT2-snord116-ic-ngn2-s4 CT2-snord116-ic-ngn2-s5 CT2-snord116-ic-ngn2-s6 CT2-snord116-ngn2-s1 CT2-snord116-ngn2-s2 CT2-snord116-ngn2-s3 CT2-snord116-ngn2-s4 CT2-snord116-ngn2-s5 CT2-snord116-ngn2-s6 H9-snord116-ic-ngn2-s1 H9-snord116-ic-ngn2-s2 H9-snord116-ic-ngn2-s5"
#Designate where the raw data is
inPATH="/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS2/fastq"
#Designate where processed data goes
outPATH="/home/FCAM/rgilmore/analysis/RNAseq/HISAT2/NGNneu/PWS2/SAM"
#Designate the genome to be used for alignment
genome="/home/FCAM/jcotney/GENOME/Homo_sapiens/UCSC/hg38/Sequence/HISAT2Index/hg38"

##################
##END USER INPUT##
##################

if [ ! -d "$outPATH" ]; then
mkdir -p $outPATH
fi

#Load modulefile
module load hisat2/2.2.1

for sample in ${SampleList}
do
	#List hostname & date for troubleshooting
	hostname
	date
	#Print to screen for script update and name confirmation
        echo Aligning reads from $sample to genome
	#Align reads to genome
		hisat2 -p 16 --fr --no-discordant --avoid-pseudogene --no-mixed -x ${genome} -1 ${inPATH}/${sample}_R1.fastq.gz -2 ${inPATH}/${sample}_R2.fastq.gz -S ${outPATH}/${sample}.sam
done
