##DESCRIPTION: This is script to create slurm scripts for aligning & sorting RNAseq reads in a high throughput manner.
##USAGE: ./makeHISAT2_align.sh
##INSTRUCTION: User must define variables up to end of user input. Go into folder containing with script to run. Make script executeable before running.

#Create variables
export email="rgilmore@uchc.edu"
export threads=16
export memory="64G"
export bash_location=/home/FCAM/rgilmore/.bashrc
export samplistPATH=/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3
#Designate where the raw data is
export inPATH=/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3/fastq
#Designate where processed data goes
export outPATH=/home/FCAM/rgilmore/analysis/RNAseq/HISAT2/NGNneu/PWS3
#Designate the genome to be used for alignment
export genome=/home/FCAM/jcotney/GENOME/Homo_sapiens/UCSC/hg38/Sequence/HISAT2Index/hg38
#Designate where slurm scripts go
export slurmDir=/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3/slurm

##################
##END USER INPUT##
##################

cat ${samplistPATH}/sample_list.txt | while read i
do
echo "#! /usr/bin/bash
#SBATCH --job-name=${i}_HISAT2align
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c ${threads}
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=${memory}
#SBATCH --mail-type=END
#SBATCH -o sbatch_%x_%j.out
#SBATCH -e sbatch_%x_%j.err
#List hostname & date for troubleshooting
hostname
date
#Load modulefiles
module load hisat2/2.2.1
module load samtools/1.9
#Align reads to genome
hisat2 -p ${threads} --fr --no-discordant --avoid-pseudogene --no-mixed -x ${genome} -1 ${inPATH}/${i}_R1.fastq.gz -2 ${inPATH}/${i}_R2.fastq.gz -S ${outPATH}/${i}.sam
#Convert sam file to sorted bams
samtools sort -@ ${threads} -o ${outPATH}/${i}.sorted.bam ${outPATH}/${i}.sam
#Remove sam file for storage reasons
rm ${outPATH}/${i}.sam" > ${slurmDir}/${i}_HISAT2align.slurm
done
