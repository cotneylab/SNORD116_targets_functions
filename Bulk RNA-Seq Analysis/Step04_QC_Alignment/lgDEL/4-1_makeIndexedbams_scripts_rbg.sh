##DESCRIPTION: This is script to create slurm scripts for indexing sorted bams (creates .bai) in a high throughput manner for use with RSeQC.
##USAGE: ./makeIndexedbams.sh
##INSTRUCTION: User must define variables up to end of user input. Go into folder containing with script to run. Make script executeable before running.

#Create variables
export email="rgilmore@uchc.edu"
export threads=4
export memory="32G"
export bash_location=/home/FCAM/rgilmore/.bashrc
#Designate where sample list is
export samplistPATH=/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3
#Designate where sorted bams are & where indexed bams go
export bamPATH=/home/FCAM/rgilmore/analysis/RNAseq/HISAT2/NGNneu/PWS3
#Designate where slurm scripts go
export slurmDir=/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3/slurm

##################
##END USER INPUT##
##################

cat ${samplistPATH}/sample_list.txt | while read i
do
echo "#! /usr/bin/bash
#SBATCH --job-name=${i}_indexbam
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
#Load modulefile
module load samtools/1.9
#Go into folder containing bams
cd ${bamPATH}
#Index sorted bams
samtools index -@ ${threads} ${i}.sorted.bam" > ${slurmDir}/${i}_indexbam.slurm
done
