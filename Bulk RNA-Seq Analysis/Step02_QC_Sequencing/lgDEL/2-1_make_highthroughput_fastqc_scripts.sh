##DESCRIPTION: This is a loop to run fastqc in a high throughput manner.
##USAGE: ./fastqc_rbg.sh
##INSTRUCTION: User must define variables up to end of user input. Go into folder containing with script to run. Make script executeable before running.

#Create variable for list of files to be processed
export email="rgilmore@uchc.edu"
export threads=8
export memory="32G"
export bash_location=/home/FCAM/rgilmore/.bashrc
export samplistPATH=/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3
export filePATH=/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3/fastq
export outPATH=/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3/fastqc
export slurmDir=/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3/slurm

##################
##END USER INPUT##
##################

cat ${samplistPATH}/sample_list.txt | while read i
do
echo "#! /usr/bin/bash
#SBATCH --job-name=fastqc_${i}_RNAseq
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
#Load module for fastqc
module load fastqc/0.11.7
#Run fastqc
fastqc -o ${outPATH} -t ${threads} ${filePATH}/${i}*.fastq.gz" > ${slurmDir}/${i}_fastqc_PWSRNAseq.slurm
done
