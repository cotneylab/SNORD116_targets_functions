##DESCRIPTION: This is a loop to concatenate files from sequencing experiments.
##USAGE: ./MakeFileConcatenate_rbg.sh
##INSTRUCTION: User must define variables up to end of user input. Go into folder containing with script to run. Make script executeable before running.

#Create variable for list of files to be processed
export email="rgilmore@uchc.edu"
export threads=4
export memory="32G"
export fastqdir=/archive/labs/Cotney/DATA/RNA-Seq/ngn_neurons
export samplistPATH=/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3
export bash_location=/home/FCAM/rgilmore/.bashrc
export outPATH=/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3/fastq
export slurmDir=/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3/slurm

##################
##END USER INPUT##
##################

cat ${samplistPATH}/sample_list.txt | while read i
do
echo "#! /usr/bin/bash
#SBATCH --job-name=fileConcat_${i}_PWSRNAseq
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
#Concatenate
module load pigz/2.2.3
zcat ${fastqdir}/${i}*R1*.fastq.gz | pigz -p ${threads} > ${outPATH}/${i}_R1.fastq.gz
zcat ${fastqdir}/${i}*R2*.fastq.gz | pigz -p ${threads} > ${outPATH}/${i}_R2.fastq.gz" > ${slurmDir}/${i}_fileconcat_PWSRNAseq.slurm
done
