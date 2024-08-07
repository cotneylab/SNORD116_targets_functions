#! /usr/bin/bash
#SBATCH --job-name=RSeQC_HISAT2_PWSRNAseq
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

##DESCRIPTION: This is script for measuring gene body coverage in RNA-seq experiments from indexed bam files.
##USAGE: sbatch RSeQC_rbg.sbatch
##INSTRUCTION: User must define variables up to end of user input. User must create a text file with list of bams (& paths) if there are other bams in the folder that don't need to be analyzed.

#Designate path to bedfile 
bedFILE="/home/FCAM/rgilmore/tools/RSeQC/hg38_Gencode_V28.bed"
#Designate folder where the .bam & .bai files are (with trailing slash!)
inBAM="/home/FCAM/rgilmore/analysis/RNAseq/HISAT2/NGNneu/PWS3/bamsforRSeQC.txt"
#Designate where output goes
outDIR="/home/FCAM/rgilmore/analysis/RNAseq/HISAT2/NGNneu/PWS3/RSeQC"
#Set prefix for output
prefix="PWSRNAseq_HISAT2"

##################
##END USER INPUT##
##################

hostname; date

if [ ! -d "$outDIR" ]; then
mkdir -p $outDIR
fi

#Source conda
source ~/.bashrc_miniconda3
#Activate RSeQC environment
conda activate rseqc

#Go to directory where you want output to go
cd ${outDIR}

geneBody_coverage.py -r ${bedFILE} -i ${inBAM}  -o ${prefix}

#Deactivate conda
conda deactivate