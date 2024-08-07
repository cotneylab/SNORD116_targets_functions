#! /usr/bin/bash
#SBATCH --job-name=RSeQC_PWSexp
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

##DESCRIPTION: This is script for measuring gene body coverage in RNA-seq experiments from .bam files.
##USAGE: sbatch RSeQC_rbg.sbatch
##INSTRUCTION: User must define variables up to end of user input.

#Designate path to bedfile 
bedFILE="/home/FCAM/rgilmore/tools/RSeQC/hg38_Gencode_V28.bed"
#Designate folder where the bam files are (with trailing slash!)
inBAM="/home/FCAM/rgilmore/analysis/RNAseq/HISAT2/NGNneu/PWS2/BAM_sorted/"
#Designate where output goes
outDIR="/home/FCAM/rgilmore/analysis/RNAseq/HISAT2/NGNneu/PWS2/RSeQC"
#Set prefix for output
prefix="PWS2_Exp"

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