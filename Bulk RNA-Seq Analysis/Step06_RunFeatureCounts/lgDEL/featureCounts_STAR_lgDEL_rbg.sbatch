#! /usr/bin/bash
#SBATCH --job-name=featureCounts_RNAseqSTAR
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mem=64G
#SBATCH --mail-user=rgilmore@uchc.edu
#SBATCH -o sbatch_%x_%j.out
#SBATCH -e sbatch_%x_%j.err

##DESCRIPTION: This is script for extracting counts with feature counts using aligned RNAseq reads converted into .bam files.
##USAGE: sbatch featureCounts_rbg.sbatch
##INSTRUCTION: User must define variables up to end of user input.

#Designate path to gtf file
gtfFile="/home/FCAM/rgilmore/annotations/gencode.v25.annotation.gtf"
#Designate where the data is
inBAM="/home/FCAM/rgilmore/analysis/RNAseq/STAR/PWS3/alignment"
#Designate where extracted counts go
outCOUNTS="/home/FCAM/rgilmore/analysis/RNAseq/STAR/PWS3/featureCounts"
#Name for output file
name="PWSRNAseq_STAR"

##################
##END USER INPUT##
##################

hostname; date

#Load modulefiles
module load subread/2.0

if [ ! -d "$outCOUNTS" ]; then
mkdir -p $outCOUNTS
fi

echo Performing featureCounts on sorted bam files
featureCounts -p -s 2 -T 4 -t exon -g gene_id -a ${gtfFile} -o ${outCOUNTS}/${name}_featureCounts.txt ${inBAM}/*Aligned.sortedByCoord.out.bam