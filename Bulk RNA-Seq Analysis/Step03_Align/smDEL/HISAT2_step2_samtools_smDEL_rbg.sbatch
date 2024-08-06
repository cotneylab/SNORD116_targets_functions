#! /usr/bin/bash
#SBATCH --job-name=Samtools_PWS2
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

##DESCRIPTION: This is script for sorting aligned RNAseq reads.
##USAGE: sbatch samtools_rbg.sbatch
##INSTRUCTION: User must define variables up to end of user input.

#Create variable for list of files to be processed
SampleList="CT2-ngn2-s1 CT2-ngn2-s2 CT2-ngn2-s3 CT2-ngn2-s4 CT2-ngn2-s5 CT2-ngn2-s6 CT2-snord116-ic-ngn2-s1 CT2-snord116-ic-ngn2-s2 CT2-snord116-ic-ngn2-s3 CT2-snord116-ic-ngn2-s4 CT2-snord116-ic-ngn2-s5 CT2-snord116-ic-ngn2-s6 CT2-snord116-ngn2-s1 CT2-snord116-ngn2-s2 CT2-snord116-ngn2-s3 CT2-snord116-ngn2-s4 CT2-snord116-ngn2-s5 CT2-snord116-ngn2-s6 H9-snord116-ic-ngn2-s1 H9-snord116-ic-ngn2-s2 H9-snord116-ic-ngn2-s5"
#Designate where the aligned data is
inSAM="/home/FCAM/rgilmore/analysis/RNAseq/HISAT2/NGNneu/PWS2/SAM"
#Designate where sorted data goes
outBAM="/home/FCAM/rgilmore/analysis/RNAseq/HISAT2/NGNneu/PWS2/BAM_sorted"

##################
##END USER INPUT##
##################

#Load modulefiles
module load samtools/1.9

if [ ! -d "$outBAM" ]; then
mkdir -p $outBAM
fi

for sample in ${SampleList}
do
	#List hostname & date for troubleshooting
		hostname
		date
	#Convert the sam files to bam files and sort
		echo Converting sam files to bam files on $sample data and sorting
		samtools sort -@ 16 -o ${outBAM}/${sample}.sorted.bam ${inSAM}/${sample}.sam
done
