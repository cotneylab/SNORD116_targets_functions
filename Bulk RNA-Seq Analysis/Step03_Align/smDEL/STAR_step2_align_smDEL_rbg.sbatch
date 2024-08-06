#! /usr/bin/bash
#SBATCH --job-name=STAR_align_PWS2
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
##USAGE: sbatch STAR_align.sbatch
##INSTRUCTION: User must define variables up to end of user input.

#Create variable for list of files to be processed
SampleList="CT2-ngn2-s1 CT2-ngn2-s2 CT2-ngn2-s3 CT2-ngn2-s4 CT2-ngn2-s5 CT2-ngn2-s6 CT2-snord116-ic-ngn2-s1 CT2-snord116-ic-ngn2-s2 CT2-snord116-ic-ngn2-s3 CT2-snord116-ic-ngn2-s4 CT2-snord116-ic-ngn2-s5 CT2-snord116-ic-ngn2-s6 CT2-snord116-ngn2-s1 CT2-snord116-ngn2-s2 CT2-snord116-ngn2-s3 CT2-snord116-ngn2-s4 CT2-snord116-ngn2-s5 CT2-snord116-ngn2-s6 H9-snord116-ic-ngn2-s1 H9-snord116-ic-ngn2-s2 H9-snord116-ic-ngn2-s5"
#Designate where the raw data is
inPATH="/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS2/fastq"
#Designate where processed data goes
outPATH="/home/FCAM/rgilmore/analysis/RNAseq/STAR/PWS_Exp/BAM_sorted"
#Designate the genome to be used for alignment
genome="/home/FCAM/jcotney/GENOME/Homo_sapiens/UCSC/hg38/Sequence/STARIndex"
#Designate path to gtf file
gtfFile="/home/FCAM/rgilmore/annotations/gencode.v25.annotation.gtf"

##################
##END USER INPUT##
##################

if [ ! -d "$outPATH" ]; then
mkdir -p $outPATH
fi

#Load modulefile
module load module load STAR/2.7.1a

for sample in ${SampleList}
do
	#List hostname & date for troubleshooting
	hostname
	date
	#Print to screen for script update and name confirmation
        echo Aligning reads from $sample to genome
	#Align reads to genome
		STAR --runThreadN 16 --outTmpDir ~/STARtemp --genomeDir ${genome} \
		--readFilesIn ${inPATH}/${sample}_R1.fastq.gz ${inPATH}/${sample}_R2.fastq.gz \
		--outFileNamePrefix ${outPATH}/${sample}_ --readFilesCommand zcat --sjdbGTFfile ${gtfFile} \
		--outSAMprimaryFlag OneBestScore --outFilterType BySJout --outFilterMultimapNmax 1 --alignSJoverhangMin 8 \
		--alignSJDBoverhangMin 2 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
		--outSAMtype BAM SortedByCoordinate --outSAMattributes Standard \
		--outWigType bedGraph --outWigStrand Stranded --outWigNorm RPM
	#Delete temporary directory
		rm -r ~/STARtemp
done
