##DESCRIPTION: This is script to create slurm scripts for aligning RNAseq reads in a high throughput manner.
##USAGE: ./STAR_align.sh
##INSTRUCTION: User must define variables up to end of user input. Go into folder containing with script to run. Make script executeable before running.

#Create variables
export email="rgilmore@uchc.edu"
export threads=16
export memory="64G"
export bash_location=/home/FCAM/rgilmore/.bashrc
#Designate where the sample list is
export samplistPATH=/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3
#Designate where the raw data is
export inPATH=/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3/fastq
#Designate where processed data goes
export outPATH=/home/FCAM/rgilmore/analysis/RNAseq/STAR/PWS3
#Designate the genome to be used for alignment
export genome=/home/FCAM/jcotney/GENOME/Homo_sapiens/UCSC/hg38/Sequence/STARIndex
#Designate path to gtf file
export gtfFile=/home/FCAM/rgilmore/annotations/gencode.v25.annotation.gtf
#Designate where slurm scripts go
export slurmDir=/home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3/slurm

##################
##END USER INPUT##
##################

cat ${samplistPATH}/sample_list.txt | while read i
do
echo "#! /usr/bin/bash
#SBATCH --job-name=${i}_STARalign
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
module load STAR/2.7.1a
#Align reads to genome
STAR --runThreadN ${threads} --outTmpDir ~/STARtemp_${i} --genomeDir ${genome} \\
--readFilesIn ${inPATH}/${i}_R1.fastq.gz ${inPATH}/${i}_R2.fastq.gz \\
--outFileNamePrefix ${outPATH}/${i}_ --readFilesCommand zcat --sjdbGTFfile ${gtfFile} \\
--outSAMprimaryFlag OneBestScore --outFilterType BySJout --outFilterMultimapNmax 1 --alignSJoverhangMin 8 \\
--alignSJDBoverhangMin 2 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \\
--outSAMtype BAM SortedByCoordinate --outSAMattributes Standard \\
--outWigType bedGraph --outWigStrand Stranded --outWigNorm RPM
#Delete temporary directory
rm -r ~/STARtemp_${i}" > ${slurmDir}/${i}_STARalign.slurm
done
