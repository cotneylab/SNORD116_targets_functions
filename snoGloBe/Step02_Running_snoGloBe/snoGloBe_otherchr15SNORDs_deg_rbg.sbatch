#! /usr/bin/bash
#SBATCH --job-name=snoGloBe_other
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=BEGIN,END
#SBATCH --mem=64G
#SBATCH --mail-user=rgilmore@uchc.edu
#SBATCH -o %x_%j_sbatch.out
#SBATCH -e %x_%j_sbatch.err

##DESCRIPTION: This is a script to run the snoGloBe snoRNA target predictor tool from the Scott Group (https://github.com/scottgroup/snoGloBe).
##USAGE: sbatch snoGloBe_rbg.sbatch

#Source conda
source /home/FCAM/jcotney/.bashrc_miniconda3
#Activate snoglobe
conda activate snoglobe
#Set paths
export PATH=/home/FCAM/jcotney/TOOLS/snoGloBe/bin/:$PATH

#Go into folder for analysis
cd /home/FCAM/rgilmore/analysis/snoGloBe
#Run snoglobe
snoglobe --verbose -n 16 -t 0.98 -m -w 3 chr15otherSNORDs.fa ENSEMBLsharedsmlgDel.txt Homo_sapiens.GRCh38.88.chr.gtf /home/FCAM/jcotney/GENOME/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/ chr15otherSNORDs_degshared

#Deactivate snoglobe
conda deactivate
