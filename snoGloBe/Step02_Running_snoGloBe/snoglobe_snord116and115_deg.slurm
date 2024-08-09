#!/bin/bash
#SBATCH --job-name=snoglobe_snord116
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cotney@uchc.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

source /home/FCAM/jcotney/.bashrc_miniconda3
conda activate snoglobe
export PATH=/home/FCAM/jcotney/TOOLS/snoGloBe/bin/:$PATH
cd /home/FCAM/jcotney/ANALYSIS/snoglobe
snoglobe --verbose -n 16 -t 0.98 -m -w 3 snord116_genes.fa ENSEMBLsharedsmlgDel.txt Homo_sapiens.GRCh38.88.chr.gtf ~/GENOME/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/ snord116_degshared
snoglobe --verbose -n 16 -t 0.98 -m -w 3 snord115_genes.fa ENSEMBLsharedsmlgDel.txt Homo_sapiens.GRCh38.88.chr.gtf ~/GENOME/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/ snord115_degshared
