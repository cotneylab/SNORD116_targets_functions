#Start interactive session
srun --pty -p general --qos=general --mem=16G -c 4 bash
#Go to folder with sorted bam files
cd /home/FCAM/rgilmore/analysis/RNAseq/HISAT2/NGNneu/PWS2/BAM_sorted
#Write script to index bam files (.bai is required input format)
ls *.bam | awk '{print "samtools index -@ 4 "$1""}' > index.sh
#Make it executeable
chmod +x index.sh
#Load samtools
module load samtools/1.9
#Run script
./index.sh