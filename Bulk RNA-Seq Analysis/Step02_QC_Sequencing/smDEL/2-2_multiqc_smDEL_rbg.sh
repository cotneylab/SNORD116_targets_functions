#Start interactive session
srun --pty -p general --qos=general --mem=16G -c 4 bash
#Load multiqc module
module load MultiQC/1.10.1
#Go to folder with fastqc files
cd /home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS2/fastqc
#Run multiqc
multiqc .