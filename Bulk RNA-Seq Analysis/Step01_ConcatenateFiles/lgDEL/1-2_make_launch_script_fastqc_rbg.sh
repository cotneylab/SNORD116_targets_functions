#Go to folder with slurm output scripts from make_highthroughput_fastqc_scripts.sh
cd /home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3/slurm
#Make script to launch files
ls *fastqc*.slurm | awk '{print "sbatch "$1}' > launch_fastqc.sh
#Make script executeable
chmod +x launch_fastqc.sh
#LAUNCH SCRIPT!!
./launch_fastqc.sh