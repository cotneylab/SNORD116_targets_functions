#Go to folder with slurm scripts from HISAT2_step1_make_align_scripts_rbg.sh
cd /home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3/slurm
#Make script to launch files
ls *HISAT2align.slurm | awk '{print "sbatch "$1}' > launch_HISAT2align.sh
#Make script executeable
chmod +x launch_HISAT2align.sh
#LAUNCH SCRIPT!!
./launch_HISAT2align.sh