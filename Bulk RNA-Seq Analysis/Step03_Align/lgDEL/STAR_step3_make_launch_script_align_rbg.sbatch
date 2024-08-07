#Go to folder with slurm scripts from STAR_step1_make_align_scripts.sh
cd /home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3/slurm
#Make script to launch files
ls *STARalign.slurm | awk '{print "sbatch "$1}' > launch_STARalign.sh
#Make script executeable
chmod +x launch_STARalign.sh
#LAUNCH SCRIPT!!
./launch_STARalign.sh