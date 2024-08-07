#Go to folder with slurm scripts from 4-1_makeIndexedbams_scripts_rbg.sh
cd /home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3/slurm
#Make script to launch files
ls *indexbam.slurm | awk '{print "sbatch "$1}' > launch_indexbam.sh
#Make script executeable
chmod +x launch_indexbam.sh
#LAUNCH SCRIPT!!
./launch_indexbam.sh