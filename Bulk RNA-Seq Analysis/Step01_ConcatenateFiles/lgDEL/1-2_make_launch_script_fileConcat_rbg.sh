#Go to folder with slurm scripts from MakeFileConcatenate_lgDEL_rbg.sh
cd /home/FCAM/rgilmore/data/rna-seq/NGNneurons/PWS3/slurm
#Make script to launch files
ls *fileconcat*.slurm | awk '{print "sbatch "$1}' > launch_fileConcat.sh
#Make script executeable
chmod +x launch_fileConcat.sh
#LAUNCH SCRIPT!!
./launch_fileConcat.sh