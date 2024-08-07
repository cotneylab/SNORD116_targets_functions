**RSeQC analysis for lgDEL alignments.**
- Step 4.1 - Use shell script to make individual slurm scripts for each sample (.bam) to be indexed. (Creating indexed BAM files - .bai)
- Step 4.2 - Create & launch shell script to launch individual slurm scripts concurrently.
- Step 4.3 - Run RSeQC on .bai files. (This step can take a while; plan accordingly.)

*Code from one representative script is shown. Outputs can be found in ../Outputs/.*
