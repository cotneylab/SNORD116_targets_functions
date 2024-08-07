**Scripts for aligning lgDEL data set.**
- For HISAT2 alignment:
  * Step 3.1 - Make individual slurm scripts for alignment of each sample.
  * Step 3.2 - Make shell script to launch individual slurm scripts concurrently.
- For STAR alignment:
  * Step 3.1 - Get annotation file (if necessary).
  * Step 3.2 - Make individual slurm scripts for alignment of each sample.
  * Step 3.3 - Make shell script to launch individual slurm scripts concurrently.

*Note: There are differences between smDEL and lgDEL alignment scripts, not out of necessity, but rather because I developed as a coder. The lgDEL alignment is a much more high-throughput way of processing data, but either style will work.*
