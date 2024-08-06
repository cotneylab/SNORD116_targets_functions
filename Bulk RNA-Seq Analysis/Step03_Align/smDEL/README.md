**Scripts for aligning smDEL data set.**
- For HISAT alignment:
  * Step 3.1 - Align using HISAT2.
  * Step 3.2 - Use samtools to sort .bam files.
- For STAR alignment:
  * Step 3.1 - Get appropriate annotation file to use in aligment.
  * Step 3.2 - Align using STAR.

*Note: There are differences between smDEL and lgDEL alignment scripts, not out of necessity, but rather because I developed as a coder. The lgDEL alignment is a much more high-throughput way of processing data, but either style will work.*
