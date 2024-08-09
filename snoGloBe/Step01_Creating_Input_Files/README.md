# Step 1 - Create Input Files for snoGloBe #
- Step 1.1 - Retrieve ENSEMBL gtf file.
- Step 1.2 - Parse the gtf file for snoRNAs and create a bed file.
- Step 1.3 - Use bedtools to convert the .bed file into a fasta (.fa) file.

An input gene list is also required. We used the list of 42 shared dysregulated genes between smDEL and lgDEL datasets.
However, any list of ENSEMBL gene identifiers should work.
While ideally one would examine predicted interactions with all genes, the authors admit doing so for even one snoRNA is time consuming and computationally intensive, which we found to be the case as well, so using a subset of genes is recommended.
