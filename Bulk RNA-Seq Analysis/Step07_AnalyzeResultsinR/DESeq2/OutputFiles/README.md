Files generated from DESeq2 analysis in R.
- Files ending in ‘.RData’ are DESeq2 objects.
- Files ending in ‘.csv’ are results tables using various ```DESeqDataSetFromMatrix(design = )``` options in DESeq2 determined from Metadata file.

'Design' used is specified in file name. It is a formula to indicate how the counts for each gene depend on the variables in the Metadata file. (For more on this topic, please reference the DESeq2 manual: https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf)
- "CplusB" means the design was ```DESeqDataSetFromMatrix(design = ~ Genetic.Background + Condition + Genetic.Background:Condition)```
- "ConditionOnly" means the design was ```DESeqDataSetFromMatrix(design = ~ Condition)```
- All other files (*unspecified in file name*) means the design was design was ```DESeqDataSetFromMatrix(design = ~ Background_Condition)```
  * This was the design used for most downstream analyses. To use this design, create a column in Metadata file binding the "Genetic.Background" and "Condition" columns with a "_". (*Shown in .R code.*) \

&emsp;*The SharedGenesAnalysis, which focused on genes shared by both genetic backgrounds and deletion types, used DEGs that were statistically significant across all three design conditions.*

Aligners and deletion types are also specified in file names.

*Note: Input/output file names may not match exactly to what is listed in .R code file.*
