# 02. Alignment

This [**snakefile_alignPE**]() takes FASTQ files, aligns them to a reference genome, indexes the resulting bams, merges sequencing replicates, and generates resulting read depth information for each sample. It simultaneously runs an IBS (Identity-by-State) analysis on the replicates to ensure that only true replicates are being merged and also to identify possible replicate individuals that were mislabeled as separate individuals (i.e. a recaptured individual that was considered new and given a new field ID).  
 <img align="left" src="Images/01.Alignment-DAG.png" width="300"> 


Note --> You shouldn’t need to modify anything in this file! Just modify the config file below.

[**snakefile_alignPE.yml**]()

Part of this pipeline also relies on an accessory R script -> [**compare_replicates.R**](https://github.com/UCDavis-MECU/Genomic-Methods/blob/main/01.Alignment-SNP-Calling/GBS-RAD/compare_replicates.R) \
You'll want to make sure you download that as well and keep it in the same directory as your snakemake and .yml files.
