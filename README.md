# Ladder-Seq benchmark

![Ladder-Seq benchmark](/Benchmark700.png)

This repository contains [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows to reproduce the benchmarking results in the Ladder-seq paper. There are workflows for (i) transcript quantification 
(ii) reference-based transcript assembly, and (iii) *de novo* assembly. 

Each pipeline simulates RNA-seq reads using program [RSEM](https://deweylab.github.io/RSEM/) (step 2) from a ground truth transcriptome with abundances and error profile
calculated by RSEM from GEUVADIS data set NA12716_7 in step 1. A matching Ladder-seq sample is obtained by separating reads in silico into 7 bands based on 
probability mass functions estimated from a mouse NPC Ladder-seq sample (step 3). 
Transcripts are quantified and assembled by our Ladder-seq tailored transcript analysis methods **kallisto-ls**, **StringTie-ls**, and **Trinity-ls** from the Ladder-seq samples, 
while their conventional counterparts are run on the corresponding RNA-seq sample (step 4). The results are compared to the ground truth to evaluate and compare their 
accuracy (step 5).

## Running the workflow
All binaries the workflows use are included in subdirectories ``ext``. Only [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) needs to be built before running the *de novo* assembly workflow and the path 
to its binary adjusted in the corresponding ``config.py`` file. 
After changing to one of the main directories that correspond to worklows (i)-(iii), run the corresponding Snakemake workflow:
```shell
snakemake -s Snakemake_runPipeline.smk -j 8
```
where parameter '-j' specifies the number of CPU cores to be used. 


