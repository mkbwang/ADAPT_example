This folder contains the codes for preprocessing and differential abundance analysis(DAA) of 16S sequencing of saliva samples from the [early childhood dental caries](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9791751/) paper. The publicly available metadata can be found [here](https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=SRP331553&o=acc_s%3Aa).

1. `seqfetch` folder contains the script for downloading the sequencing files.

2. `DADA2` folder contains the scripts for trimming, denoising and taxonomy assignment. The reference 16S rRNA sequences can be downloaded from [HOMD](https://www.homd.org/download#refseq). It also contains the codes for decontamination and sample filtering.

3. `DAA_16S_Saliva_month12.R` carry out DAA with nine methods (including ADAPT).

4. `16S_saliva_12month_summary.R` generates the plots for Figure 3.