This folder contains all the codes for preprocessing and differential abundance analysis (DAA) of the plaque shotgun metagenomics data from the [early childhood dental caries](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9791751/) paper. The publicly available metadata can be found [here](https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=SRP331553&o=acc_s%3Aa).

1. `seqfetch` folder contains the code for downloading the sequence files from SRA.

2. `metasqueeze` folder contains the scripts that run [squeezemeta](https://github.com/jtamames/SqueezeMeta) pipeline for the assembly and binning of WGS reads.

3. `DAA_plaque_metag.R` load in the phyloseq data and apply nine different DAA methods (including ADAPT).

4. `plaque_WGS_summary.R` summarizes the DAA results and generate plots for Figure 3.