# openvaccine: PERSIST-seq pipeline  

This repo accompanies the following manuscript:  

Combinatorial optimization of mRNA structure, stability, and translation for RNA-based therapeutics  

Kathrin Leppek1*, Gun Woo Byeon1*, Wipapat Kladwang2*, Hannah K. Wayment-Steele3*, Craig H. Kerr1*, Adele F. Xu1**, Do Soon Kim2**, Ved V. Topkar4, Christian Choe5, Daphna Rothschild1, Gerald C. Tiu1, Roger Wellington-Oguri6, Kotaro Fujii1, Eesha Sharma2, Andrew M. Watkins2, John J. Nicol6, Jonathan Romano6,7, Bojan Tunguz2,8, Fernando Diaz9, Hui Cai9, Pengbo Guo9, Jiewei Wu9, Fanyu Meng9, Shuai Shi9, Eterna Participants6, Philip R. Dormitzer9, Alicia Solórzano9, Maria Barna1‡, Rhiju Das2‡  

It contains the code used to process PERSIST-seq data in the manuscript.  

## Dependencies
[cutadapt](https://github.com/marcelm/cutadapt)  
[umi_tools](https://github.com/CGATOxford/UMI-tools)  
[bowtie2](https://github.com/BenLangmead/bowtie2)  
[samtools](https://github.com/samtools)  
[UMIcollapse](https://github.com/Daniel-Liu-c0deb0t/UMICollapse)  


## Inputs
Sequencing reads deposited in GEO: [GSE173083](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173083)  
counts_pipeline/index/233pool.fa: FASTA containing barcode sequences for each template ID
analysis/constructs.tsv: synthetic mRNA design constructs
analysis/samples.csv: sample sheet


## Description  
1. counts_pipeline/p1.sh + join.sh: generates counts matrix (constructs x samples)  
2. analysis/persist_p.R: does a number of things from the sequencing data  
-calculates weighted ribosome number
-fits in-cell stability degradation coefficients
-using these two values, estimates expected protein levels according to kinetic model in the manuscript  
