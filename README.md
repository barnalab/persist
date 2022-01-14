# openvaccine: PERSIST-seq

This repo accompanies the following manuscript:  

Combinatorial optimization of mRNA structure, stability, and translation for RNA-based therapeutics  

Kathrin Leppek, Gun Woo Byeon, Wipapat Kladwang, Hannah K. Wayment-Steele, Craig H. Kerr, Adele F. Xu, Do Soon Kim, Ved V. Topkar, Christian Choe, Daphna Rothschild, Gerald C. Tiu, Roger Wellington-Oguri, Kotaro Fujii, Eesha Sharma, Andrew M. Watkins, John J. Nicol, Jonathan Romano, Bojan Tunguz, Fernando Diaz, Hui Cai, Pengbo Guo, Jiewei Wu, Fanyu Meng, Shuai Shi, Eterna Participants, Philip R. Dormitzer, Alicia Solórzano, Maria Barna, Rhiju Das  

It contains the code used to process PERSIST-seq data in the manuscript.  

## Dependencies
[cutadapt](https://github.com/marcelm/cutadapt)  
[umi_tools](https://github.com/CGATOxford/UMI-tools)  
[bowtie2](https://github.com/BenLangmead/bowtie2)  
[samtools](https://github.com/samtools)  
[UMIcollapse](https://github.com/Daniel-Liu-c0deb0t/UMICollapse)  
[R](https://www.r-project.org/)
[errors](https://cran.r-project.org/web/packages/errors/index.html)
[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)  


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
