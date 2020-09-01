#!/bin/bash
#SBATCH -n 1
#SBATCH -t 12:0:0
#SBATCH --mem-per-cpu=32G
#SBATCH --account=mbarna
#SBATCH -o logs/%x.out
#SBATCH -e logs/%x.err

source /labs/mbarna/src/bashrc

unset PYTHONPATH
module add anaconda

WORKDIR=$LABDIR/covid_processed_data/20200808_233x

cd $WORKDIR

read1=$1
name=$2
outdir=$3
INDEX=index/233pool/233pool

source activate 062020

mkdir -p $outdir

#Dependencies used: cutadapt umi_tools bowtie2 samtools UMIcollapse

#cutadapt to remove Illumina adaptors at 3' end
cutadapt -j 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --discard-untrimmed --minimum-length=50 -o $outdir/$name".itrimmed.fastq.gz" $read1

#UMI_tools extract first 12 bases and add it to FASTQ record header
umi_tools extract --extract-method=string --bc-pattern=NNNNNNNNNNNN -I $outdir/$name".itrimmed.fastq.gz" -S $outdir/$name".itrimmed.umi.fastq.gz"

#cutadapt for constant regions in both 5' and 3' end (leaving 3 nt each side)
cutadapt -j 1 --discard-untrimmed --minimum-length=4 -g tcaccacggcgtgagatca...cgaagcggccgctctagaa -o $outdir/$name".itrimmed.umi.ctrimmed.fastq.gz" $outdir/$name".itrimmed.umi.fastq.gz"

#Bowtie2; should set some scoring later for better mismatch recovery; -L5 -N0; only reverse complement maps
bowtie2 -p 1 --end-to-end -L 11 -N 0 --nofw -x $INDEX -U $outdir/$name".itrimmed.umi.ctrimmed.fastq.gz" | samtools view -bS > $outdir/$name".bam"

#Sort and index bam
samtools sort $outdir/$name".bam" > $outdir/$name".sorted.bam"
samtools index $outdir/$name".sorted.bam"

#UMIcollapse to do deduplication 
./umicollapse bam -p 0.05 -i $outdir/$name".sorted.bam" -o $outdir/$name".sorted.dedup.bam"
samtools index $outdir/$name".sorted.dedup.bam"

#idxstats to make count table
samtools idxstats $outdir/$name".sorted.dedup.bam" > $outdir/$name".sorted.dedup.bam.idxstats"
