#!/bin/bash
headers="id"
cmd="paste <(cut -f1 out/fraction-spikein-1-GAATCCGT.sorted.dedup.bam.idxstats)"
for file in `ls out/*.sorted.dedup.bam.idxstats`
do 
	name=$(cut -f1 -d"." <(basename $file))
	headers+="\t"$name; cmd+=" <(cut -f3 $file)"
done
echo -e $headers > joined.txt
eval $cmd >> joined.txt
