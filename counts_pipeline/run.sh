for file in `ls fastq/*.fastq.gz`
do 

	r=$(cut -d"_" -f1 <(echo $file))
	name=$(basename $r)
	out="out"

	echo $file $name
	p1=$(sbatch --parsable --job-name=$name"_p1" ./p1.sh $file $name $out)
done
