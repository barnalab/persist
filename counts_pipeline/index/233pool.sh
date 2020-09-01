awk -F"\t" ' { print ">"$1; print "AGC"$2"CGC"  } ' 233pool.txt > 233pool.fa
