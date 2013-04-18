#!/bin/bash

hash=`git rev-parse HEAD | head -c 15`

if [ "$1" = "" ]; then
	echo ""
	echo "Usage: $(basename $0) <size of thread block in current version, WIDTHxHEIGHT>"
	echo "Example: $(basename $0) 32x3"
	echo ""
	exit 0
fi

for file in inp/*.in.cuda
do
	fn=$(basename "$file")
	fileref="${file%.in.cuda}.out.reference"
	fn=`date +"%Y.%m.%d-%H.%M"`"-cuda.${fn%.in.cuda}."`hostname`".$1.$hash"
	rm -f fout fstdout fstderr ftime
	(/usr/bin/time -f "%e" ./pairhmm-cuda $file fout > fstdout) 2> ftime
	o=`cat fstdout`
	t=`cat ftime`
	t=`echo "$t*1000" | bc`
	t="TOTAL (measured with 'time'): $t"
	e=`../tools/matcmp $fileref fout | grep "The biggest difference is: "`
	e=${e#"The biggest difference is: "}
	rm -f fout fstdout fstderr ftime
	echo "$o" > $fn
	echo "$t" >> $fn
	echo "BIGGEST ERROR: $e" >> $fn
	scp $fn eric@viete:/home/eric/Desktop/PairHMM/rep/PairHMM/impa/benchmarking
	rm $fn
done

