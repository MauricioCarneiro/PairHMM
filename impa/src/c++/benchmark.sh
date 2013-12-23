#!/bin/bash

program="pairhmm-dev"
fileinp="performances/benchmark.in"
fileout="${fileinp%.in}.out"

echo "------------------------------------------------"
echo "Benchmarking the program: "$program
echo "With input case: "$fileinp
echo "Using as reference output: "$fileout
echo " "

(/usr/bin/time -f "%e" ./$program < $fileinp > fstdout) 2> ftime

precision=`../../tools/matcmp $fileout fstdout | grep "The biggest difference is: "`
precision=${precision#"The biggest difference is: "}
tiempo=`cat ftime`
#tiempo=`echo "$tiempo*1000" | bc`
hash=`git rev-parse HEAD | head -c 15`
rm -f fstdout ftime
echo $tiempo > performances/$hash
echo $precision >> performances/$hash

echo "Precision: $precision"
echo "Time: $tiempo"
echo " "

echo "------"

revisions=`git log --pretty=format:"%H"`
for rev in $revisions
do
	rev=`echo $rev | head -c 15`
	if [ -f performances/$rev ]; then
		echo "[$rev]: "`head -1 performances/$rev`"     "`tail -1 performances/$rev`
	fi
done

echo "------------------------------------------------"

