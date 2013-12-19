#!/bin/bash

program="$1"
fileinp="performances/benchmark.in"
fileout="${fileinp%.in}.out"

echo "------------------------------------------------"
echo "Benchmarking the program: "$program
echo "With input case: "$fileinp
echo "Using as reference output: "$fileout
echo " "

./$program < $fileinp > fstdout 2> ftime
precision=`../../tools/matcmp $fileout fstdout | grep "The biggest difference is: "`
precision=${precision#"The biggest difference is: "}
tiempo=`cat ftime`
tiempo=`echo "$tiempo*1000" | bc`
hash=`git rev-parse HEAD | head -c 15`
modified=`git status --porcelain | grep -m 1 "^ M" | sed "s/^ M.*/-modified/g"`
rm -f fstdout ftime
echo $tiempo > performances/$hash$modified
echo $precision >> performances/$hash$modified

echo "Precision: $precision"
echo "Time: $tiempo"
echo " "

echo "------"

revisions=`git log --pretty=format:"%H"`
for rev in $revisions
do
	rev=`echo $rev | head -c 15`
	if [ -f performances/$rev-modified ]; then
		echo "[$rev+]: "`head -1 performances/$rev`"     "`tail -1 performances/$rev`
	fi
	if [ -f performances/$rev ]; then
		echo "[$rev ]: "`head -1 performances/$rev`"     "`tail -1 performances/$rev`
	fi
done

echo "------------------------------------------------"

