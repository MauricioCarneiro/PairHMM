#!/bin/bash

fileinp="performances/benchmark.in.cuda"
fileout="${fileinp%.in.cuda}.out"

prog32="pairhmm32x1"
./$prog32 $fileinp fout32 > fstdout32
time32=`grep COMPUTATION: fstdout32`
time32=${time32#"COMPUTATION: "}
err32=`../../tools/matcmp $fileout fout32 | grep "The biggest difference is: "`
err32=${err32#"The biggest difference is: "}

prog64="pairhmm64x1"
./$prog64 $fileinp fout64 > fstdout64
time64=`grep COMPUTATION: fstdout64`
time64=${time64#"COMPUTATION: "}
err64=`../../tools/matcmp $fileout fout64 | grep "The biggest difference is: "`
err64=${err64#"The biggest difference is: "}

prog96="pairhmm96x1"
./$prog96 $fileinp fout96 > fstdout96
time96=`grep COMPUTATION: fstdout96`
time96=${time96#"COMPUTATION: "}
err96=`../../tools/matcmp $fileout fout96 | grep "The biggest difference is: "`
err96=${err96#"The biggest difference is: "}

hash=`git rev-parse HEAD | head -c 15`"-"`hostname`
echo "TIME32: "$time32 > performances/$hash
echo "ERR32: "$err32 >> performances/$hash
echo "TIME64: "$time64 >> performances/$hash
echo "ERR64: "$err64 >> performances/$hash
echo "TIME96: "$time96 >> performances/$hash
echo "ERR96: "$err96 >> performances/$hash

rm -f fout32 fstdout32 fout64 fstdout64 fout96 fstdout96

echo " "
echo "Input: "$fileinp
echo "Reference output: "$fileout
echo " "
echo "            +-----------------+-----------------+-----------------+"
printf "            | %-15s | %-15s | %-15s |\n" $prog32 $prog64 $prog96
echo "------------+-----------------+-----------------+-----------------+"

revisions=`git log --pretty=format:"%H"`
first=1
for rev in $revisions
do
	rev=`echo $rev | head -c 15`"-"`hostname`
	if [ -f performances/$rev ]; then
		time32=" "
		time32=`grep TIME32: performances/$rev`
		time32=${time32#"TIME32: "}
		err32=" "
		err32=`grep ERR32: performances/$rev`
		err32=${err32#"ERR32: "}
		time64=" "
		time64=`grep TIME64: performances/$rev`
		time64=${time64#"TIME64: "}
		err64=" "
		err64=`grep ERR64: performances/$rev`
		err64=${err64#"ERR64: "}
		time96=" "
		time96=`grep TIME96: performances/$rev`
		time96=${time96#"TIME96: "}
		err96=" "
		err96=`grep ERR96: performances/$rev`
		err96=${err96#"ERR96: "}
		printf "Time:       | %-15s | %-15s | %-15s | %s \n" $time32 $time64 $time96 $rev
		printf "Precision:  | %-15s | %-15s | %-15s |\n" $err32 $err64 $err96
		if [ "$first" -eq "1" ]; then
			first=0
			echo "============*=================*=================*=================*"
		else
			echo "------------+-----------------+-----------------+-----------------+"
		fi
	fi
done

