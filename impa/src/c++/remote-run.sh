#!/bin/bash

if [ $# -ne 2 ] 
then
 	echo " "
	echo "Usage example:"
	echo "    $0 \"germain bresenham\" \"small medium\""
	echo " "
	exit 0
fi

make
IFS=' ' read -a testcases <<< "$2"
Host=$1

for testcase in "${testcases[@]}"
do
	ssh $Host rm -rf pairhmm-1-base pairhmm-2-omp pairhmm-3-hybrid-float-double pairhmm-4-hybrid-diagonal pairhmm-5-hybrid-diagonal-homogeneus pairhmm-6-onlythreediags matcmp reports/ firstnotacceptable
	ssh $Host mkdir ~/reports
	scp pairhmm-1-base $Host:~/
	scp pairhmm-2-omp $Host:~/
	scp pairhmm-3-hybrid-float-double $Host:~/
	scp pairhmm-4-hybrid-diagonal $Host:~/
	scp pairhmm-5-hybrid-diagonal-homogeneus $Host:~/
	scp pairhmm-6-onlythreediags $Host:~/
	scp ../../tools/matcmp $Host:~/
	scp ../../tools/firstnotacceptable $Host:~/
	inp="test_data/"$testcase"Test.in"
	ref="test_data/"$testcase"Test.out"
	#for bin in "pairhmm-1-base" "pairhmm-2-omp" "pairhmm-3-hybrid-float-double" "pairhmm-4-hybrid-diagonal" "pairhmm-5-hybrid-diagonal-homogeneus" "pairhmm-6-onlythreediags" ;
	for bin in "pairhmm-1-base" ; #"pairhmm-2-omp" "pairhmm-3-hybrid-float-double" ;
	#"pairhmm-4-hybrid-diagonal" "pairhmm-5-hybrid-diagonal-homogeneus" "pairhmm-6-onlythreediags" ;
	do
		o="reports/"$testcase"-$bin-$Host.out"
		e="reports/"$testcase"-$bin-$Host.err"
		p="reports/"$testcase"-$bin-$Host.precision"
		ssh $Host "~/$bin < $inp > $o 2> $e"
		ssh $Host "~/matcmp $ref $o > $p"
	done
	scp $Host:~/reports/*.err .
	scp $Host:~/reports/*.precision .
done
