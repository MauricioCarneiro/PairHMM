#!/bin/bash

test()
{
	program=$1
	fileinp="performances/benchmark.in"
	fileout="${fileinp%.in}.out"

	(/usr/bin/time -f "%e" ./$program < $fileinp > fstdout) 2> ftime
	e=`../../tools/matcmp $fileout fstdout | grep "The biggest difference is: "`
	e=${e#"The biggest difference is: "}
	t=`cat ftime`
	#t=`echo "$t*1000" | bc`
	rm -f fstdout ftime
	printf "%40s %10s %10s\n" $program $e $t
}

make clean
make pairhmm-7-presse
#test "pairhmm-1-base" 
#test "pairhmm-2-omp" 
#test "pairhmm-3-hybrid-float-double" 
#test "pairhmm-4-hybrid-diagonal" 
#test "pairhmm-5-hybrid-diagonal-homogeneus" 
#test "pairhmm-6-onlythreediags" 
test "pairhmm-7-presse" 
#test "pairhmm-8-sse"
#test "pairhmm-dev"

