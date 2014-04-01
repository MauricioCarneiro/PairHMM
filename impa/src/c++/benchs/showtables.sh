#!/bin/bash

IFS=' ' read -a hosts <<< "$1"
IFS=' ' read -a testcases <<< "$2"

for testcase in "${testcases[@]}"
do
	printf "\nCASE: $testcase\n\n"

		printf "PRECISION"
		for str in "pairhmm-1-base" "pairhmm-2-omp" ; # "pairhmm-3-hybrid-float-double" "pairhmm-4-hybrid-diagonal" "pairhmm-5-hybrid-diagonal-homogeneus" "pairhmm-6-onlythreediags" ;
		do
			printf ",%s" "$str"
		done
		printf "\n";

	for Host in "${hosts[@]}"
	do
		printf "%s" $Host
		for bin in "pairhmm-1-base" "pairhmm-2-omp" ; # "pairhmm-3-hybrid-float-double" "pairhmm-4-hybrid-diagonal" "pairhmm-5-hybrid-diagonal-homogeneus" "pairhmm-6-onlythreediags" ;
		do
			precfile="$testcase-$bin-$Host.precision"
			pr=""
			if [ -f $precfile ]
			then
				pr=`cat $precfile | grep "The biggest difference is: "`
				pr=${pr#"The biggest difference is: "}
				printf ",%s" $pr
			fi
		done
		printf "\n";
	done
	printf "\n";

	printf "TIME"
	for str in "pairhmm-1-base" "pairhmm-2-omp" ; #"pairhmm-3-hybrid-float-double" "pairhmm-4-hybrid-diagonal" "pairhmm-5-hybrid-diagonal-homogeneus" "pairhmm-6-onlythreediags" ;
	do
		printf ",%s" "$str"
	done
	printf "\n";

	for Host in "${hosts[@]}"
	do
		printf "%s" $Host
		for bin in "pairhmm-1-base" "pairhmm-2-omp" ; # "pairhmm-3-hybrid-float-double" "pairhmm-4-hybrid-diagonal" "pairhmm-5-hybrid-diagonal-homogeneus" "pairhmm-6-onlythreediags" ;
		do
			timefile="$testcase-$bin-$Host.err"
			ti=""
			if [ -f $timefile ]
			then
				ti=`cat $timefile | grep "COMPUTATION: "`
				ti=${ti#"COMPUTATION: "}
				ti=${ti%% "seconds"}
			fi
			printf ",%s" $ti
		done
		printf "\n";
	done
done
