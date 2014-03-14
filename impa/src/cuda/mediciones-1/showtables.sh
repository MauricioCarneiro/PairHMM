#!/bin/bash

IFS=' ' read -a hosts <<< "$1"
IFS=' ' read -a testcases <<< "$2"

for testcase in "${testcases[@]}"
do
	printf "\nCASE: $testcase\n\n"
	printf "%15s %15s %15s %15s\n" "PRECISION" "32" "64" "96"
	for Host in "${hosts[@]}"
	do
		printf "%15s " $Host
		for threads in "32" "64" "96";
		do
			precfile=$testcase"-$threads-$Host.precision"
			pr=""
			if [ -f $precfile ]
			then
				pr=`cat $precfile | grep "The biggest difference is: "`
				pr=${pr#"The biggest difference is: "}
				printf "%15s " $pr
			fi
		done
		printf "\n";
	done
	printf "\n";
	printf "%15s %15s %15s %15s\n" "TIME" "32" "64" "96"
	for Host in "${hosts[@]}"
	do
		printf "%15s " $Host
		for threads in "32" "64" "96";
		do
			timefile=$testcase"-$threads-$Host.err"
			ti=""
			if [ -f $timefile ]
			then
				ti=`cat $timefile | grep "COMPUTATION: "`
				ti=${ti#"COMPUTATION: "}
				ti=${ti%% "seconds"}
			fi
			printf "%15s " $ti
		done
		printf "\n";
	done
done
