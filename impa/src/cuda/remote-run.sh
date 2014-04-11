#!/bin/bash

launch_remote()
{
	Host=$1
	for testcase in "${testcases[@]}"
	do
		ssh $Host rm -rf pairhmm-dev-* matcmp reports/ firstnotacceptable
		ssh $Host mkdir ~/reports
		scp pairhmm-dev-* $Host:~/
		scp ../../tools/matcmp $Host:~/
		scp ../../tools/firstnotacceptable $Host:~/

		inp="test_data/"$testcase"Test.in"
		ref="test_data/"$testcase"Test.out"
		for threads in "32" "64" "96";
		#for threads in "96";
		do
			o="reports/"$testcase"-$threads-$Host.out"
			e="reports/"$testcase"-$threads-$Host.err"
			p="reports/"$testcase"-$threads-$Host.precision"
			ssh $Host "~/pairhmm-dev-$threads < $inp > $o 2> $e"
			ssh $Host "~/matcmp $ref $o > $p"
		done
		scp $Host:~/reports/*.err .
		scp $Host:~/reports/*.precision .
	done
}

launch_local()
{
	for testcase in "${testcases[@]}"
	do
		inp="../../../../test_data/"$testcase"Test.in"
		ref="../../../../test_data/"$testcase"Test.out"
		for threads in "32" "64" "96";
		do
			o="$testcase-$threads-$Host.out"
			e="$testcase-$threads-$Host.err"
			p="$testcase-$threads-$Host.precision"
			./pairhmm-dev-$threads < $inp > $o 2> $e
			../../tools/matcmp $ref $o > $p
		done
	done
}

if [ $# -ne 2 ] 
then
 	echo " "
	echo "Usage example:"
	echo "    $0 \"germain bresenham\" \"small medium\""
	echo " "
	exit 0
fi


if [ ! -f pairhmm-dev-32 ]; then
	make
fi

IFS=' ' read -a hosts <<< "$1"
IFS=' ' read -a testcases <<< "$2"

for Host in "${hosts[@]}"
do
	if [ $Host != "viete" ]
	then
		launch_remote $Host 
	else
		launch_local
	fi
done

