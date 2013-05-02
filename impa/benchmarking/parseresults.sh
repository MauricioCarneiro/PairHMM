#!/bin/bash

timeformatter="../tools/timeformatter ms "

first_line()
{
	echo -n "YEAR;MONTH;DAY;HOUR;MINUTE;LANGUAGE;TEST;HOST;GPU_MODEL;SIZE;"
	echo -n "COMMIT_ID;INIT_MEMORY;MALLOCS;COMPUTATION;KERNEL;OUTPUT;"
	echo "TOT_(measured_by_program);TOT_(measured_by_time);BIGGEST_ERROR"
}

all_lines()
{
	fn=$1

	fn1=$(basename $fn)
	fn1=`echo $fn1 | tr - .`
	fn1=`echo $fn1 | tr \: .`
	IFS='.' read -ra arr <<< "$fn1"

	year=${arr[0]}
	month=${arr[1]}
	day=${arr[2]}
	hour=${arr[3]}
	minute=${arr[4]}
	language=${arr[5]}
	test=${arr[6]}
	host=${arr[7]}
	gpu_model="n/a"
	case $host in
		"dijkstra")
			gpu_model="GeForce GTX 670"
			;;
		"germain")
			gpu_model="GeForce GTX 480"
			;;
		"bresenham")
			gpu_model="GeForce GTX 680"
			;;
		"viete")
			gpu_model="GeForce GTX 470"
			;;
		*)
			;;
	esac
	size=${arr[8]}
	commitid=${arr[9]}
	init_memory=`cat $fn | grep INIT_MEMORY`
	init_memory=${init_memory#"INIT_MEMORY: "}
	init_memory=`$timeformatter $init_memory`
	mallocs=`cat $fn | grep MALLOCS`
	mallocs=${mallocs#"MALLOCS: "}
	mallocs=`$timeformatter $mallocs`
	computation=`cat $fn | grep COMPUTATION`
	computation=${computation#"COMPUTATION: "}
	computation=`$timeformatter $computation`
	kernel=`cat $fn | grep KERNEL`
	kernel=${kernel#"KERNEL: "}
	kernel=`$timeformatter $kernel`
	output=`cat $fn | grep OUTPUT`
	output=${output#"OUTPUT: "}
	output=`$timeformatter $output`
	total_program=`cat $fn | grep "TOTAL (measured inside program):"`
	total_program=${total_program#"TOTAL (measured inside program): "}
	total_program=`$timeformatter $total_program`
	total_time=`cat $fn | grep "TOTAL (measured with 'time'):"`
	total_time=${total_time#"TOTAL (measured with 'time'): "}
	total_time=`$timeformatter $total_time`
	biggest_error=`cat $fn | grep "BIGGEST ERROR:"`
	biggest_error=${biggest_error#"BIGGEST ERROR: "}

	echo -n "$year;$month;$day;$hour;$minute;$language;$test;$host;$gpu_model;"
	echo -n "$size;$commitid;$init_memory;$mallocs;$computation;$kernel;"
	echo "$output;$total_program;$total_time;$biggest_error"
}

if [ "$1" = "" ]; then
	echo ""
	echo "Usage: "$(basename $0)" <expression describing the paths with the results>"
	echo ""
	echo "Example: "$(basename $0)" ./results/*.4ed8ee5dc1ed885"
	echo ""
	exit 0
fi

first_line
for f in $@
do
	all_lines $f
done

