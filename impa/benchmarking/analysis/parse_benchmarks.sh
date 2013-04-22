#!/bin/bash

timeformatter="../tools/timeformatter ms "

go()
{
	fn=$1

	fn1=`echo $fn | tr - .`
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

	echo -n $year";"$month";"$day";"$hour";"$minute";"$language";"$test";"$host";"$size";"$commitid";"
	echo $init_memory";"$mallocs";"$computation";"$kernel";"$output";"$total_program";"$total_time";"$biggest_error
}

#echo -n "YEAR"";""MONTH"";""DAY"";""HOUR"";""MINUTE"";""LANGUAGE"";""TEST"";""HOST"";""SIZE"";"
#echo -n "COMMIT_ID"";""INIT_MEMORY"";""MALLOCS"";""COMPUTATION"";""KERNEL"";""OUTPUT"";"
#echo "TOT_(measured_by_program)"";""TOT_(measured_by_time)"";""BIGGEST_ERROR"

for fn in *.50a27502a983f93
do
	if [ -f $fn ]; then
		go $fn
	fi
done

for fn in *.darboux
do
	if [ -f $fn ]; then
		go $fn
	fi
done

