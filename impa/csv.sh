#!/bin/bash

hash=`git rev-parse HEAD | head -c 15`

go(){
	odir="benchmark."$hash
	ofile=./sheets/sheet.$hash.csv

	echo -n ",,java.original," > $ofile
	echo -n ",,java.exact," >> $ofile
	echo -n ",,java.caching," >> $ofile
	echo -n ",,java.loglesscaching," >> $ofile
	echo -n ",,C," >> $ofile
	echo ",,CUDA," >> $ofile
	echo -n ",Tot Time, Comp. Time, Err" >> $ofile
	echo -n ",Tot Time, Comp. Time, Err" >> $ofile
	echo -n ",Tot Time, Comp. Time, Err" >> $ofile
	echo -n ",Tot Time, Comp. Time, Err" >> $ofile
	echo -n ",Tot Time, Comp. Time, Err" >> $ofile
	echo ",Tot Time, Comp. Time, Err" >> $ofile

	for file in $1/*.in.cuda
	do
		fn=$(basename "$file")
		fn="${fn%.in.cuda}"

		echo -n $fn >> $ofile
		for suffix in "java.original" "java.exact" "java.caching" "java.loglesscaching" "c" "cuda"; do
			if [ -f "$odir/$fn.$suffix.time" ]; then
				tottime=`head -1 $odir/$fn.$suffix.time`
			else
				tottime=""
			fi
			if [ -f "$odir/$fn.$suffix.stdout" ]; then 
				comptime=`cat $odir/$fn.$suffix.stdout | grep COMPUTATION_TIME`
				comptime=${comptime#"COMPUTATION_TIME"}
			else
				comptime=""
			fi
			if [ -f "$odir/$fn.$suffix.out" ]; then 
				err=`./tools/matcmp $1/$fn.out.reference $odir/$fn.$suffix.out | grep "The biggest difference is: "`
				err=${err#"The biggest difference is: "}
				#err=""
			else
				err=""
			fi

			tottime=`./tools/timeformatter s $tottime`
			comptime=`./tools/timeformatter ms $comptime`

			echo -n ,$tottime,$comptime,$err >> $ofile
		done
		echo "" >> $ofile
	done
}

usage(){
	echo " "
	echo "======================================================================"
	echo "Csv" 
	echo " "
	echo "This script generates a comma separated values sheet with the obtained"
	echo "measurements"
	echo " "
	echo "Usage: "$(basename $0)" <path to initialized dir>"
	echo " "
	echo "======================================================================"
	echo " "
}

if [ "$1" = "" ]; then
	usage
	exit 0
fi

go $1

