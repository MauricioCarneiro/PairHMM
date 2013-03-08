#!/bin/bash

hash=`git rev-parse HEAD | head -c 15`

testjava(){
	basedir=${PWD}
	odir="benchmark."$hash
	mkdir -p $odir

	jarfile=../build/libs/PairHMM-0.1.jar

	echo " "
	echo "===================================="
	echo "Running java tests..."
	echo " "
	for file in $1/*.in.txt
	do
		fn=$(basename "$file")
		fn="${fn%.in.txt}"
		echo "Testcase: =="$fn"=="

		#loglesscaching
		echo -n "    loglesscaching: "
		out=$odir/$fn.java.loglesscaching.out
		stdout=$odir/$fn.java.loglesscaching.stdout
		time=$odir/$fn.java.loglesscaching.time
		(/usr/bin/time -f "%e" java -Xmx4g -jar $jarfile $file $out --impa-mode > $stdout) 2> $time
		echo "done"

		#caching
		echo -n "    caching: "
		out=$odir/$fn.java.caching.out
		stdout=$odir/$fn.java.caching.stdout
		time=$odir/$fn.java.caching.time
		(/usr/bin/time -f "%e" java -Xmx4g -jar $jarfile $file $out --caching --impa-mode > $stdout) 2> $time
		echo "done"

		#original
		echo -n "    original: "
		out=$odir/$fn.java.original.out
		stdout=$odir/$fn.java.original.stdout
		time=$odir/$fn.java.original.time
		(/usr/bin/time -f "%e" java -Xmx4g -jar $jarfile $file $out --original --impa-mode > $stdout) 2> $time
		echo "done"

		#exact
		echo -n "    exact: "
		out=$odir/$fn.java.exact.out
		stdout=$odir/$fn.java.exact.stdout
		time=$odir/$fn.java.exact.time
		(/usr/bin/time -f "%e" java -Xmx4g -jar $jarfile $file $out --exact --impa-mode > $stdout) 2> $time
		echo "done"
	done
}

testc(){
	basedir=${PWD}
	odir="benchmark."$hash
	mkdir -p $odir

	echo " "
	echo "===================================="
	echo "Running C tests..."
	echo " "
	for file in $1/*.in.cuda
	do
		fn=$(basename "$file")
		fn="${fn%.in.cuda}"
		echo -n "Testcase: =="$fn"=="
		out=$odir/$fn.c.out
		stdout=$odir/$fn.c.stdout
		time=$odir/$fn.c.time
		(/usr/bin/time -f "%e" ./src/pairhmm-c $file $out > $stdout) 2> $time
		echo " done. "
		echo "   generated:"
		echo "      "$out
		echo "      "$stdout
		echo "      "$time
	done
}

testcuda(){
	basedir=${PWD}
	odir="benchmark."$hash
	mkdir -p $odir

	echo " "
	echo "===================================="
	echo "Running CUDA tests..."
	echo " "
	for file in $1/*.in.cuda
	do
		fn=$(basename "$file")
		fn="${fn%.in.cuda}"
		echo -n "Testcase: =="$fn"=="
		out=$odir/$fn.cuda.out
		stdout=$odir/$fn.cuda.stdout
		time=$odir/$fn.cuda.time
		(/usr/bin/time -f "%e" ./src/pairhmm-cuda $file $out > $stdout) 2> $time
		echo " done. "
		echo "   generated:"
		echo "      "$out
		echo "      "$stdout
		echo "      "$time
	done
}

echo " "
echo "======================================================================"
echo "Test" 
echo " "
echo "This script tests the different implementations and produces the files"
echo "needed for benchmarking"
echo " "
echo "Usage: "$(basename $0)" <path to initialized dir> <java | c | cuda>*"
echo " "
echo "======================================================================"
echo " "

if [ "$2" = "" ]; then
	echo "WRONG SYNTAX. SEE THE MESSAGE ABOVE"
	exit 0
fi

for v in "$@"
do
	case $v in
		"java") testjava $1 ;;
		"c") testc $1 ;;
		"cuda") testcuda $1 ;;
	esac
done

