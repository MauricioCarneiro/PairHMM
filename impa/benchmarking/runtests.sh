#!/bin/bash

showhelp()
{
	echo ""
	echo "Usage: ./runtests <ROOT_DIR>"
	echo ""
	echo "Example: ./runtests ../../../../PairHMM/"
	echo ""
	echo "This program assumes the following directory structure:"
	echo ""
	echo "   ROOT_DIR/files ===> (.tar.bz2 files)"
	echo "   ROOT_DIR/test_data"
	echo "   ROOT_DIR/PairHMM"
	echo "   ROOT_DIR/PairHMM/impa"
	echo "   ROOT_DIR/PairHMM/impa/benchmarking/results"
	echo "   ROOT_DIR/PairHMM/impa/tools"
	echo "   ROOT_DIR/PairHMM/impa/src/c++"
	echo "   ROOT_DIR/PairHMM/impa/src/cuda"
	echo ""
	echo "This script assumes that tools binaries and cuda binaries exists. It"
	echo "runs them and puts the results in ROOT_DIR/PairHMM/impa/benchmarking/results"
	echo ""
}

check_args()
{
	if [ "$ROOT_DIR" = "" ] || [ "$ROOT_DIR" = "help" ] || [ "$ROOT_DIR" = "-help" ]; then
		showhelp
		exit 0;
	fi

	if [ ! -d "$ROOT_DIR" ] || \
		[ ! -d "$ROOT_DIR/PairHMM" ] || \
		[ ! -d "$ROOT_DIR/PairHMM/impa" ] || \
		[ ! -d "$ROOT_DIR/PairHMM/impa/src" ] || \
		[ ! -d "$ROOT_DIR/PairHMM/impa/tools" ] || \
		[ ! -d "$ROOT_DIR/PairHMM/impa/benchmarking/results" ] || \
		[ ! -d "$ROOT_DIR/files" ] || \
		[ ! -d "$ROOT_DIR/test_data" ]; then
		echo "Wrong directory structure inside *** "$ROOT_DIR" *** See the help using "$(basename $ME)" help"
		exit 0;
	fi
}

run()
{
	cd $ROOT_DIR/PairHMM/impa/src/cuda
	hash=`git rev-parse HEAD | head -c 15`

	for size in "32x1" "32x3" "64x1" "64x3" "96x1" "96x3"
	do
		for file in $ROOT_DIR/test_data/*.in.cuda
		do
			fn=$(basename "$file")
			fn=${fn%.in.cuda}
			fileref="${file%.in.cuda}.out"
			fn=`date +"%Y.%m.%d-%H.%M"`"-cuda.$fn."`hostname`".$size.$hash"
			rm -f fout fstdout fstderr ftime
			(/usr/bin/time -f "%e" ./pairhmm$size $file fout > fstdout) 2> ftime
			o=`cat fstdout`
			t=`cat ftime`
			t=`echo "$t*1000" | bc`
			t="TOTAL (measured with 'time'): $t"
			e=`$ROOT_DIR/PairHMM/impa/tools/matcmp $fileref fout | grep "The biggest difference is: "`
			e=${e#"The biggest difference is: "}

			rm -f fout fstdout fstderr ftime
			echo "$o" > $fn
			echo "$t" >> $fn
			echo "BIGGEST ERROR: $e" >> $fn
			mv $fn $ROOT_DIR/PairHMM/impa/benchmarking/results/
		done
	done
}

ROOT_DIR=$1
rm -f *.log
check_args
run

