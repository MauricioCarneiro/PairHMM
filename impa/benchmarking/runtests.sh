#!/bin/bash

showhelp()
{
	echo ""
	echo "Usage: $(basename $ME) <ROOT_DIR>"
	echo ""
	echo "Example: $(basename $ME) ../../../../PairHMM/"
	echo ""
	echo "This program takes as parameter one directory path ROOT_DIR, and"
	echo "assumes the following directory structure:"
	echo ""
	echo "   ROOT_DIR/PairHMM ===> (complete repository)"
	echo "   ROOT_DIR/PairHMM/impa"
	echo "   ROOT_DIR/PairHMM/impa/benchmarking/results"
	echo "   ROOT_DIR/PairHMM/impa/tools"
	echo "   ROOT_DIR/PairHMM/impa/src"
	echo "   ROOT_DIR/files ===> (.tar.bz2 files)"
	echo "   ROOT_DIR/test_data ===> (input files, and reference outputs)"
	echo ""
	echo "This script do the following steps:"
	echo ""
	echo "1) clean_and_make"
	echo "   'make clean; make' inside ROOT_DIR/PairHMM/impa/src and ROOT_DIR/PairHMM/impa/tools"
	echo ""
	echo "2) replicate_directory_structure"
	echo "   Replicates ROOT_DIR in HOST:~/ROOT_DIR, for HOST in {dijkstra, bresenham, germain}"
	echo ""
	echo "3) run_on_hosts_and_copy_back"
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

clean_and_make()
{
	echo "STEP 1: clean_and_make"
	cd "$ROOT_DIR/PairHMM/impa/src"
	echo "    cleaning impa/src"
	make clean 1>> "$INITIAL_DIR/make.log" 2>> "$INITIAL_DIR/make.log"
	echo "    making impa/src"
	make 1>> "$INITIAL_DIR/make.log" 2>> "$INITIAL_DIR/make.log"
	cd "$ROOT_DIR/PairHMM/impa/tools"
	echo "    cleaning impa/tools"
	make clean 1>> "$INITIAL_DIR/make.log" 2>> "$INITIAL_DIR/make.log"
	echo "    making impa/tools"
	make 1>> "$INITIAL_DIR/make.log" 2>> "$INITIAL_DIR/make.log"
}

replicate_directory_structure()
{
	echo "STEP 2: replicate_directory_structure"
	for host in dijkstra bresenham germain 
	do
		echo "    replicating directory in "$host
		rsync -avz --delete $ROOT_DIR $host:~ 1>> "$INITIAL_DIR/rsync.log" 2>> "$INITIAL_DIR/rsync.log"
	done
}

run_on_hosts_and_copy_back()
{
	echo "STEP 3: run_on_hosts_and_copy_back"
	for host in "dijkstra" "bresenham" "germain" 
	do
		echo "    launching tests at "$host
		ssh eric@$host "~/"$(basename $ROOT_DIR)"/PairHMM/impa/benchmarking/"$ME" --run_test_and_copy_to_viete ~/"$(basename $ROOT_DIR)" ~/projects/PairHMM" &
	done
}

run_test_and_copy_to_viete()
{
	cd $ROOT_DIR/PairHMM/impa/src
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
			scp $fn eric@viete:$DESTINATION_DIR/PairHMM/impa/benchmarking/results
			rm $fn
		done
	done
}


ME=$0
INITIAL_DIR=`pwd`
INITIAL_DIR=${INITIAL_DIR%/}
ROOT_DIR=$1
ROOT_DIR=${ROOT_DIR%/}

if [ "$1" = "--run_test_and_copy_to_viete" ]; then
	ROOT_DIR=$2
	ROOT_DIR=${ROOT_DIR%/}
	DESTINATION_DIR=$3
	DESTINATION_DIR=${DESTINATION_DIR%/}
	run_test_and_copy_to_viete
else
	rm -f *.log
	check_args
	clean_and_make
	replicate_directory_structure
	run_on_hosts_and_copy_back
fi

