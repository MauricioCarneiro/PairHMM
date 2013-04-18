#!/bin/bash

if [ "$1" = "" ]; then
	echo ""
	echo "Usage: $(basename $0) <WIDTHxHEIGHT>"
	echo "Example: $(basename $0) 32x3"
	echo ""
	exit 0
fi

go()
{
	ssh eric@$2 'rm -rf /home/eric/PairHMM/rep/PairHMM/impa'
	scp -r ../../impa eric@$2:/home/eric/PairHMM/rep/PairHMM
	ssh eric@$2 'cd /home/eric/PairHMM/rep/PairHMM/impa/src; make cleanlocal; make cuda; make inp'
	ssh eric@$2 'cd /home/eric/PairHMM/rep/PairHMM/impa/tools; make clean; make; '
	ssh eric@$2 "cd /home/eric/PairHMM/rep/PairHMM/impa/src; ./benchmarks.sh $1"
}

for v in dijkstra germain bresenham
do
	go $1 $v &
done

