#!/bin/bash

testcase(){
	tc="medium"

	suffix=$blockwidth"x"$blockheight"-"$comprows"-"$compdiags"-"$compbuffsize
	fileinp="../../../../test_data/"$tc"Test.in.cuda"
	fileout="$tc.out.mine-"$suffix
	filestdout="$tc.stdout.mine-"$suffix
	filetime="$tc.out.time-"$suffix
	outref="../../../../test_data/"$tc"Test.out.reference"
	binfile="pairhmm-cuda-"$suffix

	echo -n "    compiling..."

	CUDA_VERSION=5.0
	CUDA=/usr/local/cuda-$CUDA_VERSION
	LIBS=$CUDA/lib64
	INCS=$CUDA/include
	NVCC=$CUDA/bin/nvcc
	$NVCC --ptxas-options=-v -L$LIBS -I$INCS pairhmm.cu common.cu -O2 -arch=sm_20 -DBLOCKWIDTH=$1 -DBLOCKHEIGHT=$2 -DCOMPROWS=$3 -DCOMPDIAGS=$4 -DCOMPBUFFSIZE=$5 -o $binfile -lm > $binfile.compilerstdout 2> $binfile.compilerstderr

	echo "done"
	echo -n "    executing..."

	(/usr/bin/time -f "%e" ./$binfile $fileinp $fileout > $filestdout) 2> $filetime

	echo "done"
	echo -n "    biggest error: "

	err=`../tools/matcmp $outref $fileout | grep "The biggest difference is: "`
	err=${err#"The biggest difference is: "}

	echo $err

	echo -n "    time: "
	ttt=`cat $filetime`
	tottime=`../tools/timeformatter s $ttt`
	echo $tottime
}

for blockwidth in 32 64 96; do
comprows=$blockwidth
for blockheight in 1 3; do
	compdiags=30
	compbuffsize=30
	echo "Testing the case BLOCKSIZE="$blockwidth"x"$blockheight", COMPROWS="$comprows", COMPDIAGS="$compdiags", COMPBUFFSIZE="$compbuffsize
	testcase $blockwidth $blockheight $comprows $compdiags $compbuffsize
done #blockheight
done #blockwidth

blockwidth=96
comprows=$blockwidth
blockheight=3
compdiags=60
compbuffsize=30
echo "Testing the case BLOCKSIZE="$blockwidth"x"$blockheight", COMPROWS="$comprows", COMPDIAGS="$compdiags", COMPBUFFSIZE="$compbuffsize
testcase $blockwidth $blockheight $comprows $compdiags $compbuffsize

