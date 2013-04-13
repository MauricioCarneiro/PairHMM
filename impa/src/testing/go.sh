#!/bin/bash

testsmall(){
	fileinp="../../../../test_data/smallTest.in.cuda"
	fileout="small.out.mine"
	outref="../../../../test_data/smallTest.out.reference"
	#t=`/usr/bin/time -f "%e" ./pairhmm-cuda-cuda $fileinp $fileout`
	./pairhmm-cuda-cuda $fileinp $fileout
	err=`../tools/matcmp $outref $fileout | grep "The biggest difference is: "`
	err=${err#"The biggest difference is: "}
	if [ $err != "0.008955" ]; then
		err=`echo $err"<================================"`
	fi
	#echo $t $err
	echo "#"$1":" $err
}

testmedium(){
	fileinp="../../../../test_data/mediumTest.in.cuda"
	fileout="medium.out.mine"
	outref="../../../../test_data/mediumTest.out.reference"
	#t=`/usr/bin/time -f "%e" ./pairhmm-cuda-cuda $fileinp $fileout`
	./pairhmm-cuda-cuda $fileinp $fileout
	err=`../tools/matcmp $outref $fileout | grep "The biggest difference is: "`
	err=${err#"The biggest difference is: "}
	if [ $err != "2.344392" ]; then
		err=`echo $err"<================================"`
	fi
	#echo $t $err
	echo "#"$1":" $err
}

testlarge(){
	fileinp="../../../../test_data/largeTest.in.cuda"
	fileout="large.out.mine"
	outref="../../../../test_data/largeTest.out.reference"
	#t=`/usr/bin/time -f "%e" ./pairhmm-cuda-cuda $fileinp $fileout`
	./pairhmm-cuda-cuda $fileinp $fileout
	err=`../tools/matcmp $outref $fileout | grep "The biggest difference is: "`
	err=${err#"The biggest difference is: "}
	if [ $err != "2.378401" ]; then
		err=`echo $err"<================================"`
	fi
	#echo $t $err
	echo "#"$1":" $err
}

if [ "$1" = "" ]; then
	echo " "
	echo "======================================================================"
	echo "Test" 
	echo " "
	echo "This script run the program and compares against the expected output"
	echo " "
	echo "Usage: "$(basename $0)" {small | medium | large}"
	echo " "
	echo "======================================================================"
	echo " "
	exit 0
fi

for v in "$@"
do
	case $v in
		"small")
			caso=1
			while [ 1 ]
			do 
				testsmall $caso
				caso=$((caso+1))
			done ;;
		"medium")
			caso=1
			while [ 1 ]
			do 
				testmedium $caso
				caso=$((caso+1))
			done ;;
		"large") 
			caso=1
			while [ 1 ]
			do 
				testlarge $caso
				caso=$((caso+1))
			done ;;
	esac
done

