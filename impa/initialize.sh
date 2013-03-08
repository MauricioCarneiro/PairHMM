#!/bin/bash

tools=./tools

go(){
	for file in $1/*.in.txt
	do
		fn=$(basename "$file")
		fn="${fn%.in.txt}"

		echo "Generating "$fn.in.cuda" and "$fn.out.reference
		$tools/curr2cuda $1/$fn.in.txt > $1/$fn.in.cuda
		$tools/extract_answers $1/$fn.in.txt > $1/$fn.out.reference
	done
}

clean_mac_files(){
	find $1 -name "._*" -delete
	find $1 -name ".DS_Store" -delete
}

echo " "
echo "======================================================================"
echo "Initialize" 
echo " "
echo "This script creates the inputs formatted for C and CUDA versions,"
echo "and the reference outputs"
echo " "
echo "Usage: "$(basename $0)" <path to directory containing the input>*"
echo " "
echo "IMPORTANT NOTE: The directory MUST contain ONLY the inputs, formatted"
echo "                for the current java version. The files MUST be named"
echo "                {something}.in.txt."
echo " "
echo "======================================================================"
echo " "

if [ "$1" = "" ]; then
	echo "WRONG SYNTAX. SEE THE MESSAGE ABOVE"
	exit 0
fi

clean_mac_files $1
go $1

