#!/bin/tcsh

if ($#argv != 1) then
	echo -n "you MUST specify the target name"
	exit 1
endif

javac RunTest.java
mkdir $1
mv *.class $1/
cd $1
java RunTest ../batch1.in && diff output.txt ../batch1.out
