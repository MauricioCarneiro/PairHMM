all:	matcmp statistics curr2cuda extract_answers timeformatter firstnotacceptable adaptcsv

adaptcsv:	.FORCE
	g++ adaptcsv.cc -O2 -W -Wall -fopenmp -pedantic -o adaptcsv -lgomp -lm

matcmp:	.FORCE
	g++ matcmp.cc -O2 -W -Wall -fopenmp -pedantic -o matcmp -lgomp -lm

statistics:	.FORCE
	gcc statistics.c common.c -O2 -W -Wall -fopenmp -pedantic -o statistics -lgomp -lm

curr2cuda:	.FORCE
	g++ curr2cuda.cc -O2 -W -Wall -fopenmp -pedantic -o curr2cuda -lgomp -lm

extract_answers:	.FORCE
	gcc extract_answers.c -O2 -W -Wall -fopenmp -pedantic -o extract_answers -lgomp -lm

timeformatter:	.FORCE
	gcc timeformatter.c -O2 -W -Wall -fopenmp -pedantic -o timeformatter -lgomp -lm

firstnotacceptable:	.FORCE
	g++ firstnotacceptable.cc -O2 -W -Wall -fopenmp -pedantic -o firstnotacceptable -lgomp -lm

maxlen: .FORCE
	g++ max-len.cc -O2 -o max-len

clean:
	rm -f max-len matcmp statistics curr2cuda extract_answers timeformatter firstnotacceptable adaptcsv

.FORCE:

