all:	matcmp statistics curr2cuda extract_answers timeformatter

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

clean:
	rm -f matcmp statistics curr2cuda extract_answers timeformatter

.FORCE:

