#ifndef INPUT_H
#define INPUT_H

#define MAX_TESTCASES_BUNCH_SIZE 100

#define VECTOR_SIZE 8

struct testcase
{
	int rslen, haplen, *q, *i, *d, *c;
	char *hap, *rs;
};

int read_testcase(testcase *);
int read_a_bunch_of_testcases(testcase *, int);

#endif

