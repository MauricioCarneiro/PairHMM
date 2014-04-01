#ifndef INPUT_H
#define INPUT_H

#define MAX_TESTCASES_BUNCH_SIZE 100

#define VECTOR_SIZE 8

struct testcase
{
	int rslen, haplen, *q, *i, *d, *c;
	char *hap, *hap_data, *rs;

  testcase() : q(NULL), i(NULL), d(NULL), c(NULL), hap(NULL), hap_data(NULL), rs(NULL)
  {
	}

  void free()
  {
	  delete[] hap_data;
	  delete[] q;
	  delete[] i;
	  delete[] d;
	  delete[] c;
		delete[] rs;

    hap_data = hap = rs = NULL;
    q = i = d = c = NULL;
    rslen = haplen = 0;
  }

	~testcase()
	{
      free();
	}
};

int read_testcase(testcase *);
int read_a_bunch_of_testcases(testcase *, int);

#endif

